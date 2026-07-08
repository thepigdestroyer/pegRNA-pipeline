"""Recoded TJ-pegRNA workflow: enumeration, gates, and the ranked-pool-walk.

This module wires together the pure building blocks (microhomology, recode_arm,
nick_search, INV-0) and the injected gate scorers (splice, importance, full
structure, complex docking) into the candidate enumeration and the
ranked-pool-walk orchestrator. Heavy scorers are passed in (fakes in CI, real
adapters on CX3); this module imports none of them at module load.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from pegrna_pipeline.hpc.allele_specificity import find_allele_specific
from pegrna_pipeline.hpc.codon_expander import _gc_content, enumerate_variants
from pegrna_pipeline.hpc.complex_dock import passes_docking_gate
from pegrna_pipeline.hpc.cryptic_splice import passes_cryptic_splice_gate
from pegrna_pipeline.hpc.full_structure import functional_region_pairing
from pegrna_pipeline.hpc.gate_thresholds import THRESHOLDS
from pegrna_pipeline.hpc.nick_search import enumerate_nick_sites, screen_nick_candidates
from pegrna_pipeline.hpc.recode_arm import recode_arm
from pegrna_pipeline.hpc.scarless import assert_inv0
from pegrna_pipeline.hpc.splice_screen import passes_splice_gate
from pegrna_pipeline.hpc.strand_orientation import orient_rcpbs2, orient_rtt_tail
from pegrna_pipeline.hpc.structure_importance import span_blocked
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD, assemble_tj_pegrna_v2
from pegrna_pipeline.hpc.tprt_sim import TPRTMismatch, simulate_tprt

__all__ = [
    "PbsRttVariant",
    "enumerate_pbs_rtt_variants",
    "passes_rcpbs2_gate",
    "run_recoded_tjpe",
    "RecodedTJPEResult",
    "NoValidPegRNA",
]


class NoValidPegRNA(RuntimeError):
    """Raised (fail-loud) when no candidate survives all gates. Carries a
    per-candidate diagnostic of which gate killed each one."""


@dataclass(frozen=True)
class PbsRttVariant:
    spacer: str
    rtt_dna: str
    pbs_length: int


def enumerate_pbs_rtt_variants(
    *,
    spacer: str,
    rtt_aa: str,
    pbs_lengths: Tuple[int, ...],
    max_rtt_variants: int = 100,
    ranking: str = "cai",
) -> List[PbsRttVariant]:
    """Enumerate (synonymous RTT) × (PBS1 length); the spacer is FIXED by site S.

    These are the candidates handed to PRIDICT2 (reverse-targeting) for ranking.
    The spacer is never varied — there is one (or a few) protospacer(s) that nick
    at S, and that placement is fixed (decision 1). Deterministic order:
    RTT-variant rank, then PBS length.

    ``ranking`` defaults to ``"cai"`` (ADR 0009): the synonymous RTT variants for the
    recoded arm are ranked by human Codon Adaptation Index so the most
    human-expressed recoding is preferred for PRIDICT2 + assembly. Pass
    ``ranking="gc_c"`` for the legacy synthesis-stability order.
    """
    rtts = enumerate_variants(rtt_aa, max_variants=max_rtt_variants, ranking=ranking)
    out: List[PbsRttVariant] = []
    for rv in rtts:
        for length in pbs_lengths:
            out.append(PbsRttVariant(spacer=spacer, rtt_dna=rv.dna, pbs_length=length))
    return out


def passes_rcpbs2_gate(
    rcpbs2: str,
    *,
    repetitive: bool,
    min_len: int = 11,
    max_len: int = 17,
    gc_lo: float = 0.40,
    gc_hi: float = 0.70,
    tm_lo: float = 37.0,
    tm_hi: float = 42.0,
) -> Tuple[bool, str]:
    """Joint rcPBS2 validation (decision 8): length, GC 40-70 %, Tm 37-42 C, not
    repetitive.

    This VALIDATES a chosen rcPBS2 (the RC-of-PBS2 element the tj-pegRNA carries).
    The length is selected upstream by :func:`_select_rcpbs2`, which scans the length
    range to hit the Tm band, so this gate should pass by construction for a selected
    rcPBS2 (it also guards a caller-supplied one). ``min_len``/``max_len`` default to
    the classic 11-17; :func:`_select_rcpbs2` uses the wider 10-18 scan from
    THRESHOLDS. Tm uses Biopython's nearest-neighbour model (imported lazily to keep
    the core import light).
    """
    if not (min_len <= len(rcpbs2) <= max_len):
        return False, f"length {len(rcpbs2)} outside [{min_len},{max_len}]"
    gc = _gc_content(rcpbs2)
    if not (gc_lo <= gc <= gc_hi):
        return False, f"GC {gc:.2f} outside [{gc_lo},{gc_hi}]"
    from Bio.SeqUtils import MeltingTemp as mt  # noqa: PLC0415

    tm = float(mt.Tm_NN(rcpbs2))
    if not (tm_lo <= tm <= tm_hi):
        return False, f"Tm {tm:.1f} outside [{tm_lo},{tm_hi}]"
    if repetitive:
        return False, "repetitive region"
    return True, "ok"


def _is_repetitive(seq: str) -> bool:
    """Cheap repetitiveness screen: a homopolymer run >= 5, or <= 2 distinct bases
    (e.g. pure AT/GC microsatellite-like), makes a poor, mis-priming PBS."""
    if len(set(seq)) <= 2:
        return True
    run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        if run >= 5:
            return True
    return False


def _select_rcpbs2(
    genome: str,
    *,
    nick_position: int,
    pegrna_strand: str,
    tm_lo: float,
    tm_hi: float,
    gc_lo: float,
    gc_hi: float,
    len_min: int,
    len_max: int,
) -> Tuple[Optional[str], str]:
    """Choose rcPBS2 by SCANNING its length to hit the Tm band, decoupled from PBS1.

    rcPBS2 is the second-strand primer-binding site carried by the tj-pegRNA:
    ``rcPBS2 = RC(TS)`` over the window immediately 5' of the ngRNA nick (ADR 0014).
    The target-strand window is the captured sense window for a reverse-targeting
    (``Rv``) pegRNA and its reverse complement for a forward (``Fw``) one; either way
    the emitted rcPBS2 is :func:`orient_rcpbs2` of that target-strand window. Its
    sequence is fixed by the genome, so the only free parameter is LENGTH. Fixing
    length (to PBS1's lengths) and then rejecting on Tm made the whole candidate pool
    collapse (all share ~one nick's window). Instead scan ``len_min..len_max`` and
    return the valid rcPBS2 whose Tm is closest to the band centre; ``(None, why)`` if
    none fit (that nick simply yields no priming-competent rcPBS2 — a per-nick miss,
    not a whole-run failure). Tm(NN)/GC/repetitiveness are reverse-complement
    invariant, so the length choice does not depend on the orientation.
    """
    from Bio.Seq import Seq  # noqa: PLC0415
    from Bio.SeqUtils import MeltingTemp as mt  # noqa: PLC0415

    center = (tm_lo + tm_hi) / 2.0
    best: Optional[Tuple[float, str, float]] = None  # (dist_to_center, rcpbs2, tm)
    considered = 0
    for length in range(len_min, len_max + 1):
        start = nick_position - length
        if start < 0:
            continue
        captured = genome[start:nick_position].upper()
        if len(captured) < length:
            continue
        if set(captured) - {"A", "C", "G", "T"}:
            # Non-ACGT window (e.g. an assembly-gap N-run) is a non-priming site;
            # Tm_NN silently ignores N, so reject the alphabet explicitly (per-nick miss).
            continue
        # Target-strand window at the capture locus, then rcPBS2 = RC(TS) (ADR 0014).
        ts_window = captured if pegrna_strand == "Rv" else str(Seq(captured).reverse_complement())
        rcpbs2 = orient_rcpbs2(ts_window)
        gc = _gc_content(rcpbs2)
        if not (gc_lo <= gc <= gc_hi):
            continue
        if _is_repetitive(rcpbs2):
            continue
        tm = float(mt.Tm_NN(rcpbs2))
        considered += 1
        if tm_lo <= tm <= tm_hi:
            d = abs(tm - center)
            if best is None or d < best[0]:
                best = (d, rcpbs2, tm)
    if best is None:
        return None, (
            f"no rcPBS2 length in [{len_min},{len_max}] hit Tm [{tm_lo},{tm_hi}] "
            f"with GC [{gc_lo},{gc_hi}] at nick {nick_position} "
            f"({considered} GC/repeat-valid lengths, none in Tm band)"
        )
    return best[1], f"Tm {best[2]:.1f} (len {len(best[1])})"


@dataclass
class RecodedTJPEResult:
    """Outcome of the recoded TJ-PE design for one target."""

    pegRNA_sequence: str
    spacer: str
    rcpbs2: str
    rtt: str
    pbs1: str
    polyu: str
    ngrna_protospacer: str
    site_s: int
    nick_position: int
    deletion_size: int
    recode_edits: int
    candidates_tried: int
    warnings: List[str] = field(default_factory=list)
    allele_specific_options: list = field(default_factory=list)
    #: The RTT decomposed into its labelled sub-parts, in RTT (5'->3') order, as a
    #: list of (label, sequence) with labels "insert" (the cargo) and
    #: "replacement_gdna" (the recoded deletion-replacement arm). Concatenating the
    #: sequences reproduces ``rtt`` exactly. Empty until set by the workflow.
    rtt_components: List[Tuple[str, str]] = field(default_factory=list)
    #: Convenience: the insert (cargo) sub-sequence as it appears in the RTT.
    insert_seq: str = ""
    #: Convenience: the replacement gDNA (recoded arm) sub-sequence in the RTT.
    replacement_gdna: str = ""
    #: Verification data (for the manual cross-check report). ``edited_allele`` is
    #: the predicted fully-edited genomic sequence; ``recode_span`` / ``locked_span``
    #: are genome-context coordinates of the synonymously-recoded and non-coding
    #: byte-exact re-templated regions. Empty/None until set by the workflow.
    edited_allele: str = ""
    recode_span: Optional[Tuple[int, int]] = None
    locked_span: Optional[Tuple[int, int]] = None
    #: Per-gate info string for the winning candidate, in pipeline order, as a list
    #: of (gate_id, info) — surfaced in the report + web flowchart so each step shows
    #: the score/verdict it produced. Empty until set by the workflow.
    gate_results: List[Tuple[str, str]] = field(default_factory=list)


def _decompose_rtt(
    rtt: str,
    *,
    cargo: str,
    arm_sense: str,
    pegrna_strand: str,
    nick_gt_s: bool,
) -> List[Tuple[str, str]]:
    """Split the RTT into its labelled sub-parts in RTT (5'->3') order.

    The RTT is ``orient_rtt_tail(sense_flap)`` where the sense flap is
    ``cargo + arm`` (nick 3' of S) or ``arm + cargo`` (nick 5' of S). Returns
    ``[(label, seq), ...]`` with labels ``"insert"`` (the cargo) and
    ``"replacement_gdna"`` (the recoded arm), ordered so the sequences concatenate
    back to ``rtt`` exactly (verified). A reverse-targeting pegRNA RC's the flap,
    which reverses the segment order; this handles both. Empty segments (e.g. a
    zero-length arm) are omitted.
    """
    sense_segments = (
        [("insert", cargo), ("replacement_gdna", arm_sense)]
        if nick_gt_s
        else [("replacement_gdna", arm_sense), ("insert", cargo)]
    )
    oriented = [
        (label, orient_rtt_tail(seq, pegrna_strand=pegrna_strand))
        for label, seq in sense_segments
        if seq
    ]
    if "".join(s for _, s in oriented) != rtt:  # Rv RC reverses order -> reorder
        oriented = list(reversed(oriented))
    if "".join(s for _, s in oriented) != rtt:  # pragma: no cover - guard
        raise AssertionError("RTT decomposition does not reconstruct the RTT")
    return oriented


def _frame_at(cds_mask: Sequence[bool], pos: int) -> int:
    """Reading-frame offset of genomic ``pos`` = (# CDS bases from CDS start to
    pos) % 3. Assumes a contiguous CDS in the local window (no intron between the
    CDS start and ``pos``), which holds for the short re-templated arm."""
    try:
        cds_start = list(cds_mask).index(True)
    except ValueError:
        return 0
    return max(0, pos - cds_start) % 3


def _importance_metrics(scorer, residues: Sequence[int]) -> Dict[str, Any]:
    """Summarise the footprint residues' AlphaFold metrics vs the block cutoffs
    (worst RSA / max contact / min pLDDT + any special-residue flag) for the
    flowchart. Defensive: a scorer failure yields an empty (non-numeric) summary."""
    try:
        ms = [scorer.metrics(i) for i in residues]
    except Exception:  # noqa: BLE001
        return {"metric": "buried RSA", "value": None, "note": "metrics unavailable"}
    if not ms:
        return {"metric": "buried RSA", "value": None}
    min_rsa = min(m.rsa for m in ms)
    max_contact = max(m.contact_number for m in ms)
    special = any(
        m.is_catalytic or m.is_disulfide or m.is_positive_phi_gly or m.is_structural_pro
        for m in ms
    )
    # Primary flagged value = worst (lowest) RSA vs the buried cutoff 0.25; a
    # residue passes when it is NOT buried (RSA >= 0.25) and not special.
    return {
        "metric": "min RSA", "value": round(min_rsa, 2), "cutoff": 0.25, "cmp": ">=",
        "extra": f"max contacts {max_contact} (≤18)" + (" · special residue!" if special else ""),
    }


def _insertion_site_footprint(
    site_s: int, cds_mask: Sequence[bool], *, context: int = 1, protein_length: Optional[int] = None
) -> List[int]:
    """Residue indices flanking the cargo insertion point ``site_s`` (ADR 0013).

    The cargo (tag + linker + T7-leftover junk) is inserted at ``site_s``; the
    coding residues it physically abuts -- the residue at ``site_s - 1`` (5') and at
    ``site_s`` (3', if coding), plus ``context`` residues each side -- are the ones
    whose structural context the insert could disrupt. This is the correct target of
    the structural-importance gate: forbid inserting into a buried / interaction-prone
    region. It is NOT the synonymous recode arm (INV-0 guarantees the protein there is
    unchanged, so the arm residues' burial is irrelevant). Returns ``[]`` when the
    insertion point has no flanking CDS. Residue index = (pos - cds_start) // 3.
    """
    mask = list(cds_mask)
    if not any(mask):
        return []
    cds0 = mask.index(True)
    residues = set()
    for p in (site_s - 1, site_s):
        if 0 <= p < len(mask) and mask[p]:
            residues.add((p - cds0) // 3)
    if not residues:
        return []
    # Map local-window residue indices to ABSOLUTE protein residues. The local
    # design window rarely starts at protein residue 1, so ``(pos-cds0)//3`` is only
    # a window-relative index; querying the AFDB structure with it reads the wrong
    # residue (a C-terminal knock-in mis-reads as an N-terminal core -> false
    # "buried" block). Given ``protein_length``, anchor the window's LAST coding
    # residue to the protein C-terminus (0-based last residue = protein_length-1).
    if protein_length:
        last_true = len(mask) - 1 - mask[::-1].index(True)
        offset = (protein_length - 1) - ((last_true - cds0) // 3)
        residues = {r + offset for r in residues}
    lo, hi = min(residues), max(residues)
    span = range(max(0, lo - context), hi + context + 1)
    if protein_length:
        return [i for i in span if 0 <= i <= protein_length - 1]
    return list(span)


def _derive_pbs1(genome: str, *, site_s: int, pbs_length: int, pegrna_strand: str) -> str:
    """PBS1 primes off the 3'-OH the pegRNA nick creates (the pegRNA nick is at S,
    decision 1). It base-pairs antiparallel to the L genomic bases ending at the
    nick, so PBS1 = reverse_complement(genome[S-L:S]) for a reverse-targeting
    pegRNA (literal window for forward). This is the mechanistically-correct
    derivation from the genome — NOT a slice of the spacer. In production the exact
    PBS comes from PRIDICT; this is the fallback when designing from genome.
    """
    from Bio.Seq import Seq  # noqa: PLC0415

    window = genome[max(0, site_s - pbs_length) : site_s].upper()
    if pegrna_strand == "Rv":
        return str(Seq(window).reverse_complement())
    return window


def _assemble_edit(
    *, genome: str, site_s: int, nick_position: int, arm_sense: str, cargo: str
) -> Tuple[str, str, int, int]:
    """Build the sense flap + predicted edited allele for a cargo-at-S insert.

    The cargo is ALWAYS inserted at ``site_s``. The recode/deletion-replacement
    arm sits between S and the nick:

    * nick 3' of S  -> flap = ``cargo + arm`` ; arm replaces ``genome[S:nick]``.
    * nick 5' of S  -> flap = ``arm + cargo`` ; arm replaces ``genome[nick:S]``
      (C-terminal recode-into-CDS: the last codons, 5' of S, are recoded).

    Returns ``(sense_flap, edited, tprt_site, tprt_hi)`` where ``tprt_site``/
    ``tprt_hi`` are the boundaries handed to :func:`simulate_tprt` so the flap
    replaces ``genome[tprt_site:tprt_hi]``.
    """
    if nick_position > site_s:
        sense_flap = cargo + arm_sense
        edited = genome[:site_s] + cargo + arm_sense + genome[nick_position:]
        return sense_flap, edited, site_s, nick_position
    sense_flap = arm_sense + cargo
    edited = genome[:nick_position] + arm_sense + cargo + genome[site_s:]
    return sense_flap, edited, nick_position, site_s


def _arm_protein_equal(original_arm: str, recoded_arm: str, frame: int) -> bool:
    """True iff the recoded arm encodes the SAME protein as the original over its
    codon-aligned interior (independent synonymity check for INV-0)."""
    from Bio.Seq import Seq  # noqa: PLC0415

    pad = (3 - frame) % 3
    a = original_arm[pad:]
    b = recoded_arm[pad:]
    a = a[: (len(a) // 3) * 3]
    b = b[: (len(b) // 3) * 3]
    return str(Seq(a).translate()) == str(Seq(b).translate())


@dataclass(frozen=True)
class ArmBuild:
    """Result of partitioning the re-templated arm ``genome[lo:hi)`` by cds_mask.

    CDS runs are synonymously recoded; non-coding runs are kept byte-exact
    (Decision 15). ``sense_arm`` is the rebuilt arm (5'->3', sense/window frame).
    ``recode_span`` / ``locked_span`` are genome-coordinate spans handed to
    :func:`assert_inv0`; a non-coding (locked) run sets ``microhomology_flag``.
    """

    sense_arm: str
    recode_span: Optional[Tuple[int, int]]
    locked_span: Optional[Tuple[int, int]]
    locked_expected: Optional[str]
    n_edits: int
    max_run_after: int
    microhomology_flag: bool
    protein_ok: bool


def _partition_arm(
    genome: str, *, lo: int, hi: int, cds_mask: Sequence[bool], max_run: int
) -> Optional[ArmBuild]:
    """Partition ``genome[lo:hi)`` into maximal CDS / non-coding runs and build the
    re-templated arm: recode CDS runs (synonymous, microhomology ``<= max_run``),
    keep non-coding runs byte-exact. Returns ``None`` (fail-loud upstream) if a CDS
    run has no synonymous solution or the arm crosses more than one coding boundary
    (> 2 runs; not a C-terminal geometry)."""

    def _is_cds(p: int) -> bool:
        return bool(cds_mask[p]) if p < len(cds_mask) else False

    runs: List[Tuple[bool, int, int]] = []
    i = lo
    while i < hi:
        cls = _is_cds(i)
        j = i + 1
        while j < hi and _is_cds(j) == cls:
            j += 1
        runs.append((cls, i, j))
        i = j
    if len(runs) > 2:
        return None  # non-monotone arm (CDS<->UTR<->CDS) unsupported

    parts: List[str] = []
    n_edits = 0
    max_run_after = 0
    recode_span: Optional[Tuple[int, int]] = None
    locked_span: Optional[Tuple[int, int]] = None
    microhomology = False
    for cls, a, b in runs:
        seg = genome[a:b].upper()
        if cls:
            frame = _frame_at(cds_mask, a)
            rec = recode_arm(seg, genome=seg, frame=frame, max_run=max_run)
            if rec is None or len(rec.recoded) != len(seg):
                return None
            if not _arm_protein_equal(seg, rec.recoded, frame):
                return None
            parts.append(rec.recoded)
            n_edits += rec.n_edits
            max_run_after = max(max_run_after, rec.max_run_after)
            recode_span = (
                (a, b) if recode_span is None else (min(recode_span[0], a), max(recode_span[1], b))
            )
        else:
            parts.append(seg)
            locked_span = (
                (a, b) if locked_span is None else (min(locked_span[0], a), max(locked_span[1], b))
            )
            microhomology = True
    locked_expected = genome[locked_span[0] : locked_span[1]].upper() if locked_span else None
    return ArmBuild(
        sense_arm="".join(parts),
        recode_span=recode_span,
        locked_span=locked_span,
        locked_expected=locked_expected,
        n_edits=n_edits,
        max_run_after=max_run_after,
        microhomology_flag=microhomology,
        protein_ok=True,
    )


def run_recoded_tjpe(
    *,
    spacer: str,
    rtt_aa: str,
    cargo: str,
    genome: str,
    site_s: int,
    cds_mask: Sequence[bool],
    pbs_lengths: Tuple[int, ...] = (13, 15, 17),
    ngrna_distance_min: int = 40,
    ngrna_distance_max: int = 70,
    pegrna_strand: str = "Rv",
    arm_side: str = "3prime",
    max_run: int = 8,
    max_rtt_variants: int = 20,
    # injected scorers (fakes in CI, real adapters on CX3)
    splice_scorer,
    full_structure_scorer,
    importance_scorer,
    dock_scorer,
    offtarget_scorer: Callable[[List[str]], Dict[str, Tuple[float, int, bool]]],
    maxent_scorer=None,
    cas9: str = "nCas9",
    allele_hunt: bool = False,
    het_snps: Optional[Sequence[dict]] = None,
    structure_threshold: float = 0.5,
    max_interface_drop: float = 0.15,
    insertion_importance_context: int = 1,
    protein_length: Optional[int] = None,
    on_gate: Optional[Callable[..., None]] = None,
) -> "RecodedTJPEResult":
    """Ranked-pool-walk orchestrator (decisions 0-15).

    Builds a ranked candidate pool once; for each candidate runs the gates in
    cheap->expensive order (recode + INV-0 + rcPBS2 + ViennaRNA-class structure;
    then splice; then AlphaFold importance; then RiNALMo 2D structure; then
    Boltz-2 docking LAST). The first candidate passing every gate wins. The
    insert site ``site_s`` NEVER moves. Exhausting the pool raises
    :class:`NoValidPegRNA` with a per-candidate diagnostic.
    """
    if arm_side not in ("3prime", "5prime", "auto"):
        raise ValueError(f"arm_side must be '3prime'/'5prime'/'auto'; got {arm_side!r}")

    def _emit(cand: int, gate_id: str, info: str, passed: bool, metrics: Optional[dict] = None) -> None:
        # Live progress hook (per-gate, per-candidate). Never let it break a run.
        # ``metrics`` (raw value + cutoff for the flowchart) is passed when the
        # callback accepts it; older 4-arg callbacks are called without it.
        if on_gate:
            try:
                on_gate(cand, gate_id, info, passed, metrics)
            except TypeError:
                try:
                    on_gate(cand, gate_id, info, passed)
                except Exception:  # noqa: BLE001
                    pass
            except Exception:  # noqa: BLE001
                pass

    het_snps = list(het_snps or [])
    # Advisory only — must never crash the design run on a malformed SNP record.
    try:
        allele_options = find_allele_specific(het_snps) if allele_hunt else []
    except Exception:  # noqa: BLE001
        allele_options = []

    # Nick candidates on the opposite strand, off-target screened (incl. other-allele).
    nicks = enumerate_nick_sites(
        genome,
        site_s=site_s,
        pegrna_strand=pegrna_strand,
        distance_min=ngrna_distance_min,
        distance_max=ngrna_distance_max,
    )
    raw_nick_count = len(nicks)
    nicks = screen_nick_candidates(nicks, scorer=offtarget_scorer)

    variants = enumerate_pbs_rtt_variants(
        spacer=spacer,
        rtt_aa=rtt_aa,
        pbs_lengths=pbs_lengths,
        max_rtt_variants=max_rtt_variants,
    )

    # Warm the apo baseline once; the real Boltz-2 adapter caches it internally so
    # repeated passes_docking_gate calls across candidates reuse it. Guard it: a
    # scorer that cannot run must not crash the whole walk before any candidate.
    diagnostics: List[str] = []
    try:
        dock_scorer.apo_interface_iptm(cas9)
    except Exception as exc:  # noqa: BLE001 - scorer/runtime failures are diagnostics
        raise NoValidPegRNA(f"docking apo baseline failed to compute: {exc}") from exc
    tried = 0

    # Diagnosability: an empty nick pool means the candidate loop never runs, so the
    # final NoValidPegRNA would otherwise carry a blank diagnostic. Say WHY here.
    if not nicks:
        if raw_nick_count == 0:
            diagnostics.append(
                f"0 nick sites: no NGG-PAM protospacer in the "
                f"[{ngrna_distance_min},{ngrna_distance_max}] bp window on the strand "
                f"opposite the pegRNA ({pegrna_strand})"
            )
        else:
            diagnostics.append(
                f"0 nicks after off-target screen: all {raw_nick_count} PAM candidate(s) "
                f"were rejected (each needs a unique on-target and 0 off-targets within the "
                f"gate mismatch tier)"
            )

    for nick in nicks:
        pos = nick.nick_position
        if pos == site_s:
            diagnostics.append(f"nick {pos}: coincides with S")
            continue
        # arm_side selects which side of S the recode/replacement arm sits on.
        # "3prime" (default) preserves the historical contract (arm 3' of S);
        # "5prime" puts the arm into the CDS 5' of S (C-terminal recode-into-CDS).
        if arm_side == "3prime" and pos < site_s:
            diagnostics.append(f"nick {pos}: 5'-of-S excluded (arm_side=3prime)")
            continue
        if arm_side == "5prime" and pos > site_s:
            diagnostics.append(f"nick {pos}: 3'-of-S excluded (arm_side=5prime)")
            continue
        # Re-templated span between S and the nick (the deletion-replacement arm).
        lo, hi = sorted((site_s, pos))
        # Partition the arm by cds_mask: recode CDS runs, keep non-coding runs
        # byte-exact (Decision 15). Replaces the unconditional whole-arm recode.
        arm = _partition_arm(genome, lo=lo, hi=hi, cds_mask=cds_mask, max_run=max_run)
        if arm is None:
            diagnostics.append(
                f"nick {nick.nick_position}: arm partition failed "
                "(no synonymous CDS solution or non-monotone arm)"
            )
            continue

        # rcPBS2 depends only on the nick (not the RTT variant): choose its LENGTH
        # once per nick to hit the Tm band (decision 8), decoupled from PBS1's
        # length. A per-nick miss here skips this nick, not the whole run.
        rcpbs2, rcpbs2_why = _select_rcpbs2(
            genome,
            nick_position=pos,
            pegrna_strand=pegrna_strand,
            tm_lo=THRESHOLDS.pbs2_tm_lo,
            tm_hi=THRESHOLDS.pbs2_tm_hi,
            gc_lo=THRESHOLDS.pbs2_gc_lo,
            gc_hi=THRESHOLDS.pbs2_gc_hi,
            len_min=THRESHOLDS.pbs2_len_min,
            len_max=THRESHOLDS.pbs2_len_max,
        )
        if rcpbs2 is None:
            diagnostics.append(f"nick {pos}: rcPBS2 select: {rcpbs2_why}")
            continue

        for variant in variants:
            tried += 1
            # Sense flap + predicted allele, ordered by which side of S the arm
            # sits (cargo always at S). Then orient per strand (Rv -> rev-comp).
            sense_flap, edited, tprt_site, tprt_hi = _assemble_edit(
                genome=genome,
                site_s=site_s,
                nick_position=pos,
                arm_sense=arm.sense_arm,
                cargo=cargo,
            )
            rtt = orient_rtt_tail(sense_flap, pegrna_strand=pegrna_strand)
            # PBS2 already selected per-nick above (length scanned to hit Tm band).
            # INV-0 on the predicted product. recode_span (CDS) changes must be
            # protein-silent (verified in _partition_arm); locked_span (non-coding)
            # must be byte-exact to the targeted allele (Decision 15).
            reference = genome

            def _protein_check(_r, _e, _ok=arm.protein_ok):
                return _ok

            try:
                assert_inv0(
                    reference=reference,
                    edited=edited,
                    site_s=site_s,
                    cargo=cargo,
                    recode_span=arm.recode_span if arm.recode_span is not None else (lo, lo),
                    protein_check=_protein_check,
                    allowed_silent=None,
                    locked_span=arm.locked_span,
                    locked_expected=arm.locked_expected,
                )
            except ValueError as exc:
                diagnostics.append(f"cand{tried}: INV-0: {exc}")
                continue
            # Strand-correctness: simulate TPRT from the ORIENTED RTT tail and
            # assert it reconstructs the intended sense allele. Catches any
            # complement/reverse/reverse_complement error loudly.
            try:
                simulate_tprt(
                    genome=genome,
                    site_s=tprt_site,
                    hi=tprt_hi,
                    rtt_tail_on_peg=rtt,
                    pegrna_strand=pegrna_strand,
                    expected_allele=edited,
                )
            except TPRTMismatch as exc:
                diagnostics.append(f"cand{tried}: TPRT strand mismatch: {exc}")
                continue
            # PBS1 derived from the genome at the pegRNA nick (RC of the nicked
            # 3' end) — NOT a slice of the spacer.
            pbs1 = _derive_pbs1(
                genome, site_s=site_s, pbs_length=variant.pbs_length, pegrna_strand=pegrna_strand
            )
            peg = assemble_tj_pegrna_v2(
                spacer=spacer,
                pbs2=rcpbs2,
                rtt=rtt,
                pbs1=pbs1,
                polyu="UUU",
                scaffold=STANDARD_SGRNA_SCAFFOLD,
            )
            # Per-gate info for the winning candidate (surfaced in the report + web
            # flowchart so each step shows the score/verdict it produced).
            gate_info: Dict[str, str] = {}
            # --- expensive gates last; each is GUARDED so a scorer crash skips the
            # candidate (with a diagnostic) instead of aborting the whole walk. ---
            # PRIMARY splice gate: MaxEntScan cryptic-splice-site scan of the edited
            # construct (recode span + cassette + junctions). Context-free, so it is
            # the right instrument for a large exogenous-cassette insertion (ADR 0011).
            if maxent_scorer is not None:
                ms_metrics: Dict[str, Any] = {}
                try:
                    ms_ok, ms_why = passes_cryptic_splice_gate(
                        reference, edited, scorer=maxent_scorer, metrics=ms_metrics
                    )
                except Exception as exc:  # noqa: BLE001
                    diagnostics.append(f"cand{tried}: cryptic-splice scorer error: {exc}")
                    continue
                _emit(tried, "cryptic", ms_why, ms_ok, ms_metrics)
                if not ms_ok:
                    diagnostics.append(f"cand{tried}: cryptic splice: {ms_why}")
                    continue
                gate_info["cryptic"] = ms_why

            # ADVISORY cross-check: SpliceAI needs >=10 kb of real genomic context. When
            # a caller supplies it, a concordant disruption call also fails the candidate;
            # on short context it raises ValueError and is skipped (NOT fail-loud) so the
            # MaxEnt gate above remains the primary check (ADR 0011).
            sp_metrics: Dict[str, Any] = {}
            try:
                sp_ok, sp_why = passes_splice_gate(
                    reference, edited, scorer=splice_scorer, metrics=sp_metrics
                )
            except ValueError:
                gate_info["splice"] = "advisory skip (needs >=10kb context; MaxEnt is primary)"
                _emit(tried, "splice", "advisory skip (<10kb context)", True,
                      {"metric": "SpliceAI Δ", "value": None, "note": "advisory skip (<10kb context)"})
            except Exception as exc:  # noqa: BLE001
                diagnostics.append(f"cand{tried}: splice cross-check error: {exc}")
                continue
            else:
                _emit(tried, "splice", sp_why, sp_ok, sp_metrics)
                if not sp_ok:
                    diagnostics.append(f"cand{tried}: splice (SpliceAI cross-check): {sp_why}")
                    continue
                gate_info["splice"] = sp_why

            # Insertion-site structural tolerance (ADR 0013): the cargo (tag + linker
            # + T7-leftover junk) is inserted at S and must NOT land in a structured /
            # interaction-prone region, or the insert disrupts the fold. Judge the
            # residues flanking the insertion point -- NOT the synonymous recode arm
            # (INV-0 guarantees the protein there is unchanged, so the arm residues'
            # burial is irrelevant). Retains the special-residue hard block.
            ins_footprint = _insertion_site_footprint(
                site_s, cds_mask, context=insertion_importance_context, protein_length=protein_length
            )
            if ins_footprint:
                try:
                    blocked = span_blocked(importance_scorer, ins_footprint)
                except Exception as exc:  # noqa: BLE001
                    diagnostics.append(f"cand{tried}: importance scorer error: {exc}")
                    continue
                imp_metrics = _importance_metrics(importance_scorer, ins_footprint)
                _emit(
                    tried,
                    "importance",
                    f"residues {ins_footprint} " + ("blocked (buried/special)" if blocked else "tolerated (surface/flexible)"),
                    not blocked,
                    imp_metrics,
                )
                if blocked:
                    diagnostics.append(
                        f"cand{tried}: insertion site S residue(s) {ins_footprint} "
                        "structurally important (buried/packed/special) -- inserting a "
                        "tag/linker here risks disrupting the protein; choose a "
                        "surface/loop/terminal insertion site"
                    )
                    continue
                gate_info["importance"] = f"residues {ins_footprint} tolerated (surface/flexible)"
            else:
                gate_info["importance"] = "no flanking CDS at insertion site (n/a)"
                _emit(tried, "importance", "no flanking CDS (n/a)", True)

            regions = _regions_for(spacer, rcpbs2, rtt, pbs1)
            try:
                pairing = functional_region_pairing(
                    peg,
                    regions=regions,
                    scorer=full_structure_scorer,
                    threshold=structure_threshold,
                )
            except Exception as exc:  # noqa: BLE001
                diagnostics.append(f"cand{tried}: structure scorer error: {exc}")
                continue
            _emit(
                tried,
                "structure",
                f"{pairing.cross_region_pairs} cross-region pairs (p>={structure_threshold})",
                pairing.cross_region_pairs == 0,
                {"metric": "cross-region pairs", "value": pairing.cross_region_pairs,
                 "cutoff": 0, "cmp": "==", "extra": f"base-pair prob ≥ {structure_threshold}"},
            )
            if pairing.cross_region_pairs > 0:
                diagnostics.append(f"cand{tried}: RiNALMo cross-region pairing")
                continue
            gate_info["structure"] = (
                f"{pairing.cross_region_pairs} cross-region pairs "
                f"(threshold p>={structure_threshold})"
            )

            dk_metrics: Dict[str, Any] = {}
            try:
                dock_ok, dock_why = passes_docking_gate(
                    peg, cas9, scorer=dock_scorer, max_drop=max_interface_drop,
                    metrics=dk_metrics,
                )
            except Exception as exc:  # noqa: BLE001
                diagnostics.append(f"cand{tried}: docking scorer error: {exc}")
                continue
            _emit(tried, "docking", dock_why, dock_ok, dk_metrics)
            if not dock_ok:
                diagnostics.append(f"cand{tried}: docking: {dock_why}")
                continue
            gate_info["docking"] = dock_why
            # First survivor wins.
            warnings: List[str] = []
            if arm.max_run_after > max_run:
                warnings.append("microhomology advisory: residual run above target")
            if arm.microhomology_flag:
                warnings.append(
                    "microhomology advisory: non-coding locked span re-templated "
                    "byte-exact (cannot recode; Decision 15)"
                )
            rtt_components = _decompose_rtt(
                rtt,
                cargo=cargo,
                arm_sense=arm.sense_arm,
                pegrna_strand=pegrna_strand,
                nick_gt_s=pos > site_s,
            )
            parts_map = dict(rtt_components)
            gate_info["offtarget"] = f"nick {nick.nick_position}: unique on-target, 0 off-targets in tier"
            gate_info["recode"] = f"{arm.n_edits} synonymous codon edits (microhomology broken)"
            gate_info["inv0"] = "0 bp endogenous deletion; recode protein-silent (INV-0 held)"
            # Order gate_info to match the pipeline flow for display.
            _order = ["offtarget", "recode", "inv0", "cryptic", "splice", "importance", "structure", "docking"]
            gate_results = [(k, gate_info[k]) for k in _order if k in gate_info]
            return RecodedTJPEResult(
                pegRNA_sequence=peg,
                spacer=spacer,
                rcpbs2=rcpbs2,
                rtt=rtt,
                pbs1=pbs1,
                polyu="UUU",
                ngrna_protospacer=nick.protospacer,
                site_s=site_s,
                nick_position=nick.nick_position,
                deletion_size=0,
                recode_edits=arm.n_edits,
                candidates_tried=tried,
                warnings=warnings,
                allele_specific_options=allele_options,
                rtt_components=rtt_components,
                insert_seq=parts_map.get("insert", ""),
                replacement_gdna=parts_map.get("replacement_gdna", ""),
                edited_allele=edited,
                recode_span=arm.recode_span,
                locked_span=arm.locked_span,
                gate_results=gate_results,
            )

    raise NoValidPegRNA(
        f"no candidate survived all gates after {tried} tried; "
        f"diagnostics: {'; '.join(diagnostics[:10])}"
    )


def _regions_for(spacer: str, rcpbs2: str, rtt: str, pbs1: str) -> Dict[str, Tuple[int, int]]:
    """Region offsets in the assembled pegRNA (authoritative order, scaffold excluded
    from cross-region scoring by name)."""
    scaf = len(STANDARD_SGRNA_SCAFFOLD)
    i0 = 0
    i1 = i0 + len(spacer)
    i2 = i1 + scaf
    i3 = i2 + len(rcpbs2)
    i4 = i3 + len(rtt)
    i5 = i4 + len(pbs1)
    return {
        "spacer": (i0, i1),
        "scaffold": (i1, i2),
        "rcpbs2": (i2, i3),
        "rtt": (i3, i4),
        "pbs1": (i4, i5),
    }
