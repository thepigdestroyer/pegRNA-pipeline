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
from pegrna_pipeline.hpc.complex_dock import DockingServiceError, passes_docking_gate
from pegrna_pipeline.hpc.cryptic_splice import passes_cryptic_splice_gate
from pegrna_pipeline.hpc.full_structure import (
    functional_region_pairing,
    scaffold_extension_pairing,
)
from pegrna_pipeline.hpc.gate_thresholds import THRESHOLDS
from pegrna_pipeline.hpc.guide_sanity import assess_spacer
from pegrna_pipeline.hpc.nick_search import (
    enumerate_nick_sites,
    filter_nick_sanity,
    screen_nick_candidates,
)
from pegrna_pipeline.hpc.pridict_ngrna import InjectedNgRNA
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
    per-candidate diagnostic of which gate killed each one (``.diagnostics``), which the
    INSERT-lever wrapper inspects to decide whether a synonymous cargo swap could help."""

    def __init__(self, message: str, diagnostics: Optional[Sequence[str]] = None) -> None:
        super().__init__(message)
        self.diagnostics: List[str] = list(diagnostics or [])


@dataclass(frozen=True)
class PbsRttVariant:
    spacer: str
    rtt_dna: str
    pbs_length: int


def enumerate_pbs_rtt_variants(
    *,
    spacer: str,
    rtt_aa: Optional[str] = None,
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

    ``rtt_aa`` is OPTIONAL and does not change any pegRNA: the walk uses only a variant's PBS
    length (the RTT is the cargo plus the genome-derived arm), so five different anchor residues
    produce byte-identical designs. Without it there is one variant per PBS length.
    """
    rtts = (
        enumerate_variants(rtt_aa, max_variants=max_rtt_variants, ranking=ranking)
        if rtt_aa
        else [None]
    )
    out: List[PbsRttVariant] = []
    for rv in rtts:
        for length in pbs_lengths:
            out.append(
                PbsRttVariant(spacer=spacer, rtt_dna=rv.dna if rv else "", pbs_length=length)
            )
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
        m.is_catalytic or m.is_disulfide or m.is_positive_phi_gly or m.is_structural_pro for m in ms
    )
    # Primary flagged value = worst (lowest) RSA vs the buried cutoff 0.25; a
    # residue passes when it is NOT buried (RSA >= 0.25) and not special.
    return {
        "metric": "min RSA",
        "value": round(min_rsa, 2),
        "cutoff": 0.25,
        "cmp": ">=",
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


def _resolve_protein_length(protein_length: Optional[int], importance_scorer: Any) -> Optional[int]:
    """Absolute protein length used to map the local insertion-site footprint onto the
    AlphaFold structure.

    An explicit ``protein_length`` wins. Otherwise ask the importance scorer for its own
    structure's residue count (:meth:`AFDBImportanceScorer.protein_length`) -- the
    authoritative source, since that is the exact structure the gate indexes -- so a real
    run maps local->absolute WITHOUT the caller having to supply the length. Returns
    ``None`` when neither is available (a frame-agnostic fake scorer), leaving the
    footprint window-local; that is harmless because such a scorer ignores the residue
    index. This removes the silent mis-index where a C-terminal knock-in was read as an
    N-terminal core because the local window index was queried against the full protein.
    """
    if protein_length is not None:
        return protein_length
    getter = getattr(importance_scorer, "protein_length", None)
    if callable(getter):
        try:
            value = getter()
        except Exception:  # noqa: BLE001
            return None
        return int(value) if value else None
    return None


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


def _advisory_splice_text(
    sp_ok: bool, metrics: Dict[str, Any], where: str, *, primary_passed: bool
) -> str:
    """What the report says about a SpliceAI number computed without real genomic context.

    Without that context the number cannot tell a clean large insert from a bad one (ADR 0028), so it never
    rejects. Under the cutoff it is reported with where it came from. At or above it, it is a warning that
    says plainly that it is advisory, why, and what the primary gate found.
    """
    value = metrics.get("value")
    shown = f"{value:.3f}" if isinstance(value, (int, float)) else "?"
    where = where or "no genomic context"
    if sp_ok:
        return f"splice delta {shown} ok (advisory only: {where})"
    text = (
        f"SpliceAI delta {shown} is at or above the {metrics.get('cutoff')} cutoff, but this is "
        f"advisory only ({where}): without real flanks SpliceAI cannot tell a clean large insert "
        "from a bad one, so the number is reported and does not reject."
    )
    if primary_passed:
        text += " The primary splice gate (MaxEntScan) passed."
    return text


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


@dataclass(frozen=True)
class _NickPlan:
    """One second-nick option fed to the pool walk. ``rcpbs2`` is pre-derived for the
    PRIDICT (injected) path; ``None`` means the genome path selects it per-nick via
    :func:`_select_rcpbs2`."""

    nick_position: int
    protospacer: str
    rcpbs2: Optional[str] = None
    #: Off-target sites the screen counted for this guide; ``None`` if it was not screened.
    offtargets: Optional[int] = None


def _run_recoded_tjpe_single(
    *,
    spacer: str,
    rtt_aa: Optional[str] = None,
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
    structure_threshold: Optional[float] = None,
    max_interface_drop: float = 0.15,
    insertion_importance_context: int = 1,
    protein_length: Optional[int] = None,
    on_gate: Optional[Callable[..., None]] = None,
    injected_ngrnas: Optional[Sequence[InjectedNgRNA]] = None,
    ngrna_max_offtargets: int = 0,
) -> "RecodedTJPEResult":
    """Ranked-pool-walk orchestrator (decisions 0-15).

    ``ngrna_max_offtargets`` is how many off-target sites (extra exact copies plus near matches
    within the scorer's mismatch tier) the second-nick guide may have anywhere in the genome.
    The default, 0, is ADR 0010's zero tolerance. It applies to the genome-enumerated path
    here; the PRIDICT path applies the same number where it screens
    (:func:`pegrna_pipeline.hpc.pridict_ngrna.build_pridict_ngrnas`). Whatever count the chosen
    guide actually has is reported on the design, so a tolerated off-target is never silent.

    Builds a ranked candidate pool once; for each candidate runs the gates in
    cheap->expensive order (recode + INV-0 + rcPBS2 + ViennaRNA-class structure;
    then splice; then AlphaFold importance; then RiNALMo 2D structure; then
    Boltz-2 docking LAST). The first candidate passing every gate wins. The
    insert site ``site_s`` NEVER moves. Exhausting the pool raises
    :class:`NoValidPegRNA` with a per-candidate diagnostic.
    """
    if arm_side not in ("3prime", "5prime", "auto"):
        raise ValueError(f"arm_side must be '3prime'/'5prime'/'auto'; got {arm_side!r}")
    if int(ngrna_max_offtargets) != ngrna_max_offtargets or ngrna_max_offtargets < 0:
        raise ValueError(
            f"ngrna_max_offtargets must be an integer >= 0; got {ngrna_max_offtargets!r}"
        )

    def _emit(
        cand: int, gate_id: str, info: str, passed: bool, metrics: Optional[dict] = None
    ) -> None:
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

    # Second-nick candidates. AUTO (PRIDICT) path: the pipeline already selected +
    # off-target-screened the ngRNA from PRIDICT nicking guides and derived each rcPBS2
    # (pridict_ngrna.build_pridict_ngrnas), so inject them verbatim -- iterating them IS
    # the NICK lever. GENOME (fallback / explicit-spacer) path: enumerate PAM-anchored
    # opposite-strand nicks and off-target screen them here (unchanged).
    unusable_nicks: list = []
    if injected_ngrnas is not None:
        plans = [
            _NickPlan(
                nick_position=inj.nick_position,
                protospacer=inj.protospacer,
                rcpbs2=inj.rcpbs2,
                offtargets=inj.offtargets,
            )
            for inj in injected_ngrnas
        ]
        raw_nick_count = len(plans)
    else:
        nicks = enumerate_nick_sites(
            genome,
            site_s=site_s,
            pegrna_strand=pegrna_strand,
            distance_min=ngrna_distance_min,
            distance_max=ngrna_distance_max,
        )
        raw_nick_count = len(nicks)
        # The same hard sequence rules the PRIDICT path applies, before the genome scan.
        nicks, unusable_nicks = filter_nick_sanity(nicks)
        nicks = screen_nick_candidates(
            nicks, scorer=offtarget_scorer, max_offtargets=ngrna_max_offtargets
        )
        plans = [
            _NickPlan(
                nick_position=n.nick_position, protospacer=n.protospacer, offtargets=n.offtargets
            )
            for n in nicks
        ]

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
        # The length-calibrated gate (scorers with per-sample ipTMs) compares each
        # candidate with references measured at its own length and never uses the apo
        # scaffold, so skip that call: it costs a NIM request and measures nothing it uses.
        if not hasattr(dock_scorer, "interface_iptm_samples"):
            dock_scorer.apo_interface_iptm(cas9)
    except DockingServiceError:
        raise
    except Exception as exc:  # noqa: BLE001 - scorer/runtime failures are diagnostics
        raise NoValidPegRNA(f"docking apo baseline failed to compute: {exc}") from exc
    tried = 0

    # Diagnosability: an empty nick pool means the candidate loop never runs, so the
    # final NoValidPegRNA would otherwise carry a blank diagnostic. Say WHY here.
    if not plans:
        if injected_ngrnas is not None:
            diagnostics.append(
                "0 injected PRIDICT ngRNA candidates (none passed the distance window / "
                "opposite-strand / off-target gates upstream)"
            )
        elif raw_nick_count == 0:
            diagnostics.append(
                f"0 nick sites: no NGG-PAM protospacer in the "
                f"[{ngrna_distance_min},{ngrna_distance_max}] bp window on the strand "
                f"opposite the pegRNA ({pegrna_strand})"
            )
        else:
            screened = raw_nick_count - len(unusable_nicks)
            parts = []
            if unusable_nicks:
                parts.append(
                    f"{len(unusable_nicks)} failed the sequence rules "
                    f"({'; '.join(sorted({r.split(':')[0] for _, rs in unusable_nicks for r in rs}))})"
                )
            if screened:
                parts.append(
                    f"{screened} failed the off-target screen (each needs a unique on-target and "
                    f"at most {ngrna_max_offtargets} off-target site(s) within the gate mismatch "
                    "tier)"
                )
            # The lead phrase is the one ADR 0010 documents; the breakdown says which rule.
            diagnostics.append(
                f"0 nicks after off-target screen: all {raw_nick_count} PAM candidate(s) were "
                "rejected: " + ", ".join(parts)
            )

    for plan in plans:
        pos = plan.nick_position
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
                f"nick {pos}: arm partition failed "
                "(no synonymous CDS solution or non-monotone arm)"
            )
            continue

        # rcPBS2. PRIDICT path: use the genome-derived rcPBS2 the pipeline already built
        # (derive_rc_pbs2, plan D4). GENOME path: choose its LENGTH once per nick to hit
        # the Tm band (decision 8), decoupled from PBS1's length; a per-nick miss skips
        # this nick, not the whole run.
        if plan.rcpbs2 is not None:
            rcpbs2 = plan.rcpbs2
        else:
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

        # The RTT variant never enters the assembled pegRNA: only its PBS length does (the RTT is
        # the cargo plus the genome-derived arm). So every variant of one PBS length at this nick
        # is the SAME design, and the pool of 6 variants x 3 PBS lengths is really 3 designs. A
        # design that failed once fails identically again (every gate is deterministic), and the
        # repeats were not free: each re-ran the docking gate, a ~40 s paid Boltz-2 request.
        seen_pbs_lengths: set = set()
        for variant in variants:
            if variant.pbs_length in seen_pbs_lengths:
                continue
            seen_pbs_lengths.add(variant.pbs_length)
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
            # MaxEnt gate above remains the primary check (ADR 0011). A scorer that had to
            # invent its context (N-padding) says so with ``binding = False``: its number is
            # then reported and, if high, warned about, but it does not reject (ADR 0028).
            sp_metrics: Dict[str, Any] = {}
            splice_advisory: Optional[str] = None
            try:
                sp_ok, sp_why = passes_splice_gate(
                    reference, edited, scorer=splice_scorer, metrics=sp_metrics
                )
            except ValueError:
                gate_info["splice"] = "advisory skip (needs >=10kb context; MaxEnt is primary)"
                _emit(
                    tried,
                    "splice",
                    "advisory skip (<10kb context)",
                    True,
                    {
                        "metric": "SpliceAI Δ",
                        "value": None,
                        "note": "advisory skip (<10kb context)",
                    },
                )
            except Exception as exc:  # noqa: BLE001
                diagnostics.append(f"cand{tried}: splice cross-check error: {exc}")
                continue
            else:
                binding = bool(getattr(splice_scorer, "binding", True))
                where = str(getattr(splice_scorer, "context_note", "") or "")
                if binding:
                    _emit(tried, "splice", sp_why, sp_ok, sp_metrics)
                    if not sp_ok:
                        diagnostics.append(f"cand{tried}: splice (SpliceAI cross-check): {sp_why}")
                        continue
                    gate_info["splice"] = f"{sp_why} ({where})" if where else sp_why
                else:
                    advisory = _advisory_splice_text(
                        sp_ok, sp_metrics, where, primary_passed="cryptic" in gate_info
                    )
                    gate_info["splice"] = advisory
                    if not sp_ok:
                        splice_advisory = advisory
                    _emit(
                        tried,
                        "splice",
                        advisory,
                        True,
                        {
                            "metric": "SpliceAI Δ",
                            "value": sp_metrics.get("value"),
                            "note": f"advisory, not binding: {where}",
                        },
                    )

            # Insertion-site structural tolerance (ADR 0013): the cargo (tag + linker
            # + T7-leftover junk) is inserted at S and must NOT land in a structured /
            # interaction-prone region, or the insert disrupts the fold. Judge the
            # residues flanking the insertion point -- NOT the synonymous recode arm
            # (INV-0 guarantees the protein there is unchanged, so the arm residues'
            # burial is irrelevant). Retains the special-residue hard block.
            ins_footprint = _insertion_site_footprint(
                site_s,
                cds_mask,
                context=insertion_importance_context,
                # Prefer the real scorer's own residue count over the (often omitted)
                # caller protein_length, so a C-terminal knock-in maps onto the correct
                # absolute residues instead of silently indexing local-window residues
                # into the full-protein structure.
                protein_length=_resolve_protein_length(protein_length, importance_scorer),
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
                    f"residues {ins_footprint} "
                    + ("blocked (buried/special)" if blocked else "tolerated (surface/flexible)"),
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
            structure_threshold = resolve_structure_threshold(
                structure_threshold, full_structure_scorer
            )
            try:
                # Fold once: the gate and the scaffold advisory read the same probability map
                # (a GPU folder would otherwise run twice per candidate).
                pair_map = full_structure_scorer.pair_probs(peg)
                pairing = functional_region_pairing(
                    peg,
                    regions=regions,
                    scorer=full_structure_scorer,
                    threshold=structure_threshold,
                    pairs=pair_map,
                )
                scaffold_pairing = scaffold_extension_pairing(
                    peg,
                    regions=regions,
                    scorer=full_structure_scorer,
                    threshold=structure_threshold,
                    pairs=pair_map,
                )
            except Exception as exc:  # noqa: BLE001
                diagnostics.append(f"cand{tried}: structure scorer error: {exc}")
                continue
            _emit(
                tried,
                "structure",
                f"{pairing.cross_region_pairs} cross-region pairs (p>={structure_threshold})",
                pairing.cross_region_pairs == 0,
                {
                    "metric": "cross-region pairs",
                    "value": pairing.cross_region_pairs,
                    "cutoff": 0,
                    "cmp": "==",
                    "extra": (
                        f"base-pair prob ≥ {structure_threshold} · extension pairs with "
                        f"{scaffold_pairing.paired_nt} scaffold nt"
                    ),
                },
            )
            if pairing.cross_region_pairs > 0:
                diagnostics.append(f"cand{tried}: RiNALMo cross-region pairing")
                continue
            gate_info["structure"] = (
                f"{pairing.cross_region_pairs} cross-region pairs "
                f"(threshold p>={structure_threshold}); extension pairs with "
                f"{scaffold_pairing.paired_nt} scaffold nt"
            )

            dk_metrics: Dict[str, Any] = {}
            try:
                dock_ok, dock_why = passes_docking_gate(
                    peg,
                    cas9,
                    scorer=dock_scorer,
                    max_drop=max_interface_drop,
                    metrics=dk_metrics,
                )
            except DockingServiceError:
                raise  # an outage is not a verdict on this candidate (see its docstring)
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
            if scaffold_pairing.paired_nt >= THRESHOLDS.scaffold_pairing_advisory_nt:
                # Advisory only. The folding gate leaves scaffold pairs out and the docking gate
                # sees only some of them, so nothing else flags this (ADR 0026).
                warnings.append(
                    f"the extension pairs with {scaffold_pairing.paired_nt} scaffold nucleotides "
                    f"(p >= {structure_threshold}; flagged from "
                    f"{THRESHOLDS.scaffold_pairing_advisory_nt}). Ordinary designs pair with far "
                    "fewer; designs built to tie up the scaffold pair with 18 or more. Neither "
                    "the folding gate nor the docking gate reliably catches this, so check the "
                    "predicted structure before ordering"
                )
            if splice_advisory:
                warnings.append("splice cross-check: " + splice_advisory)
            if dk_metrics.get("borderline"):
                warnings.append("docking gate: " + dk_metrics["note"])
            if dk_metrics.get("advisory"):
                # Outside the calibrated length range the docking gate can only report,
                # not judge. The design is accepted on the other gates, and says so.
                warnings.append("docking gate is advisory-only for this design: " + dock_why)
            if arm.max_run_after > max_run:
                warnings.append("microhomology advisory: residual run above target")
            if arm.microhomology_flag:
                warnings.append(
                    "microhomology advisory: non-coding locked span re-templated "
                    "byte-exact (cannot recode; Decision 15)"
                )
            # Spacer sequence findings. The ngRNA cannot fail a hard rule (it was filtered), but
            # the pegRNA spacer is the caller's or PRIDICT's choice, so a failure there is
            # reported rather than silently rejected.
            for label, guide in (("pegRNA spacer", spacer), ("ngRNA spacer", plan.protospacer)):
                verdict = assess_spacer(guide)
                warnings.extend(f"{label}: {f}" for f in verdict.failures)
                warnings.extend(f"{label} advisory: {a}" for a in verdict.advisories)
            tier = getattr(offtarget_scorer, "count_mismatch_threshold", None)
            tier_txt = f"within {tier} mismatches" if tier is not None else "in the mismatch tier"
            n_off = plan.offtargets or 0
            if n_off > 0:
                warnings.append(
                    f"second-nick guide has {n_off} off-target site(s) {tier_txt} "
                    f"(ngrna_max_offtargets={ngrna_max_offtargets}). Each is a single-strand "
                    "nick, not a double-strand break, but check them before ordering"
                )
            rtt_components = _decompose_rtt(
                rtt,
                cargo=cargo,
                arm_sense=arm.sense_arm,
                pegrna_strand=pegrna_strand,
                nick_gt_s=pos > site_s,
            )
            parts_map = dict(rtt_components)
            how = (
                "PRIDICT ngRNA, off-target-screened upstream"
                if plan.rcpbs2 is not None
                else "unique on-target"
            )
            gate_info["offtarget"] = (
                f"nick {plan.nick_position}: {how}, 0 off-target sites {tier_txt}"
                if n_off == 0
                else f"nick {plan.nick_position}: {how}, {n_off} off-target site(s) {tier_txt} "
                f"(tolerance {ngrna_max_offtargets})"
            )
            gate_info["recode"] = f"{arm.n_edits} synonymous codon edits (microhomology broken)"
            gate_info["inv0"] = "0 bp endogenous deletion; recode protein-silent (INV-0 held)"
            # Order gate_info to match the pipeline flow for display.
            _order = [
                "offtarget",
                "recode",
                "inv0",
                "cryptic",
                "splice",
                "importance",
                "structure",
                "docking",
            ]
            gate_results = [(k, gate_info[k]) for k in _order if k in gate_info]
            return RecodedTJPEResult(
                pegRNA_sequence=peg,
                spacer=spacer,
                rcpbs2=rcpbs2,
                rtt=rtt,
                pbs1=pbs1,
                polyu="UUU",
                ngrna_protospacer=plan.protospacer,
                site_s=site_s,
                nick_position=plan.nick_position,
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
        f"diagnostics: {'; '.join(diagnostics[:10])}",
        diagnostics=diagnostics,
    )


#: Diagnostic substrings that mean the failure is INSERT-side (the cargo SEQUENCE affects
#: the gate): a synonymous cargo swap can change these. Everything else is cargo-independent
#: -- importance (genomic insertion-site residues), off-target (spacer/PBS/ngRNA), INV-0 /
#: TPRT / arm / nick -- so a swap cannot help and must not trigger a retry.
_INSERT_SIDE_MARKERS = ("splice", "cryptic", "cross-region", "rinalmo", "structure")


def _is_insert_side_failure(diagnostics: Sequence[str]) -> bool:
    text = " ".join(diagnostics).lower()
    return any(marker in text for marker in _INSERT_SIDE_MARKERS)


def _cargo_variants(cargo: str, max_variants: int) -> List[str]:
    """The original cargo first, then up to ``max_variants``-1 CAI-ranked SYNONYMOUS
    variants of the whole cargo (the INSERT lever, Stage B/F). Variants are produced ONLY
    when the cargo is a clean in-frame CDS (length a multiple of 3, translatable, no internal
    stop); otherwise a synonymous swap is undefined and the lever no-ops (original only).
    A trailing stop codon is preserved verbatim. ``site_s`` and the spacer are unaffected by
    the swap (the spacer targets the genomic PAM, not the cargo)."""
    out = [cargo]
    if max_variants <= 1 or len(cargo) % 3 != 0:
        return out
    try:
        from Bio.Seq import Seq  # noqa: PLC0415

        aa = str(Seq(cargo).translate())
    except Exception:  # noqa: BLE001
        return out
    core = aa[:-1] if aa.endswith("*") else aa
    if not core or "*" in core:  # internal stop -> not a single clean CDS
        return out
    tail = cargo[len(core) * 3 :]  # trailing stop codon (if present), kept byte-exact
    from pegrna_pipeline.hpc.codon_expander import rank_insert_variants_by_cai  # noqa: PLC0415

    seen = {cargo}
    for cv in rank_insert_variants_by_cai(core, max_variants=max_variants * 4):
        dna = cv.dna + tail
        if dna in seen:
            continue
        seen.add(dna)
        out.append(dna)
        if len(out) >= max_variants:
            break
    return out


def run_recoded_tjpe(
    *, cargo: str, max_insert_variants: int = 6, **kwargs: Any
) -> "RecodedTJPEResult":
    """Canonical recoded tj-PE design workflow (the CLI/web entrypoint).

    Delegates to :func:`_run_recoded_tjpe_single` (the full gate cascade; see its docstring
    for every parameter) and adds the **INSERT lever** (unified plan Stage B/F): on an
    INSERT-side failure -- splice, cryptic-splice, or RNA structure, the gates a cargo
    sequence can change -- it retries with the next CAI-ranked SYNONYMOUS cargo variant. The
    edit site ``site_s`` and the spacer never change and the swap is protein-silent. A
    cargo-INDEPENDENT failure (importance / off-target / INV-0 / nick / arm) does NOT trigger
    a retry, because a synonymous swap cannot help it. ``max_insert_variants`` bounds the
    lever (>= 1; 1 disables it; the lever also no-ops unless the cargo is a clean in-frame
    CDS). Raises the last :class:`NoValidPegRNA` when every attempt is exhausted.
    """
    last: Optional[NoValidPegRNA] = None
    for cargo_variant in _cargo_variants(cargo, max_insert_variants):
        try:
            return _run_recoded_tjpe_single(cargo=cargo_variant, **kwargs)
        except NoValidPegRNA as exc:
            last = exc
            if not _is_insert_side_failure(exc.diagnostics):
                break  # a synonymous cargo swap cannot rescue a cargo-independent failure
    assert last is not None  # _cargo_variants always yields >= 1 attempt
    raise last


def resolve_structure_threshold(explicit: Optional[float], scorer: Any) -> float:
    """Cross-region pair-probability cutoff: explicit > the scorer's own > generic.

    A base-pair-probability boundary only means something relative to the folder that
    produced it, so the default lives on the scorer. An explicit value (config or
    argument) always wins.
    """
    if explicit is not None:
        return float(explicit)
    rec = getattr(scorer, "recommended_threshold", None)
    return float(rec) if rec is not None else THRESHOLDS.rinalmo_pair_prob


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
