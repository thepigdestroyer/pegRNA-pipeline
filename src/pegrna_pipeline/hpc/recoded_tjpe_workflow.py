"""Recoded TJ-pegRNA workflow: enumeration, gates, and the ranked-pool-walk.

This module wires together the pure building blocks (microhomology, recode_arm,
nick_search, INV-0) and the injected gate scorers (splice, importance, full
structure, complex docking) into the candidate enumeration and the
ranked-pool-walk orchestrator. Heavy scorers are passed in (fakes in CI, real
adapters on CX3); this module imports none of them at module load.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from pegrna_pipeline.hpc.allele_specificity import find_allele_specific
from pegrna_pipeline.hpc.codon_expander import _gc_content, enumerate_variants
from pegrna_pipeline.hpc.complex_dock import passes_docking_gate
from pegrna_pipeline.hpc.full_structure import functional_region_pairing
from pegrna_pipeline.hpc.nick_search import enumerate_nick_sites, screen_nick_candidates
from pegrna_pipeline.hpc.recode_arm import recode_arm
from pegrna_pipeline.hpc.scarless import assert_inv0
from pegrna_pipeline.hpc.splice_screen import passes_splice_gate
from pegrna_pipeline.hpc.strand_orientation import orient_pbs2, orient_rtt_tail
from pegrna_pipeline.hpc.structure_importance import span_blocked
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD, assemble_tj_pegrna_v2
from pegrna_pipeline.hpc.tprt_sim import TPRTMismatch, simulate_tprt

__all__ = [
    "PbsRttVariant",
    "enumerate_pbs_rtt_variants",
    "passes_pbs2_gate",
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
) -> List[PbsRttVariant]:
    """Enumerate (synonymous RTT) × (PBS1 length); the spacer is FIXED by site S.

    These are the candidates handed to PRIDICT2 (reverse-targeting) for ranking.
    The spacer is never varied — there is one (or a few) protospacer(s) that nick
    at S, and that placement is fixed (decision 1). Deterministic order:
    RTT-variant rank, then PBS length.
    """
    rtts = enumerate_variants(rtt_aa, max_variants=max_rtt_variants)
    out: List[PbsRttVariant] = []
    for rv in rtts:
        for length in pbs_lengths:
            out.append(PbsRttVariant(spacer=spacer, rtt_dna=rv.dna, pbs_length=length))
    return out


def passes_pbs2_gate(
    pbs2: str,
    *,
    repetitive: bool,
    min_len: int = 11,
    max_len: int = 17,
    gc_lo: float = 0.40,
    gc_hi: float = 0.70,
    tm_lo: float = 37.0,
    tm_hi: float = 42.0,
) -> Tuple[bool, str]:
    """Joint PBS2 gate (decision 8): length 11-17, GC 40-70 %, Tm 37-42 C, not
    repetitive.

    Length and Tm are coupled — a short PBS2 rarely reaches 37 C unless GC-rich —
    so the feasible set tends to collapse to ~14-17 nt; that is expected. Tm uses
    Biopython's nearest-neighbour model (imported lazily so the core import stays
    light).
    """
    if not (min_len <= len(pbs2) <= max_len):
        return False, f"length {len(pbs2)} outside [{min_len},{max_len}]"
    gc = _gc_content(pbs2)
    if not (gc_lo <= gc <= gc_hi):
        return False, f"GC {gc:.2f} outside [{gc_lo},{gc_hi}]"
    from Bio.SeqUtils import MeltingTemp as mt  # noqa: PLC0415

    tm = float(mt.Tm_NN(pbs2))
    if not (tm_lo <= tm <= tm_hi):
        return False, f"Tm {tm:.1f} outside [{tm_lo},{tm_hi}]"
    if repetitive:
        return False, "repetitive region"
    return True, "ok"


@dataclass
class RecodedTJPEResult:
    """Outcome of the recoded TJ-PE design for one target."""

    pegRNA_sequence: str
    spacer: str
    pbs2: str
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


def _captured_window(genome: str, *, nick_position: int, pbs_length: int) -> str:
    """Return the genomic window (5'->3', sense) immediately 5' of the nick that
    the new flap must capture. ``orient_pbs2`` then orients it per pegRNA strand."""
    start = max(0, nick_position - pbs_length)
    return genome[start:nick_position].upper()


def _frame_at(cds_mask: Sequence[bool], pos: int) -> int:
    """Reading-frame offset of genomic ``pos`` = (# CDS bases from CDS start to
    pos) % 3. Assumes a contiguous CDS in the local window (no intron between the
    CDS start and ``pos``), which holds for the short re-templated arm."""
    try:
        cds_start = list(cds_mask).index(True)
    except ValueError:
        return 0
    return max(0, pos - cds_start) % 3


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
    max_run: int = 8,
    max_rtt_variants: int = 20,
    # injected scorers (fakes in CI, real adapters on CX3)
    splice_scorer,
    full_structure_scorer,
    importance_scorer,
    dock_scorer,
    offtarget_scorer: Callable[[List[str]], Dict[str, Tuple[float, int, bool]]],
    cas9: str = "nCas9",
    allele_hunt: bool = False,
    het_snps: Optional[Sequence[dict]] = None,
    structure_threshold: float = 0.5,
    max_interface_drop: float = 0.15,
) -> "RecodedTJPEResult":
    """Ranked-pool-walk orchestrator (decisions 0-15).

    Builds a ranked candidate pool once; for each candidate runs the gates in
    cheap->expensive order (recode + INV-0 + PBS2 + ViennaRNA-class structure;
    then splice; then AlphaFold importance; then RiNALMo 2D structure; then
    Boltz-2 docking LAST). The first candidate passing every gate wins. The
    insert site ``site_s`` NEVER moves. Exhausting the pool raises
    :class:`NoValidPegRNA` with a per-candidate diagnostic.
    """
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

    for nick in nicks:
        # Re-templated span between S and the nick (the deletion-replacement arm).
        lo, hi = sorted((site_s, nick.nick_position))
        if lo != site_s:
            # Arm-placement contract: the pegRNA nick is at S, so the arm starts
            # at S. A nick 5' of S is not yet supported (would duplicate the arm).
            diagnostics.append(f"nick {nick.nick_position}: 5'-of-S placement unsupported")
            continue
        arm_genomic = genome[lo:hi].upper()
        frame = _frame_at(cds_mask, lo)
        # Recode the arm ONCE per nick (independent of the PBS-length variants).
        recode = recode_arm(arm_genomic, genome=arm_genomic, frame=frame, max_run=max_run)
        if recode is None:
            diagnostics.append(f"nick {nick.nick_position}: recode found no synonymous solution")
            continue
        if len(recode.recoded) != len(arm_genomic):
            diagnostics.append(f"nick {nick.nick_position}: recode changed arm length")
            continue
        # Independent synonymity check (feeds INV-0): recoded arm == same protein.
        if not _arm_protein_equal(arm_genomic, recode.recoded, frame):
            diagnostics.append(f"nick {nick.nick_position}: recode not protein-silent")
            continue

        for variant in variants:
            tried += 1
            # RTT tail = cargo + recoded arm, in the SENSE frame; orient it onto
            # the pegRNA per strand (Rv -> reverse_complement, Fw -> literal).
            sense_flap = cargo + recode.recoded
            rtt = orient_rtt_tail(sense_flap, pegrna_strand=pegrna_strand)
            # PBS2 = the captured-strand window 5' of the nick, oriented per strand.
            captured = _captured_window(
                genome, nick_position=nick.nick_position, pbs_length=variant.pbs_length
            )
            pbs2 = orient_pbs2(captured, pegrna_strand=pegrna_strand)
            ok, why = passes_pbs2_gate(pbs2, repetitive=False)
            if not ok:
                diagnostics.append(f"cand{tried}: PBS2 gate: {why}")
                continue
            # INV-0 on the predicted product, with a REAL protein-silence check
            # (recoded arm must translate to the same protein in-frame).
            reference = genome
            edited = genome[:site_s] + cargo + recode.recoded + genome[hi:]

            def _protein_check(_r, _e, _arm=arm_genomic, _rec=recode.recoded, _f=frame):
                return _arm_protein_equal(_arm, _rec, _f)

            try:
                assert_inv0(
                    reference=reference,
                    edited=edited,
                    site_s=site_s,
                    cargo=cargo,
                    recode_span=(lo, hi),
                    protein_check=_protein_check,
                    allowed_silent=None,
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
                    site_s=site_s,
                    hi=hi,
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
                pbs2=pbs2,
                rtt=rtt,
                pbs1=pbs1,
                polyu="UUU",
                scaffold=STANDARD_SGRNA_SCAFFOLD,
            )
            # --- expensive gates last; each is GUARDED so a scorer crash skips the
            # candidate (with a diagnostic) instead of aborting the whole walk. ---
            try:
                sp_ok, sp_why = passes_splice_gate(reference, edited, scorer=splice_scorer)
            except Exception as exc:  # noqa: BLE001
                diagnostics.append(f"cand{tried}: splice scorer error: {exc}")
                continue
            if not sp_ok:
                diagnostics.append(f"cand{tried}: splice: {sp_why}")
                continue

            # AlphaFold importance: map the recoded NT footprint to RESIDUE indices.
            res_lo = (lo - (list(cds_mask).index(True) if any(cds_mask) else 0)) // 3
            n_res = (len(recode.recoded) + 2) // 3
            footprint = list(range(res_lo, res_lo + max(1, n_res)))
            try:
                blocked = span_blocked(importance_scorer, footprint)
            except Exception as exc:  # noqa: BLE001
                diagnostics.append(f"cand{tried}: importance scorer error: {exc}")
                continue
            if blocked:
                diagnostics.append(f"cand{tried}: AlphaFold importance blocked a residue")
                continue

            regions = _regions_for(spacer, pbs2, rtt, pbs1)
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
            if pairing.cross_region_pairs > 0:
                diagnostics.append(f"cand{tried}: RiNALMo cross-region pairing")
                continue

            try:
                dock_ok, dock_why = passes_docking_gate(
                    peg, cas9, scorer=dock_scorer, max_drop=max_interface_drop
                )
            except Exception as exc:  # noqa: BLE001
                diagnostics.append(f"cand{tried}: docking scorer error: {exc}")
                continue
            if not dock_ok:
                diagnostics.append(f"cand{tried}: docking: {dock_why}")
                continue
            # First survivor wins.
            warnings: List[str] = []
            if recode.max_run_after > max_run:
                warnings.append("microhomology advisory: residual run above target")
            return RecodedTJPEResult(
                pegRNA_sequence=peg,
                spacer=spacer,
                pbs2=pbs2,
                rtt=rtt,
                pbs1=pbs1,
                polyu="UUU",
                ngrna_protospacer=nick.protospacer,
                site_s=site_s,
                nick_position=nick.nick_position,
                deletion_size=0,
                recode_edits=recode.n_edits,
                candidates_tried=tried,
                warnings=warnings,
                allele_specific_options=allele_options,
            )

    raise NoValidPegRNA(
        f"no candidate survived all gates after {tried} tried; "
        f"diagnostics: {'; '.join(diagnostics[:10])}"
    )


def _regions_for(spacer: str, pbs2: str, rtt: str, pbs1: str) -> Dict[str, Tuple[int, int]]:
    """Region offsets in the assembled pegRNA (authoritative order, scaffold excluded
    from cross-region scoring by name)."""
    scaf = len(STANDARD_SGRNA_SCAFFOLD)
    i0 = 0
    i1 = i0 + len(spacer)
    i2 = i1 + scaf
    i3 = i2 + len(pbs2)
    i4 = i3 + len(rtt)
    i5 = i4 + len(pbs1)
    return {
        "spacer": (i0, i1),
        "scaffold": (i1, i2),
        "pbs2": (i2, i3),
        "rtt": (i3, i4),
        "pbs1": (i4, i5),
    }
