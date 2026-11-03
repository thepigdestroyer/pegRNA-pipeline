"""Tethered-junction prime editing (TJ-PE) workflow — shared primitives + legacy orchestrator.

.. note:: SINGLE CANONICAL DESIGN WORKFLOW (Stage H, ADR 0021).
   :func:`pegrna_pipeline.hpc.recoded_tjpe_workflow.run_recoded_tjpe` is the ONE canonical
   tj-PE design workflow used by the CLI + web app on the DSI box (auto-spacer via PRIDICT2,
   the full gate cascade, the recode arm, and the het-edit pair). :func:`run_tjpe_workflow`
   below is the LEGACY variant-ranking design orchestrator (the HPC ``driver``/PBS path); it
   is retained for backward compatibility and for the shared primitives this module exports
   (``STANDARD_SGRNA_SCAFFOLD``, ``TJPEResult``, :func:`rank_variants_by_pe2`,
   :func:`select_tjpe_ngrna`, :func:`derive_rc_pbs2`, :func:`assemble_tj_pegrna_v2`). Both
   paths assemble through :func:`assemble_tj_pegrna_v2` in the authoritative
   ``Spacer.Scaffold.rcPBS2.RTT.PBS1.PolyU`` layout (ADR 0014/0020), so they cannot drift.

Legacy entry point: :func:`run_tjpe_workflow`. See Plan 2 v1.1 for the pipeline
contract. Tasks here:

* :func:`rank_variants_by_pe2` — T6
* :func:`select_tjpe_ngrna`   — T7
* :func:`derive_rc_pbs2`      — T8
* :func:`assemble_tj_pegrna`  — T9
* :class:`TJPEResult`         — T10
* :func:`run_tjpe_workflow`   — T10
"""

from __future__ import annotations

import re
import warnings as _warnings  # aliased: run_tjpe_workflow has a local `warnings` list
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence, Union

from Bio.Seq import Seq

from pegrna_pipeline.hpc.data_types import (
    DEFAULT_SCORE_CELL_LINE,
    GeneRecord,
    NgRNADesign,
    PegRNADesign,
    VariantScore,
    _score_attr,
)
from pegrna_pipeline.hpc.guide_sanity import spacer_passes_sanity
from pegrna_pipeline.hpc.pridict_parser import (
    parse_ngrna_csv,
    parse_pegrna_csv,
)
from pegrna_pipeline.hpc.scarless import (
    assert_inv0,
    build_restored_segment,
    predict_edited_allele,
)

OPPOSITE_STRAND = {"Fw": "Rv", "Rv": "Fw"}


#: Standard 76-nt sgRNA scaffold sequence (tracrRNA fusion used in
#: prime-editing constructs).
STANDARD_SGRNA_SCAFFOLD = (
    "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
)


_ACGT_ONLY = re.compile(r"^[ACGT]+$")


# ------------------------------------------------------------------ T6


def rank_variants_by_pe2(variants: List[VariantScore]) -> List[VariantScore]:
    """Return a new list, sorted by ``top_score`` descending, stably.

    Stable sort ensures callers can rely on insertion order for ties.
    """
    return sorted(variants, key=lambda v: v.top_score, reverse=True)


# ------------------------------------------------------------------ T7


def select_tjpe_ngrna(
    candidates: List[NgRNADesign],
    *,
    min_distance: int = 40,
    max_distance: int = 100,
    pegrna_strand: Optional[str] = None,
) -> Optional[NgRNADesign]:
    """Pick the best (highest ``deepcas9_score``) surviving ngRNA.

    Parameters
    ----------
    candidates:
        List of :class:`NgRNADesign` rows, typically parsed from a PRIDICT
        nicking-guides CSV.
    min_distance, max_distance:
        Inclusive bounds on ``|distance_to_edit|``.
    pegrna_strand:
        If given (``"Fw"`` or ``"Rv"``), only keep opposite-strand
        candidates: ``{"Fw": "Rv", "Rv": "Fw"}``. Raises :class:`KeyError`
        for any other value — unknown strands are always a bug, not a
        silent filter.

    Returns
    -------
    The candidate with the maximum ``deepcas9_score`` that passes every
    filter, or ``None`` if no candidate survives.
    """
    required_strand: Optional[str] = None
    if pegrna_strand is not None:
        required_strand = OPPOSITE_STRAND[pegrna_strand]  # raises KeyError if unknown

    survivors: List[NgRNADesign] = []
    for c in candidates:
        d = abs(c.distance_to_edit)
        if d < min_distance or d > max_distance:
            continue
        if required_strand is not None and c.target_strand != required_strand:
            continue
        if not spacer_passes_sanity(c.spacer):  # hard rules only; see guide_sanity
            continue
        survivors.append(c)

    if not survivors:
        return None
    return max(survivors, key=lambda c: c.deepcas9_score)


# ------------------------------------------------- scarless helpers (Tasks 6/7)


def require_b_at_s(*, distance_to_edit: int, pegrna_strand: str) -> bool:
    """Scarless strand rule: a reverse-strand pegRNA requires the nicking-sgRNA
    nick to fall **exactly at the edit site S** (signed ``distance_to_edit == 0``)
    so no endogenous base is lost. For a forward-strand pegRNA this constraint
    does not apply (returns ``True``).
    """
    if pegrna_strand == "Rv":
        return distance_to_edit == 0
    if pegrna_strand == "Fw":
        return True
    raise ValueError(f"pegrna_strand must be 'Fw'/'Rv'; got {pegrna_strand!r}")


def combine_two_end_scores(*, insert_score: float, pbs2_score: float) -> float:
    """Bottleneck combination of the two-end PRIDICT decomposition: a design is
    only as good as its weaker end (insert/initiation vs PBS2/second-strand)."""
    return min(float(insert_score), float(pbs2_score))


# ------------------------------------------------------------------ T8


def derive_rc_pbs2(
    top_strand_sequence: str,
    *,
    ngrna_nick_position: int,
    pbs2_length: int = 17,
) -> str:
    """Return the reverse complement of the PBS2 window.

    PBS2 is the primer-binding site used by the **downstream** pegRNA in
    a TJ-PE construct. It's extracted from the top strand immediately
    5' of the ngRNA nick, then reverse-complemented so it anneals to the
    opposite (nicked) strand.

    Parameters
    ----------
    top_strand_sequence:
        mRNA-sense sequence (the same strand as ``GeneRecord.mrna_sequence``).
    ngrna_nick_position:
        0-indexed position on ``top_strand_sequence`` where the ngRNA cuts.
    pbs2_length:
        Number of nt upstream of the nick to use as PBS2 (default 17).

    Raises
    ------
    ValueError
        If the window ``[nick - pbs2_length, nick)`` falls outside the
        sequence.
    """
    start = ngrna_nick_position - pbs2_length
    if start < 0:
        raise ValueError(
            f"PBS2 window start {start} < 0 "
            f"(nick={ngrna_nick_position}, pbs2_length={pbs2_length})"
        )
    if ngrna_nick_position > len(top_strand_sequence):
        raise ValueError(
            f"ngrna_nick_position {ngrna_nick_position} exceeds sequence length "
            f"{len(top_strand_sequence)}"
        )
    window = top_strand_sequence[start:ngrna_nick_position]
    return str(Seq(window).reverse_complement())


# ------------------------------------------------------------------ T9


def _check_acgt(label: str, value: str) -> None:
    if not value or not _ACGT_ONLY.match(value):
        raise ValueError(f"{label} must be non-empty uppercase ACGT; got {value!r}")


def assemble_tj_pegrna(
    *,
    spacer: str,
    insert_dna: str,
    pbs1: str,
    rc_pbs2: str,
    restored_segment: str = "",
    scaffold: str = STANDARD_SGRNA_SCAFFOLD,
) -> str:
    """Return the full TJ-PE pegRNA DNA string.

    Current layout (5' -> 3'):

        spacer + scaffold + insert_dna + restored_segment + pbs1 + rc_pbs2

    .. note:: AUTHORITATIVE tjPE pegRNA STRUCTURE (locked 2026-06-18, B1).
       The recoded-redesign target order is::

           Spacer · Scaffold · PBS2 · RTT(Insert + Deletion-Replacement-Arm) · PBS1 · PolyU(2-3 U)

       i.e. **PBS2 occupies the Homology-Arm slot** (HA is *replaced* by PBS2 in
       tjPE) immediately 3' of the scaffold and 5' of the RTT; PBS1 is INTERNAL;
       the construct ends in a **PolyU (2-3 U) tail that binds the La protein of
       the PE7 editor**. Mechanism: PBS1 (3'-proximal) primes first-strand
       synthesis; the RT copies the RTT; PBS2 (5' of the extension) is copied
       last and becomes the 3' jumping end for the template jump. This function
       is rewritten to that signature in
       ``docs/superpowers/plans/2026-06-18-recoded-tjpe-pegrna-redesign.md``
       Task 11b. The historical 3'-terminal ``rc_pbs2`` arg below predates that
       decision; the PBS2-slot literal-vs-reverse-complement convention is
       asserted against the predicted product in the recoded workflow (Task 13).

    The optional ``restored_segment`` carries the re-templated A→B genomic span
    used by the **scarless** mode (zero endogenous deletion); when empty the
    layout is the standard deleting tjPE construct, so existing callers are
    unaffected. All non-empty sequence components (plus the scaffold) must be
    uppercase ACGT. Any lowercase letter, ``U``, or ambiguity code raises
    :class:`ValueError`.
    """
    checks = [
        ("spacer", spacer),
        ("insert_dna", insert_dna),
        ("pbs1", pbs1),
        ("rc_pbs2", rc_pbs2),
        ("scaffold", scaffold),
    ]
    if restored_segment:
        checks.append(("restored_segment", restored_segment))
    for label, value in checks:
        _check_acgt(label, value)
    return spacer + scaffold + insert_dna + restored_segment + pbs1 + rc_pbs2


def assemble_tj_pegrna_v2(
    *,
    spacer: str,
    pbs2: str,
    rtt: str,
    pbs1: str,
    polyu: str = "UUU",
    scaffold: str = STANDARD_SGRNA_SCAFFOLD,
) -> str:
    """Assemble the full tjPE pegRNA in the AUTHORITATIVE 5'->3' order (decision 0).

    Layout::

        spacer + scaffold + PBS2 + RTT + PBS1 + PolyU

    * ``pbs2`` occupies the Homology-Arm slot (HA is REPLACED by PBS2 in tjPE),
      immediately 3' of the scaffold and 5' of the RTT. Its RT copy becomes the
      3' jumping end that hybridises to the second nick (template jump).
    * ``rtt`` = insert + deletion-replacement-arm (recoded; passed pre-assembled).
    * ``pbs1`` primes first-strand synthesis; it is INTERNAL (5' of PolyU).
    * ``polyu`` = 2-3 uridines, the 3' free tail bound by the La protein of the
      PE7 editor. Only ``U`` permitted.

    DNA components (spacer/pbs2/rtt/pbs1/scaffold) must be uppercase ACGT; the
    PolyU tail must be 2-3 ``U``. This is a SEPARATE function from the legacy
    :func:`assemble_tj_pegrna` so the old scarless workflow + its tests are
    unaffected (see Task 11b deviation note).

    .. note:: The ``pbs2`` orientation (literal vs reverse-complement) at this 5'
       slot is asserted against the predicted product by the recoded workflow
       (Task 13), not here.
    """
    for label, value in (
        ("spacer", spacer),
        ("scaffold", scaffold),
        ("pbs2", pbs2),
        ("rtt", rtt),
        ("pbs1", pbs1),
    ):
        _check_acgt(label, value)
    if not polyu or set(polyu.upper()) - {"U"}:
        raise ValueError(f"polyu must be 2-3 uridines (U only); got {polyu!r}")
    if not (2 <= len(polyu) <= 3):
        raise ValueError(f"polyu length must be 2-3; got {len(polyu)}")
    return spacer + scaffold + pbs2 + rtt + pbs1 + polyu


# ------------------------------------------------------------------ T10


@dataclass
class TJPEResult:
    """Final output of :func:`run_tjpe_workflow` — v1.1 field names.

    Field names match the corrections supplement exactly. Renames vs v1.0:

    * ``first_nick_position``  -> ``edit_position_in_mrna``
    * ``second_nick_position`` -> ``ngrna_nick_in_mrna``
    * ``ngrna_distance``       -> ``ngrna_distance_to_edit``
    """

    winning_variant: Optional[VariantScore]
    full_insert_dna: str
    tjpe_pegRNA_sequence: str
    spacer: str
    pbs1: str
    rc_pbs2: str
    rc_pbs2_length: int
    ngrna_spacer: str
    ngrna_distance_to_edit: int
    ngrna_deepcas9_score: float
    pegrna_pe2_score: float  # PE2 score actually used for ranking
    pegrna_pe2_cell_line: str  # "k562" | "hek293t" — cell line that score came from
    edit_position_in_mrna: int
    ngrna_nick_in_mrna: int
    # Informational PE3 composite — pegRNA PE2 score × ngRNA DeepCas9 score
    # / 100. Reported for transparency; NOT used for ranking.
    joint_pe3_score: float = 0.0
    all_evaluated_variants: List[VariantScore] = field(default_factory=list)
    # ---- scarless mode fields (defaults keep back-compat for deleting tjPE) ----
    scarless: bool = False
    restored_segment: str = ""
    silent_pam_edit: Optional[dict] = None
    predicted_edited_allele: str = ""
    deletion_size: int = 0
    insert_end_score: float = 0.0
    pbs2_end_score: float = 0.0
    combined_score: float = 0.0
    microhomology_flag: bool = False
    recut_flag: bool = False
    structure_norm_mfe: float = 0.0
    warnings: List[str] = field(default_factory=list)


class ScarlessNoValidDesign(RuntimeError):
    """Raised when scarless mode finds no candidate passing all hard gates.

    Per the adversarial-review risk register (R4), scarless mode must **fail
    loudly with a reason** rather than silently return an empty/garbage design.
    """


def _load_variants_from_fixtures(fixtures_dir: Path) -> List[VariantScore]:
    """Discover ``*_pegRNA_Pridict_full.csv`` + ``*_nicking_guides.csv`` pairs
    in ``fixtures_dir`` and build :class:`VariantScore` objects.

    Variant name = the common prefix (filename up to the ``_pegRNA_...``
    suffix). The ngRNA file must exist with the same prefix; a missing
    partner is a :class:`FileNotFoundError`.
    """
    fixtures_dir = Path(fixtures_dir)
    variants: List[VariantScore] = []
    for peg_path in sorted(fixtures_dir.glob("*_pegRNA_Pridict_full.csv")):
        prefix = peg_path.name.replace("_pegRNA_Pridict_full.csv", "")
        ng_path = fixtures_dir / f"{prefix}_nicking_guides.csv"
        if not ng_path.exists():
            raise FileNotFoundError(
                f"no matching nicking guides CSV for {peg_path.name} " f"(expected {ng_path.name})"
            )
        pegrnas = parse_pegrna_csv(peg_path)
        ngrnas = parse_ngrna_csv(ng_path)
        variants.append(
            VariantScore(
                sequence_name=prefix,
                codon_variant_dna="",  # insert_full_dna is passed separately
                pegRNA_designs=pegrnas,
                ngRNA_candidates=ngrnas,
            )
        )
    return variants


def run_tjpe_workflow(
    *,
    gene: GeneRecord,
    insert_full_dna: str,
    variant_fixtures_dir: Optional[Union[str, Path]] = None,
    variant_scores: Optional[Sequence[VariantScore]] = None,
    enumeration_window_bp: int = 24,  # noqa: ARG001  # Plan 3 forward-compat
    ngrna_distance_min: int = 40,
    ngrna_distance_max: int = 100,
    pbs2_length: int = 17,
    score_cell_line: str = DEFAULT_SCORE_CELL_LINE,
    scarless: bool = False,
) -> TJPEResult:
    """Run the LEGACY TJ-PE design workflow end-to-end for one gene + insertion.

    .. deprecated:: Stage H (ADR 0021)
       Use :func:`pegrna_pipeline.hpc.recoded_tjpe_workflow.run_recoded_tjpe`, the single
       canonical design workflow (auto-spacer + full gate cascade, used by the CLI/web).
       This orchestrator is retained for the HPC ``driver`` path and the shared primitives.

    Exactly one of ``variant_fixtures_dir`` or ``variant_scores`` must be
    supplied. The former is used by the current test / fixture-driven
    path; the latter is the programmatic entry point Plan 3 will wire up
    (variants enumerated in-memory by the codon expander).

    Nick arithmetic follows Plan 2 v1.1 exactly::

        edit_position_in_mrna = gene.stop_codon_start
        ngrna_nick_in_mrna    = edit_position_in_mrna + ng.distance_to_edit
        rc_pbs2               = derive_rc_pbs2(
                                    gene.mrna_sequence,
                                    ngrna_nick_position=ngrna_nick_in_mrna,
                                    pbs2_length=pbs2_length,
                                )
    """
    _warnings.warn(
        "run_tjpe_workflow is the LEGACY tj-PE design orchestrator; the single canonical "
        "design workflow is recoded_tjpe_workflow.run_recoded_tjpe (auto-spacer + full gate "
        "cascade, used by the CLI/web). This path is retained for the HPC driver + shared "
        "primitives only (ADR 0021).",
        DeprecationWarning,
        stacklevel=2,
    )
    # -------- input source contract --------
    if (variant_fixtures_dir is None) == (variant_scores is None):
        raise ValueError("exactly one of variant_fixtures_dir or variant_scores must be provided")

    if variant_fixtures_dir is not None:
        variants = _load_variants_from_fixtures(Path(variant_fixtures_dir))
    else:
        variants = list(variant_scores or [])

    # Stamp the requested cell line on every variant so best_pegRNA /
    # top_score select the right PE2 column.
    for v in variants:
        v.score_cell_line = score_cell_line

    # -------- rank & pick best pegRNA --------
    ranked = rank_variants_by_pe2(variants)
    winning: Optional[VariantScore] = ranked[0] if ranked else None
    best_peg: Optional[PegRNADesign] = winning.best_pegRNA if winning else None

    edit_position_in_mrna = gene.stop_codon_start

    if winning is None or best_peg is None:
        if scarless:
            raise ScarlessNoValidDesign(
                f"no valid scarless design for {gene.gene_name}: no ranked pegRNA "
                f"variant available (0 of {len(variants)} variants usable)"
            )
        # Degenerate case: caller asked us to design with no variants.  We
        # return a TJPEResult with zeroed strings so downstream diagnostics
        # can still inspect all_evaluated_variants.
        return TJPEResult(
            winning_variant=winning,
            full_insert_dna=insert_full_dna,
            tjpe_pegRNA_sequence="",
            spacer="",
            pbs1="",
            rc_pbs2="",
            rc_pbs2_length=0,
            ngrna_spacer="",
            ngrna_distance_to_edit=0,
            ngrna_deepcas9_score=0.0,
            pegrna_pe2_score=0.0,
            pegrna_pe2_cell_line=score_cell_line,
            edit_position_in_mrna=edit_position_in_mrna,
            ngrna_nick_in_mrna=edit_position_in_mrna,
            all_evaluated_variants=list(variants),
        )

    # -------- pick best ngRNA (opposite strand) --------
    best_ng = select_tjpe_ngrna(
        winning.ngRNA_candidates,
        min_distance=ngrna_distance_min,
        max_distance=ngrna_distance_max,
        pegrna_strand=best_peg.target_strand,
    )

    if best_ng is None:
        if scarless:
            raise ScarlessNoValidDesign(
                f"no valid scarless design for {gene.gene_name}: no safe nicking "
                f"sgRNA in the {ngrna_distance_min}-{ngrna_distance_max} bp window "
                f"(strand-opposite of {best_peg.target_strand})"
            )
        # No usable ngRNA for this variant — surface via empty ngRNA fields
        # so callers can decide to fail vs retry.  A full-on exception would
        # be wrong here because the caller may well want to iterate through
        # fallback variants.
        return TJPEResult(
            winning_variant=winning,
            full_insert_dna=insert_full_dna,
            tjpe_pegRNA_sequence="",
            spacer=best_peg.spacer,
            pbs1=best_peg.pbs,
            rc_pbs2="",
            rc_pbs2_length=0,
            ngrna_spacer="",
            ngrna_distance_to_edit=0,
            ngrna_deepcas9_score=0.0,
            pegrna_pe2_score=float(getattr(best_peg, _score_attr(score_cell_line))),
            pegrna_pe2_cell_line=score_cell_line,
            edit_position_in_mrna=edit_position_in_mrna,
            ngrna_nick_in_mrna=edit_position_in_mrna,
            all_evaluated_variants=list(variants),
        )

    # -------- derive rc_PBS2 from mRNA around ngRNA nick --------
    ngrna_nick_in_mrna = edit_position_in_mrna + best_ng.distance_to_edit
    rc_pbs2 = derive_rc_pbs2(
        gene.mrna_sequence,
        ngrna_nick_position=ngrna_nick_in_mrna,
        pbs2_length=pbs2_length,
    )

    # -------- scarless construction (zero endogenous deletion) --------
    restored_segment = ""
    predicted_allele = ""
    deletion_size = 0
    warnings: List[str] = []
    if scarless:
        # Reverse-strand pegRNA must have its nicking nick exactly at S.
        if not require_b_at_s(
            distance_to_edit=best_ng.distance_to_edit,
            pegrna_strand=best_peg.target_strand,
        ):
            raise ScarlessNoValidDesign(
                f"no valid scarless design for {gene.gene_name}: reverse-strand "
                f"pegRNA requires the nicking nick exactly at S (distance 0), "
                f"got {best_ng.distance_to_edit}"
            )
        a_nick = edit_position_in_mrna
        b_nick = ngrna_nick_in_mrna
        restored_segment = build_restored_segment(
            gene.mrna_sequence,
            a_nick=a_nick,
            b_nick=b_nick,
            pegrna_strand=best_peg.target_strand,
        )
        predicted_allele = predict_edited_allele(
            gene.mrna_sequence,
            site_s=edit_position_in_mrna,
            cargo=insert_full_dna,
            a_nick=a_nick,
            b_nick=b_nick,
            pegrna_strand=best_peg.target_strand,
            silent_pam_edit=None,
        )
        # INV-0 HARD GATE — zero endogenous deletion (raises on failure).
        assert_inv0(
            reference=gene.mrna_sequence,
            edited=predicted_allele,
            site_s=edit_position_in_mrna,
            cargo=insert_full_dna,
            allowed_silent=None,
        )
        # Advisory warnings (A-only: surfaced, never blocking).
        warnings.append(
            "microhomology: restored segment is 100% genome-homologous " "(see spec §10 R1)"
        )
        if abs(best_ng.distance_to_edit) < 70:
            warnings.append(
                f"distance {abs(best_ng.distance_to_edit)} bp < 70 bp validated "
                "minimum (R3, untested)"
            )

    # -------- assemble final pegRNA (authoritative layout, ADR 0014 / decision 0) --------
    # The legacy 3'-terminal-rc_pbs2 assemble_tj_pegrna (v1) predates the locked
    # tjPE structure; converge onto assemble_tj_pegrna_v2 so this path emits the same
    # authoritative order as the recoded path: Spacer.Scaffold.rcPBS2.RTT.PBS1.PolyU.
    # rc_pbs2 (already RC(TS) via derive_rc_pbs2) fills the HA slot; the RTT is the
    # insert plus the scarless deletion-replacement arm.
    tjpe_pegrna = assemble_tj_pegrna_v2(
        spacer=best_peg.spacer,
        pbs2=rc_pbs2,
        rtt=insert_full_dna + restored_segment,
        pbs1=best_peg.pbs,
    )

    pe2_value = float(getattr(best_peg, _score_attr(score_cell_line)))
    return TJPEResult(
        winning_variant=winning,
        full_insert_dna=insert_full_dna,
        tjpe_pegRNA_sequence=tjpe_pegrna,
        spacer=best_peg.spacer,
        pbs1=best_peg.pbs,
        rc_pbs2=rc_pbs2,
        rc_pbs2_length=len(rc_pbs2),
        ngrna_spacer=best_ng.spacer,
        ngrna_distance_to_edit=best_ng.distance_to_edit,
        ngrna_deepcas9_score=best_ng.deepcas9_score,
        pegrna_pe2_score=pe2_value,
        pegrna_pe2_cell_line=score_cell_line,
        edit_position_in_mrna=edit_position_in_mrna,
        ngrna_nick_in_mrna=ngrna_nick_in_mrna,
        joint_pe3_score=pe2_value * float(best_ng.deepcas9_score) / 100.0,
        all_evaluated_variants=list(variants),
        scarless=scarless,
        restored_segment=restored_segment,
        predicted_edited_allele=predicted_allele,
        deletion_size=deletion_size,
        microhomology_flag=bool(scarless),
        recut_flag=False,
        warnings=warnings,
    )
