"""Tethered-junction prime editing (TJ-PE) workflow.

Entry point: :func:`run_tjpe_workflow`. See Plan 2 v1.1 for the pipeline
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
from pegrna_pipeline.hpc.pridict_parser import (
    parse_ngrna_csv,
    parse_pegrna_csv,
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


def _spacer_passes_sanity(spacer: str) -> bool:
    """Reject unusable spacers on their own sequence content.

    * ``TTTT`` is a RNA-pol III terminator — guide will truncate.
    * A homopolymer run of >= 5 bp kills expression.
    * GC content outside [0.30, 0.75] gives poor DNA binding.
    """
    s = spacer.upper()
    if "TTTT" in s:
        return False
    # Homopolymer >= 5
    run = 1
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            run += 1
            if run >= 5:
                return False
        else:
            run = 1
    if not s:
        return False
    gc = (s.count("G") + s.count("C")) / len(s)
    if gc < 0.30 or gc > 0.75:
        return False
    return True


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
        if not _spacer_passes_sanity(c.spacer):
            continue
        survivors.append(c)

    if not survivors:
        return None
    return max(survivors, key=lambda c: c.deepcas9_score)


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
    scaffold: str = STANDARD_SGRNA_SCAFFOLD,
) -> str:
    """Return the full TJ-PE pegRNA DNA string.

    Layout (5' -> 3'):

        spacer + scaffold + insert_dna + pbs1 + rc_pbs2

    All four sequence components (plus the scaffold) must be non-empty
    uppercase ACGT. Any lowercase letter, ``U``, or ambiguity code raises
    :class:`ValueError`.
    """
    for label, value in (
        ("spacer", spacer),
        ("insert_dna", insert_dna),
        ("pbs1", pbs1),
        ("rc_pbs2", rc_pbs2),
        ("scaffold", scaffold),
    ):
        _check_acgt(label, value)
    return spacer + scaffold + insert_dna + pbs1 + rc_pbs2


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
) -> TJPEResult:
    """Run the TJ-PE workflow end-to-end for one gene + insertion.

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

    # -------- assemble final pegRNA --------
    tjpe_pegrna = assemble_tj_pegrna(
        spacer=best_peg.spacer,
        insert_dna=insert_full_dna,
        pbs1=best_peg.pbs,
        rc_pbs2=rc_pbs2,
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
    )
