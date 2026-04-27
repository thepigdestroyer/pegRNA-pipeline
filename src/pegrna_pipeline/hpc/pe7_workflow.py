"""Standard PE7 / PE3-strategy prime-editing workflow.

Entry point: :func:`run_pe7_workflow`.

PE7 is the canonical "single edit + nick on the opposite strand" prime-editing
strategy. Compared with TJ-PE / PA, it's structurally the simplest mode in
this pipeline:

* The caller supplies the **exact** insertion DNA (already codon-optimised).
  No first-N-amino-acid codon enumeration happens here — it's one variant in,
  one pegRNA out.
* PRIDICT 2.0 designs the pegRNA in a single batch call; no flap splitting,
  no rc_PBS2 derivation, no custom assembly.
* The PE3 nicking guide is picked from PRIDICT's ``--nicking`` output via the
  same opposite-strand selector TJ-PE uses
  (:func:`select_tjpe_ngrna`).

Hard limit: the insert must be ≤ ``PRIDICT_MAX_INSERT_BP`` (40 bp). Larger
edits should go through TJ-PE or PA — see ``tjpe_workflow`` / ``pa_workflow``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence

from pegrna_pipeline.hpc.data_types import (
    DEFAULT_SCORE_CELL_LINE,
    GeneRecord,
    NgRNADesign,
    PegRNADesign,
    VariantScore,
)
from pegrna_pipeline.hpc.tjpe_workflow import (
    rank_variants_by_pe2,
    select_tjpe_ngrna as select_pe3_ngrna,
)

#: Reliable upper bound for PRIDICT 2.0 insertion edits, per ground rule #3.
PRIDICT_MAX_INSERT_BP = 40


@dataclass
class PE7Result:
    """Final output of :func:`run_pe7_workflow`.

    Field names mirror :class:`TJPEResult` where the meaning is identical
    (``edit_position_in_mrna``, ``ngrna_*``) so report writers and the
    off-target scorer can treat the two modes uniformly.
    """

    winning_variant: Optional[VariantScore]
    full_insert_dna: str
    pegRNA_full: str
    spacer: str
    pbs: str
    rt: str
    pegrna_strand: str  # "Fw" | "Rv"
    pegrna_pe2_score: float  # PE2 score actually used for ranking
    pegrna_pe2_cell_line: str  # "k562" | "hek293t" — cell line that score came from
    ngrna_spacer: str
    ngrna_distance_to_edit: int
    ngrna_strand: str  # always opposite of pegrna_strand
    ngrna_deepcas9_score: float
    edit_position_in_mrna: int
    ngrna_nick_in_mrna: int
    all_evaluated_variants: List[VariantScore] = field(default_factory=list)


def _empty_result(
    *,
    gene: GeneRecord,
    insert_full_dna: str,
    variants: Sequence[VariantScore],
    score_cell_line: str,
    winning: Optional[VariantScore] = None,
    best_peg: Optional[PegRNADesign] = None,
) -> PE7Result:
    """Return a PE7Result with empty ngRNA fields for the degenerate paths."""
    from pegrna_pipeline.hpc.data_types import _score_attr

    edit_position = gene.stop_codon_start
    pe2_value = float(getattr(best_peg, _score_attr(score_cell_line))) if best_peg else 0.0
    return PE7Result(
        winning_variant=winning,
        full_insert_dna=insert_full_dna,
        pegRNA_full=best_peg.pegrna_full if best_peg else "",
        spacer=best_peg.spacer if best_peg else "",
        pbs=best_peg.pbs if best_peg else "",
        rt=best_peg.rt if best_peg else "",
        pegrna_strand=best_peg.target_strand if best_peg else "",
        pegrna_pe2_score=pe2_value,
        pegrna_pe2_cell_line=score_cell_line,
        ngrna_spacer="",
        ngrna_distance_to_edit=0,
        ngrna_strand="",
        ngrna_deepcas9_score=0.0,
        edit_position_in_mrna=edit_position,
        ngrna_nick_in_mrna=edit_position,
        all_evaluated_variants=list(variants),
    )


def run_pe7_workflow(
    *,
    gene: GeneRecord,
    insert_full_dna: str,
    variant_scores: Sequence[VariantScore],
    ngrna_distance_min: int = 40,
    ngrna_distance_max: int = 100,
    score_cell_line: str = DEFAULT_SCORE_CELL_LINE,
) -> PE7Result:
    """Run the PE7 / PE3-strategy workflow for one gene + insertion.

    Parameters
    ----------
    gene:
        Resolved :class:`GeneRecord` for the target locus.
    insert_full_dna:
        The DNA payload to insert at ``gene.stop_codon_start``. Must be
        non-empty uppercase ACGT with ``len <= PRIDICT_MAX_INSERT_BP``.
    variant_scores:
        Output of PRIDICT for **one** input editseq (typically a single
        :class:`VariantScore`). The list form is preserved so the driver can
        feed exactly one entry without an awkward unwrap.
    ngrna_distance_min, ngrna_distance_max:
        PE3 nicking-window bounds on ``|distance_to_edit|``.

    Returns
    -------
    PE7Result
        Contains the winning pegRNA + ngRNA pair (or empty fields if either
        could not be selected — e.g. PRIDICT returned no candidates passing
        the spacer sanity filters).

    Raises
    ------
    ValueError
        If ``insert_full_dna`` is empty or longer than
        :data:`PRIDICT_MAX_INSERT_BP`.
    """
    if not insert_full_dna:
        raise ValueError("insert_full_dna must be non-empty for PE7 mode")
    if len(insert_full_dna) > PRIDICT_MAX_INSERT_BP:
        raise ValueError(
            f"PE7 mode supports inserts ≤ {PRIDICT_MAX_INSERT_BP} bp "
            f"(got {len(insert_full_dna)} bp); use TJ-PE or PA for larger edits"
        )

    # Propagate score_cell_line to every variant so best_pegRNA / top_score
    # use the requested column.
    variants = [_with_score_cell_line(v, score_cell_line) for v in variant_scores]
    ranked = rank_variants_by_pe2(variants)
    winning: Optional[VariantScore] = ranked[0] if ranked else None
    best_peg: Optional[PegRNADesign] = winning.best_pegRNA if winning else None

    if winning is None or best_peg is None:
        return _empty_result(
            gene=gene,
            insert_full_dna=insert_full_dna,
            variants=variants,
            score_cell_line=score_cell_line,
            winning=winning,
            best_peg=None,
        )

    best_ng: Optional[NgRNADesign] = select_pe3_ngrna(
        winning.ngRNA_candidates,
        min_distance=ngrna_distance_min,
        max_distance=ngrna_distance_max,
        pegrna_strand=best_peg.target_strand,
    )

    if best_ng is None:
        return _empty_result(
            gene=gene,
            insert_full_dna=insert_full_dna,
            variants=variants,
            score_cell_line=score_cell_line,
            winning=winning,
            best_peg=best_peg,
        )

    edit_position_in_mrna = gene.stop_codon_start
    ngrna_nick_in_mrna = edit_position_in_mrna + best_ng.distance_to_edit
    from pegrna_pipeline.hpc.data_types import _score_attr

    pe2_value = float(getattr(best_peg, _score_attr(score_cell_line)))

    return PE7Result(
        winning_variant=winning,
        full_insert_dna=insert_full_dna,
        pegRNA_full=best_peg.pegrna_full,
        spacer=best_peg.spacer,
        pbs=best_peg.pbs,
        rt=best_peg.rt,
        pegrna_strand=best_peg.target_strand,
        pegrna_pe2_score=pe2_value,
        pegrna_pe2_cell_line=score_cell_line,
        ngrna_spacer=best_ng.spacer,
        ngrna_distance_to_edit=best_ng.distance_to_edit,
        ngrna_strand=best_ng.target_strand,
        ngrna_deepcas9_score=best_ng.deepcas9_score,
        edit_position_in_mrna=edit_position_in_mrna,
        ngrna_nick_in_mrna=ngrna_nick_in_mrna,
        all_evaluated_variants=variants,
    )


def _with_score_cell_line(variant: VariantScore, score_cell_line: str) -> VariantScore:
    """Return a copy of ``variant`` with ``score_cell_line`` set.

    Builds a fresh VariantScore so callers' originals stay untouched (the
    field is mutable on the dataclass since ``frozen=False``).
    """
    return VariantScore(
        sequence_name=variant.sequence_name,
        codon_variant_dna=variant.codon_variant_dna,
        pegRNA_designs=variant.pegRNA_designs,
        ngRNA_candidates=variant.ngRNA_candidates,
        score_cell_line=score_cell_line,
    )
