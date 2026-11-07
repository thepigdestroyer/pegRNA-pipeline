"""Standard PE7 / PE3-strategy prime-editing workflow.

Entry point: :func:`run_pe7_workflow`.

PE7 is the canonical "single edit + nick on the opposite strand" prime-editing
strategy. Compared with TJ-PE / PA, it's structurally the simplest mode in
this pipeline:

* The caller supplies the **exact** insertion DNA (already codon-optimised).
  No first-N-amino-acid codon enumeration happens here — it's one variant in,
  many ranked (pegRNA, ngRNA) pairs out.
* PRIDICT 2.0 designs the pegRNA in a single batch call; no flap splitting,
  no rc_PBS2 derivation, no custom assembly.
* The PE3 nicking guide is picked from PRIDICT's ``--nicking`` output via the
  same opposite-strand selector TJ-PE uses
  (:func:`select_tjpe_ngrna`).
* We return up to ``top_n`` ranked pegRNA + ngRNA pairs (default 20) so the
  user can inspect alternatives — not just the rank-1 winner.

Hard limit: the insert must be ≤ ``PRIDICT_MAX_INSERT_BP`` (40 bp). Larger
edits should go through TJ-PE (``tjpe_workflow``) or Prime Assembly (``pegrna design-pa``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence

from pegrna_pipeline.hpc.data_types import (
    DEFAULT_SCORE_CELL_LINE,
    GeneRecord,
    VariantScore,
    _score_attr,
)
from pegrna_pipeline.hpc.tjpe_workflow import (
    rank_variants_by_pe2,
)
from pegrna_pipeline.hpc.tjpe_workflow import (
    select_tjpe_ngrna as select_pe3_ngrna,
)

#: Reliable upper bound for PRIDICT 2.0 insertion edits, per ground rule #3.
PRIDICT_MAX_INSERT_BP = 40

#: Default top-N pairs to return.
DEFAULT_TOP_N = 20


@dataclass(frozen=True)
class PE7DesignPair:
    """One ranked (pegRNA, ngRNA) PE7 design.

    ``rank`` is 1-based — rank 1 is the highest-scoring pair by
    ``pegrna_pe2_score`` (with K562 as the default cell-line column,
    overridable via ``score_cell_line``).
    """

    rank: int
    variant_sequence_name: str
    pegRNA_full: str
    spacer: str
    pbs: str
    rt: str
    pegrna_strand: str  # "Fw" | "Rv"
    pegrna_pe2_score: float
    pegrna_pe2_cell_line: str
    ngrna_spacer: str
    ngrna_distance_to_edit: int
    ngrna_strand: str
    ngrna_deepcas9_score: float
    ngrna_nick_in_mrna: int
    # Informational composite — pegRNA PE2 score × ngRNA DeepCas9 score / 100.
    joint_pe3_score: float


@dataclass
class PE7Result:
    """Final output of :func:`run_pe7_workflow`.

    Holds up to ``top_n`` ranked (pegRNA, ngRNA) pairs in ``top_pairs``.
    The legacy top-level fields (``spacer``, ``ngrna_spacer``, etc.) mirror
    ``top_pairs[0]`` for backward compatibility with callers that only
    care about the rank-1 design.
    """

    winning_variant: Optional[VariantScore]
    full_insert_dna: str
    pegRNA_full: str
    spacer: str
    pbs: str
    rt: str
    pegrna_strand: str
    pegrna_pe2_score: float
    pegrna_pe2_cell_line: str
    ngrna_spacer: str
    ngrna_distance_to_edit: int
    ngrna_strand: str
    ngrna_deepcas9_score: float
    edit_position_in_mrna: int
    ngrna_nick_in_mrna: int
    joint_pe3_score: float = 0.0
    top_pairs: List[PE7DesignPair] = field(default_factory=list)
    all_evaluated_variants: List[VariantScore] = field(default_factory=list)


def _empty_result(
    *,
    gene: GeneRecord,
    insert_full_dna: str,
    variants: Sequence[VariantScore],
    score_cell_line: str,
) -> PE7Result:
    """Return a PE7Result with empty fields for the no-design path."""
    edit_position = gene.stop_codon_start
    return PE7Result(
        winning_variant=None,
        full_insert_dna=insert_full_dna,
        pegRNA_full="",
        spacer="",
        pbs="",
        rt="",
        pegrna_strand="",
        pegrna_pe2_score=0.0,
        pegrna_pe2_cell_line=score_cell_line,
        ngrna_spacer="",
        ngrna_distance_to_edit=0,
        ngrna_strand="",
        ngrna_deepcas9_score=0.0,
        edit_position_in_mrna=edit_position,
        ngrna_nick_in_mrna=edit_position,
        joint_pe3_score=0.0,
        top_pairs=[],
        all_evaluated_variants=list(variants),
    )


def _build_top_pairs(
    *,
    variants: Sequence[VariantScore],
    score_cell_line: str,
    ngrna_distance_min: int,
    ngrna_distance_max: int,
    edit_position: int,
    top_n: int,
) -> List[PE7DesignPair]:
    """Walk every pegRNA across every variant in score order; pair each with
    its best valid PE3 ngRNA; keep the first ``top_n`` that have a partner.
    """
    score_attr = _score_attr(score_cell_line)
    pairs: List[PE7DesignPair] = []

    # Flatten (variant, peg) across all variants and sort by peg score desc.
    flat: List[tuple] = []
    for v in variants:
        for peg in v.pegRNA_designs:
            flat.append((v, peg, float(getattr(peg, score_attr))))
    flat.sort(key=lambda t: t[2], reverse=True)

    rank = 0
    for v, peg, peg_score in flat:
        ng = select_pe3_ngrna(
            v.ngRNA_candidates,
            min_distance=ngrna_distance_min,
            max_distance=ngrna_distance_max,
            pegrna_strand=peg.target_strand,
        )
        if ng is None:
            continue
        rank += 1
        pairs.append(
            PE7DesignPair(
                rank=rank,
                variant_sequence_name=v.sequence_name,
                pegRNA_full=peg.pegrna_full,
                spacer=peg.spacer,
                pbs=peg.pbs,
                rt=peg.rt,
                pegrna_strand=peg.target_strand,
                pegrna_pe2_score=peg_score,
                pegrna_pe2_cell_line=score_cell_line,
                ngrna_spacer=ng.spacer,
                ngrna_distance_to_edit=ng.distance_to_edit,
                ngrna_strand=ng.target_strand,
                ngrna_deepcas9_score=float(ng.deepcas9_score),
                ngrna_nick_in_mrna=edit_position + ng.distance_to_edit,
                joint_pe3_score=peg_score * float(ng.deepcas9_score) / 100.0,
            )
        )
        if rank >= top_n:
            break
    return pairs


def run_pe7_workflow(
    *,
    gene: GeneRecord,
    insert_full_dna: str,
    variant_scores: Sequence[VariantScore],
    ngrna_distance_min: int = 40,
    ngrna_distance_max: int = 100,
    score_cell_line: str = DEFAULT_SCORE_CELL_LINE,
    top_n: int = DEFAULT_TOP_N,
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
    score_cell_line:
        ``"k562"`` (default) or ``"hek293t"`` — which PRIDICT 2.0 PE2
        column drives ranking.
    top_n:
        How many ranked (pegRNA, ngRNA) pairs to return (default 20).

    Returns
    -------
    PE7Result
        ``top_pairs`` carries up to ``top_n`` ranked pairs; legacy top-
        level fields mirror ``top_pairs[0]``. Both are empty if PRIDICT
        returned no usable candidates.

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

    # Propagate score_cell_line so VariantScore.best_pegRNA / top_score use
    # the requested column.
    variants = [_with_score_cell_line(v, score_cell_line) for v in variant_scores]
    ranked = rank_variants_by_pe2(variants)
    winning: Optional[VariantScore] = ranked[0] if ranked else None

    edit_position_in_mrna = gene.stop_codon_start
    pairs = _build_top_pairs(
        variants=variants,
        score_cell_line=score_cell_line,
        ngrna_distance_min=ngrna_distance_min,
        ngrna_distance_max=ngrna_distance_max,
        edit_position=edit_position_in_mrna,
        top_n=top_n,
    )

    if not pairs:
        return _empty_result(
            gene=gene,
            insert_full_dna=insert_full_dna,
            variants=variants,
            score_cell_line=score_cell_line,
        )

    rank1 = pairs[0]
    return PE7Result(
        winning_variant=winning,
        full_insert_dna=insert_full_dna,
        pegRNA_full=rank1.pegRNA_full,
        spacer=rank1.spacer,
        pbs=rank1.pbs,
        rt=rank1.rt,
        pegrna_strand=rank1.pegrna_strand,
        pegrna_pe2_score=rank1.pegrna_pe2_score,
        pegrna_pe2_cell_line=rank1.pegrna_pe2_cell_line,
        ngrna_spacer=rank1.ngrna_spacer,
        ngrna_distance_to_edit=rank1.ngrna_distance_to_edit,
        ngrna_strand=rank1.ngrna_strand,
        ngrna_deepcas9_score=rank1.ngrna_deepcas9_score,
        edit_position_in_mrna=edit_position_in_mrna,
        ngrna_nick_in_mrna=rank1.ngrna_nick_in_mrna,
        joint_pe3_score=rank1.joint_pe3_score,
        top_pairs=pairs,
        all_evaluated_variants=variants,
    )


def _with_score_cell_line(variant: VariantScore, score_cell_line: str) -> VariantScore:
    """Return a copy of ``variant`` with ``score_cell_line`` set."""
    return VariantScore(
        sequence_name=variant.sequence_name,
        codon_variant_dna=variant.codon_variant_dna,
        pegRNA_designs=variant.pegRNA_designs,
        ngRNA_candidates=variant.ngRNA_candidates,
        score_cell_line=score_cell_line,
    )
