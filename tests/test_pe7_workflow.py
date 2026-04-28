"""Unit tests for :func:`run_pe7_workflow`.

PE7 mirrors TJ-PE for the parts that overlap (rank-by-pe2, opposite-strand
ngRNA selection) and adds:

* a hard ≤40 bp insert-size guard,
* the ``insert_full_dna`` validation (non-empty, ACGT not enforced here —
  caller responsibility, same as TJ-PE / PA),
* no second-nick PBS2 derivation (output uses PRIDICT's ``pegRNA_full``
  string verbatim).
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.data_types import (
    GeneRecord,
    NgRNADesign,
    PegRNADesign,
    VariantScore,
)
from pegrna_pipeline.hpc.pe7_workflow import (
    PE7Result,
    PRIDICT_MAX_INSERT_BP,
    run_pe7_workflow,
)


def _gene() -> GeneRecord:
    """Synthetic gene with the stop codon at position 1000."""
    upstream = "ACGT" * 250
    downstream = "TGCA" * 250
    mrna = upstream + "TAA" + downstream
    return GeneRecord(
        gene_name="FAKE",
        refseq_accession="NM_FAKE.1",
        mrna_sequence=mrna,
        cds_start=0,
        cds_end=1003,
        stop_codon_start=1000,
        stop_codon="TAA",
        upstream_flank=upstream[-100:],
        downstream_flank=downstream[:100],
        flanking_bp=100,
    )


def _peg(
    *,
    name: str = "pe7_var_0",
    spacer: str = "ACGTACGTACGTACGTACGT",
    score: float = 0.7,
    strand: str = "Fw",
    pegrna_full: str = "PEGFULL",
) -> PegRNADesign:
    # K562 is the default ranking column, so we make ``score`` the K562
    # value directly — assertions then read straight off ``score``.
    return PegRNADesign(
        sequence_name=name,
        spacer=spacer,
        pbs="ACGTACGTACGTACGTA",
        rt="GCATGCATGCATGCATGCATGCATGCATGCAT",
        pbs_length=17,
        rt_length=32,
        target_strand=strand,  # type: ignore[arg-type]
        editing_position=50,
        pe2_hek293t=score + 0.1,  # arbitrary; not used under K562 ranking
        pe2_k562=score,
        pegrna_full=pegrna_full,
    )


def _ng(
    *,
    spacer: str = "TTACGTACGTACGTACGTAC",
    distance: int = 50,
    strand: str = "Rv",
    score: float = 0.85,
) -> NgRNADesign:
    return NgRNADesign(
        spacer=spacer,
        distance_to_edit=distance,
        target_strand=strand,  # type: ignore[arg-type]
        deepcas9_score=score,
        pe3b=False,
        pam_disrupt=False,
    )


def _variant(
    pegs: list[PegRNADesign] | None = None,
    ngs: list[NgRNADesign] | None = None,
) -> VariantScore:
    return VariantScore(
        sequence_name="pe7_var_0",
        codon_variant_dna="",
        pegRNA_designs=pegs if pegs is not None else [_peg()],
        ngRNA_candidates=ngs if ngs is not None else [_ng()],
    )


# -------------------------------------------------------- happy path


def test_run_pe7_picks_top_pegrna_and_ngrna():
    insert = "ATG" * 5  # 15 bp
    pegs = [
        _peg(name="lo", score=0.10, spacer="GAATTCGTGTTCCTAAGCGT"),
        _peg(name="hi", score=0.80, spacer="ACGTACGTACGTACGTACGT"),
    ]
    ngs = [
        _ng(distance=50, strand="Rv", score=0.40),
        _ng(distance=70, strand="Rv", score=0.95, spacer="GACGTACGTACGTACGTACG"),
    ]
    result = run_pe7_workflow(
        gene=_gene(),
        insert_full_dna=insert,
        variant_scores=[_variant(pegs=pegs, ngs=ngs)],
    )
    assert isinstance(result, PE7Result)
    assert result.full_insert_dna == insert
    assert result.pegrna_pe2_score == pytest.approx(0.80)
    assert result.spacer == "ACGTACGTACGTACGTACGT"
    assert result.ngrna_deepcas9_score == pytest.approx(0.95)
    # ngRNA must be opposite strand of pegRNA
    assert result.pegrna_strand == "Fw"
    assert result.ngrna_strand == "Rv"
    # nick arithmetic mirrors TJ-PE
    assert result.edit_position_in_mrna == 1000
    assert result.ngrna_nick_in_mrna == 1000 + result.ngrna_distance_to_edit
    # PRIDICT's pegRNA_full string passes through unmodified
    assert result.pegRNA_full == "PEGFULL"
    # Per-variant tracking is preserved
    assert len(result.all_evaluated_variants) == 1


def test_run_pe7_picks_opposite_strand_when_pegrna_is_rv():
    insert = "ATG" * 5
    pegs = [_peg(strand="Rv")]
    ngs = [
        _ng(strand="Fw", distance=50, score=0.6),
        _ng(strand="Rv", distance=50, score=0.9),  # same strand → must be filtered out
    ]
    result = run_pe7_workflow(
        gene=_gene(),
        insert_full_dna=insert,
        variant_scores=[_variant(pegs=pegs, ngs=ngs)],
    )
    assert result.ngrna_strand == "Fw"
    assert result.ngrna_deepcas9_score == pytest.approx(0.6)


# -------------------------------------------------------- size guard


def test_run_pe7_rejects_oversized_insert():
    too_big = "A" * (PRIDICT_MAX_INSERT_BP + 1)
    with pytest.raises(ValueError, match=r"PE7 mode supports inserts"):
        run_pe7_workflow(
            gene=_gene(),
            insert_full_dna=too_big,
            variant_scores=[_variant()],
        )


def test_run_pe7_accepts_max_size_insert():
    insert = "A" * PRIDICT_MAX_INSERT_BP
    result = run_pe7_workflow(
        gene=_gene(),
        insert_full_dna=insert,
        variant_scores=[_variant()],
    )
    assert result.full_insert_dna == insert


def test_run_pe7_rejects_empty_insert():
    with pytest.raises(ValueError, match=r"non-empty"):
        run_pe7_workflow(
            gene=_gene(),
            insert_full_dna="",
            variant_scores=[_variant()],
        )


# -------------------------------------------------------- degenerate paths


def test_run_pe7_with_no_variants_returns_empty_fields():
    result = run_pe7_workflow(
        gene=_gene(),
        insert_full_dna="ATG",
        variant_scores=[],
    )
    assert result.spacer == ""
    assert result.ngrna_spacer == ""
    assert result.pegrna_pe2_score == 0.0
    assert result.edit_position_in_mrna == 1000


def test_run_pe7_with_no_usable_ngrna_returns_empty_pairs():
    """All ngRNAs fail the PE3 distance window — no usable PE7 pair exists.

    Under top-N semantics, a PE7 design needs both a pegRNA AND a valid
    PE3 ngRNA. So this case yields an empty result (top_pairs == []),
    not a pegRNA-only one.
    """
    pegs = [_peg(score=0.5)]
    ngs = [_ng(distance=10, score=0.9)]  # below the 40 bp PE3 floor
    result = run_pe7_workflow(
        gene=_gene(),
        insert_full_dna="ATG",
        variant_scores=[_variant(pegs=pegs, ngs=ngs)],
    )
    assert result.top_pairs == []
    assert result.spacer == ""
    assert result.ngrna_spacer == ""
    assert result.pegrna_pe2_score == 0.0
    assert result.ngrna_deepcas9_score == 0.0
