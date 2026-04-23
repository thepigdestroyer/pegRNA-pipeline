"""Edge case: no compatible pegRNA pairs survive the filter.

Plan 2 v1.1 requires :func:`run_pa_workflow` to still return a
:class:`PAResult` (not raise) — with ``winning_pair=None`` and counts
that make the failure explicit.
"""

from __future__ import annotations

from pegrna_pipeline.hpc.data_types import (
    GeneRecord,
    PegRNADesign,
    VariantScore,
)
from pegrna_pipeline.hpc.pa_workflow import PAResult, run_pa_workflow


def _gene() -> GeneRecord:
    upstream = "ACGT" * 250
    downstream = "TGCA" * 250
    return GeneRecord(
        gene_name="X",
        refseq_accession="NM_X",
        mrna_sequence=upstream + "TAA" + downstream,
        cds_start=0,
        cds_end=1003,
        stop_codon_start=1000,
        stop_codon="TAA",
        upstream_flank=upstream[-100:],
        downstream_flank=downstream[:100],
        flanking_bp=100,
    )


def _peg(name: str, strand: str, position: int) -> PegRNADesign:
    return PegRNADesign(
        sequence_name=name,
        spacer="A" * 20,
        pbs="C" * 13,
        rt="G" * 30,
        pbs_length=13,
        rt_length=30,
        target_strand=strand,
        editing_position=position,
        pe2_hek293t=0.5,
        pe2_k562=0.4,
        pegrna_full="A" * 140,
    )


def test_same_strand_only_yields_empty_result():
    # All Fw — no opposite-strand pairings possible.
    up = VariantScore(
        sequence_name="up",
        codon_variant_dna="GCA",
        pegRNA_designs=[_peg("u1", "Fw", 50)],
    )
    dn = VariantScore(
        sequence_name="dn",
        codon_variant_dna="GCC",
        pegRNA_designs=[_peg("d1", "Fw", 100)],
    )
    result = run_pa_workflow(
        gene=_gene(),
        upstream_variant_scores=[up],
        downstream_variant_scores=[dn],
        donor_dna="ACGT" * 10,
    )
    assert isinstance(result, PAResult)
    assert result.winning_pair is None
    assert result.top_pairs == []
    assert result.pairs_after_compat_filter == 0
    assert result.total_pairs_evaluated == 1


def test_all_pairs_fail_distance_filter():
    up = VariantScore(
        sequence_name="up",
        codon_variant_dna="GCA",
        pegRNA_designs=[_peg("u1", "Fw", 50)],
    )
    dn = VariantScore(
        sequence_name="dn",
        codon_variant_dna="GCC",
        pegRNA_designs=[_peg("d1", "Rv", 500)],  # too far
    )
    result = run_pa_workflow(
        gene=_gene(),
        upstream_variant_scores=[up],
        downstream_variant_scores=[dn],
        donor_dna="ACGT" * 10,
    )
    assert result.winning_pair is None


def test_empty_pool_returns_empty_pair_list():
    result = run_pa_workflow(
        gene=_gene(),
        upstream_variant_scores=[],
        downstream_variant_scores=[],
        donor_dna="ACGT" * 10,
    )
    assert result.winning_pair is None
    assert result.total_pairs_evaluated == 0
    assert result.pairs_after_compat_filter == 0
