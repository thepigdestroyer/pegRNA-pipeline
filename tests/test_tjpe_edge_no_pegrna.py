"""T11: the selected variant has no pegRNA designs → degenerate TJPEResult.

Run the workflow with a single variant that PRIDICT couldn't score. We
expect a :class:`TJPEResult` whose winning_variant points to the empty
variant, with zero-length pegRNA / PBS fields, and
``pegrna_pe2_hek293t == 0.0``.
"""

from __future__ import annotations

from pegrna_pipeline.hpc.data_types import GeneRecord, VariantScore
from pegrna_pipeline.hpc.tjpe_workflow import run_tjpe_workflow


def _gene() -> GeneRecord:
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


def test_variant_without_pegrna_designs_returns_degenerate_result():
    empty_variant = VariantScore(
        sequence_name="empty_v1",
        codon_variant_dna="ACG" * 10,
        pegRNA_designs=[],
        ngRNA_candidates=[],
    )
    result = run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna="ATGCATGC",
        variant_scores=[empty_variant],
    )
    # Variant is still reported — caller needs to see it to diagnose.
    assert result.winning_variant is empty_variant
    assert result.all_evaluated_variants == [empty_variant]
    # pegRNA fields are empty since we had nothing to assemble.
    assert result.tjpe_pegRNA_sequence == ""
    assert result.spacer == ""
    assert result.pbs1 == ""
    assert result.rc_pbs2 == ""
    assert result.rc_pbs2_length == 0
    assert result.pegrna_pe2_hek293t == 0.0


def test_empty_variant_list_returns_none_winner():
    result = run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna="ATGCATGC",
        variant_scores=[],
    )
    assert result.winning_variant is None
    assert result.tjpe_pegRNA_sequence == ""
    assert result.all_evaluated_variants == []
