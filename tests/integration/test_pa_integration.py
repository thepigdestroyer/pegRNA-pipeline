"""PA integration test verifying donor 5' == upstream_flap, 3' == downstream_flap."""

from __future__ import annotations

from pathlib import Path

from pegrna_pipeline.hpc.data_types import GeneRecord
from pegrna_pipeline.hpc.pa_workflow import PAResult, run_pa_workflow

FIXTURES = Path(__file__).parent.parent / "fixtures" / "pridict"
UPSTREAM_DIR = FIXTURES / "pa_upstream"
DOWNSTREAM_DIR = FIXTURES / "pa_downstream"


UPSTREAM_FLAP = "GCAGCAGCAGCAGCAGCAGCAGCAGCAGCA"
DOWNSTREAM_FLAP = "GCCGCCGCCGCCGCCGCCGCCGCCGCCGCC"
BULK = "AAACCCGGGTTTAAACCCGGGTTT"


def _gene() -> GeneRecord:
    upstream = "ACGT" * 250
    downstream = "TGCA" * 250
    return GeneRecord(
        gene_name="CTCF",
        refseq_accession="NM_006565.4",
        mrna_sequence=upstream + "TAA" + downstream,
        cds_start=0,
        cds_end=1003,
        stop_codon_start=1000,
        stop_codon="TAA",
        upstream_flank=upstream[-100:],
        downstream_flank=downstream[:100],
        flanking_bp=100,
    )


def test_donor_has_flaps_at_its_ends():
    result = run_pa_workflow(
        gene=_gene(),
        upstream_fixtures_dir=UPSTREAM_DIR,
        downstream_fixtures_dir=DOWNSTREAM_DIR,
        upstream_flap_dna=UPSTREAM_FLAP,
        downstream_flap_dna=DOWNSTREAM_FLAP,
        bulk_dna=BULK,
    )
    assert isinstance(result, PAResult)
    # 5' end == upstream flap.
    assert result.linear_donor_dna[: len(UPSTREAM_FLAP)] == UPSTREAM_FLAP
    # 3' end == downstream flap.
    assert result.linear_donor_dna[-len(DOWNSTREAM_FLAP) :] == DOWNSTREAM_FLAP
    # Middle is bulk.
    assert result.linear_donor_dna[len(UPSTREAM_FLAP) : -len(DOWNSTREAM_FLAP)] == BULK


def test_integrated_product_contains_donor_and_mrna_flanks():
    gene = _gene()
    result = run_pa_workflow(
        gene=gene,
        upstream_fixtures_dir=UPSTREAM_DIR,
        downstream_fixtures_dir=DOWNSTREAM_DIR,
        upstream_flap_dna=UPSTREAM_FLAP,
        downstream_flap_dna=DOWNSTREAM_FLAP,
        bulk_dna=BULK,
    )
    assert result.full_integrated_product.startswith(gene.upstream_flank)
    assert result.full_integrated_product.endswith(gene.downstream_flank)
    assert result.linear_donor_dna in result.full_integrated_product


def test_winning_pair_opposite_strands():
    result = run_pa_workflow(
        gene=_gene(),
        upstream_fixtures_dir=UPSTREAM_DIR,
        downstream_fixtures_dir=DOWNSTREAM_DIR,
        upstream_flap_dna=UPSTREAM_FLAP,
        downstream_flap_dna=DOWNSTREAM_FLAP,
        bulk_dna=BULK,
    )
    pair = result.winning_pair
    assert pair is not None
    assert pair.upstream.target_strand != pair.downstream.target_strand


def test_top_pairs_sorted_by_joint_score():
    result = run_pa_workflow(
        gene=_gene(),
        upstream_fixtures_dir=UPSTREAM_DIR,
        downstream_fixtures_dir=DOWNSTREAM_DIR,
        upstream_flap_dna=UPSTREAM_FLAP,
        downstream_flap_dna=DOWNSTREAM_FLAP,
        bulk_dna=BULK,
    )
    scores = [p.joint_score for p in result.top_pairs]
    assert scores == sorted(scores, reverse=True)
