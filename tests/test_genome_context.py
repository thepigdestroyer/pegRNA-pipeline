"""Task 7: diploid GenomicContext + coordinate conversions + fail-loud guard."""

import pytest

from pegrna_pipeline.hpc.genome_context import (
    ExonBlock,
    GenomicContext,
    Haplotype,
    MissingGenomicContextError,
    mrna_to_genomic,
)


def _ctx():
    # two exons: mRNA[0:10) <- genome[100:110); mRNA[10:20) <- genome[200:210)
    exons = (ExonBlock(100, 110, 0, 10), ExonBlock(200, 210, 10, 20))
    hap1 = Haplotype(name="hap1", genomic_sequence=None, exons=exons)
    hap2 = Haplotype(name="hap2", genomic_sequence=None, exons=exons)
    return GenomicContext(
        genome_build="T2T-KOLF2.1Jv1.1",
        chrom="chr1",
        strand=1,
        haplotypes=(hap1, hap2),
        het_snps=({"pos": 105, "ref": "A", "alt": "G", "phased": True},),
    )


def test_mrna_to_genomic_maps_across_exons():
    ctx = _ctx()
    assert mrna_to_genomic(ctx, 5) == 105
    assert mrna_to_genomic(ctx, 12) == 202


def test_mrna_to_genomic_out_of_range_raises():
    with pytest.raises(ValueError):
        mrna_to_genomic(_ctx(), 25)


def test_missing_context_guard():
    with pytest.raises(MissingGenomicContextError):
        mrna_to_genomic(None, 5)


def test_diploid_carries_two_haplotypes():
    ctx = _ctx()
    assert len(ctx.haplotypes) == 2
    assert ctx.haplotypes[0].name == "hap1"
    assert ctx.het_snps[0]["phased"] is True


def test_to_grch38_identity_when_no_chain():
    # With no chain loaded, to_grch38 is the identity map but flags it.
    ctx = _ctx()
    assert ctx.to_grch38(105) == 105  # placeholder identity until a real chain is wired


def test_per_haplotype_mapping_uses_specific_exons():
    ctx = _ctx()
    assert mrna_to_genomic(ctx, 5, haplotype=0) == 105
    assert mrna_to_genomic(ctx, 5, haplotype=1) == 105


def test_minus_strand_maps_from_g_end_downward():
    # strand=-1: mRNA 5' end -> highest genomic coord; increasing m_pos -> lower.
    exons = (ExonBlock(200, 209, 0, 10),)  # 10-nt exon, g 200..209 inclusive
    hap = Haplotype(name="h1", genomic_sequence=None, exons=exons)
    hap2 = Haplotype(name="h2", genomic_sequence=None, exons=exons)
    ctx = GenomicContext(
        genome_build="GRCh38",
        chrom="chr1",
        strand=-1,
        haplotypes=(hap, hap2),
        het_snps=(),
    )
    assert mrna_to_genomic(ctx, 0) == 209  # first mRNA base -> g_end
    assert mrna_to_genomic(ctx, 9) == 200  # last mRNA base -> g_start
    # plus-strand sanity: same exon as plus reads upward from g_start
    ctxp = GenomicContext(
        genome_build="GRCh38",
        chrom="chr1",
        strand=1,
        haplotypes=(hap, hap2),
        het_snps=(),
    )
    assert mrna_to_genomic(ctxp, 0) == 200
    assert mrna_to_genomic(ctxp, 9) == 209
