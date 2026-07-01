"""Real genomic-context ingest (Ensembl GRCh38 + optional KOLF VCF overlay).

The HTTP layer is injected so CI tests run offline against canned Ensembl-shaped
responses; the real path uses ``requests`` on a networked host.
"""

import pytest

from pegrna_pipeline.hpc.genome_context import GenomicContext, MissingGenomicContextError
from pegrna_pipeline.hpc.genome_ingest import build_genomic_context, fetch_genomic_context

# Canned Ensembl /lookup (expand=1) + /sequence responses for a 2-exon transcript.
_LOOKUP = {
    "seq_region_name": "1",
    "strand": 1,
    "Exon": [
        {"start": 1001, "end": 1010, "strand": 1},
        {"start": 2001, "end": 2010, "strand": 1},
    ],
}
_SEQ = {"seq": "ACGTACGTAC" "TTTTGGGGCC"}  # exon1 (10) + exon2 (10), spliced mRNA


def _fake_http(url, **kw):
    if "/lookup/" in url:
        return _LOOKUP
    if "/sequence/" in url:
        return _SEQ
    raise AssertionError(f"unexpected url {url}")


def test_build_genomic_context_from_ensembl_payload():
    ctx = build_genomic_context(lookup=_LOOKUP, mrna_seq=_SEQ["seq"], het_snps=())
    assert isinstance(ctx, GenomicContext)
    assert ctx.genome_build.startswith("GRCh38")
    assert ctx.chrom == "1"
    assert ctx.strand == 1
    # two exon blocks, mRNA 0..10 -> genome 1001.., 10..20 -> genome 2001..
    assert len(ctx.haplotypes[0].exons) == 2
    assert ctx.haplotypes[0].exons[0].g_start == 1001
    # diploid: two haplotypes present even without variants (identical to ref)
    assert len(ctx.haplotypes) == 2


def test_fetch_genomic_context_uses_injected_http():
    ctx = fetch_genomic_context("ENST00000TEST", http_get=_fake_http)
    assert ctx.chrom == "1"
    assert ctx.haplotypes[0].exons[1].g_start == 2001


def test_vcf_overlay_creates_distinct_alt_haplotype():
    # a het SNP at genomic 1005 (A->G) makes hap2 differ from hap1 at that base
    het = [{"chrom": "1", "pos": 1005, "ref": "A", "alt": "G", "phased": True}]
    ctx = build_genomic_context(lookup=_LOOKUP, mrna_seq=_SEQ["seq"], het_snps=het)
    assert ctx.het_snps and ctx.het_snps[0]["pos"] == 1005
    assert ctx.het_snps[0]["phased"] is True


def test_missing_transcript_raises():
    def _empty_http(url, **kw):
        return {}

    with pytest.raises((MissingGenomicContextError, KeyError, ValueError)):
        fetch_genomic_context("ENST00000MISSING", http_get=_empty_http)
