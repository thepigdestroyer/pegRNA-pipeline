"""Coordinate-safe diploid materialization (scope A). RED-first.

`apply_phased_variants` builds a per-haplotype sequence by applying KOLF2.1J phased
SUBSTITUTIONS onto a GRCh38 region sequence, staying in GRCh38 coordinates so the
GRCh38-trained gates remain valid (`to_grch38` identity holds). Substitutions only:
indels are skipped for the coordinate-anchored sequence (indel-as-truth = scope B).
Fail-loud on a ref-base mismatch (wrong build/coordinates).
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.genome_context import GenomicContext
from pegrna_pipeline.hpc.genome_ingest import build_genomic_context
from pegrna_pipeline.hpc.variant_apply import apply_phased_variants

REF = "ACGTACGT"  # genomic region, forward strand
START = 100  # 1-based genomic coordinate of REF[0]


def test_het_snv_applies_to_one_haplotype_only():
    # pos 102 (index 2) G->T, present on haplotype 1 only.
    variants = [{"pos": 102, "ref": "G", "alt": "T", "haplotype": 1, "phased": True}]
    hap0 = apply_phased_variants(REF, variants, region_start=START, haplotype=0)
    hap1 = apply_phased_variants(REF, variants, region_start=START, haplotype=1)
    assert hap0 == REF  # untouched haplotype == reference
    assert hap1 == "ACTTACGT"  # index 2 flipped G->T


def test_homozygous_snv_applies_to_both_haplotypes():
    variants = [{"pos": 100, "ref": "A", "alt": "C", "homozygous": True}]
    hap0 = apply_phased_variants(REF, variants, region_start=START, haplotype=0)
    hap1 = apply_phased_variants(REF, variants, region_start=START, haplotype=1)
    assert hap0 == "CCGTACGT"
    assert hap1 == "CCGTACGT"


def test_ref_base_mismatch_is_fail_loud():
    # Claims ref 'A' at pos 102 where REF actually has 'G' -> wrong build/coords.
    bad = [{"pos": 102, "ref": "A", "alt": "T", "haplotype": 0, "phased": True}]
    with pytest.raises(ValueError):
        apply_phased_variants(REF, bad, region_start=START, haplotype=0)


def test_variant_outside_window_ignored():
    variants = [{"pos": 999, "ref": "A", "alt": "T", "haplotype": 0, "phased": True}]
    assert apply_phased_variants(REF, variants, region_start=START, haplotype=0) == REF


def test_indels_skipped_in_coordinate_safe_mode():
    ins = [{"pos": 101, "ref": "C", "alt": "CGG", "haplotype": 0, "phased": True}]
    dele = [{"pos": 101, "ref": "CGT", "alt": "C", "haplotype": 0, "phased": True}]
    assert apply_phased_variants(REF, ins, region_start=START, haplotype=0) == REF
    assert apply_phased_variants(REF, dele, region_start=START, haplotype=0) == REF


def test_empty_variants_returns_reference():
    assert apply_phased_variants(REF, [], region_start=START, haplotype=0) == REF
    assert apply_phased_variants(REF, [], region_start=START, haplotype=1) == REF


def test_multiple_snvs_on_one_haplotype_all_applied():
    variants = [
        {"pos": 100, "ref": "A", "alt": "T", "haplotype": 0, "phased": True},
        {"pos": 104, "ref": "A", "alt": "C", "haplotype": 0, "phased": True},
    ]
    hap0 = apply_phased_variants(REF, variants, region_start=START, haplotype=0)
    assert hap0 == "TCGTCCGT"  # index0 A->T, index4 A->C


def test_build_genomic_context_materializes_real_diploid():
    lookup = {
        "Exon": [{"start": 100, "end": 107}],
        "seq_region_name": "1",
        "strand": 1,
    }
    variants = [{"pos": 103, "ref": "T", "alt": "A", "haplotype": 1, "phased": True}]
    ctx = build_genomic_context(lookup=lookup, mrna_seq=REF, het_snps=(), phased_variants=variants)
    assert isinstance(ctx, GenomicContext)
    hap1_seq = ctx.haplotypes[0].genomic_sequence
    hap2_seq = ctx.haplotypes[1].genomic_sequence
    assert hap1_seq == REF  # reference haplotype
    assert hap2_seq == "ACGAACGT"  # index 3 T->A on hap2
    assert hap1_seq != hap2_seq  # a REAL diploid, not two copies of reference
    # Coordinate-safe: GRCh38 frame preserved, identity to_grch38 still valid.
    assert "GRCh38" in ctx.genome_build
    assert ctx.to_grch38(103) == 103


def test_build_genomic_context_backward_compatible_without_variants():
    lookup = {"Exon": [{"start": 100, "end": 107}], "seq_region_name": "1", "strand": 1}
    ctx = build_genomic_context(lookup=lookup, mrna_seq=REF)
    assert ctx.haplotypes[0].genomic_sequence == REF
    assert ctx.haplotypes[1].genomic_sequence == REF  # unchanged legacy behaviour
