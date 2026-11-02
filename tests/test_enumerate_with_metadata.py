"""Tests for enumerate_with_metadata returning EnumerationResult.

EnumerationResult has fields:
  - variants: list[CodonVariant]
  - theoretical_max: int -- the number of synonymous combinations without filters
  - was_capped: bool -- True if more valid variants existed than max_variants
"""

import math
from dataclasses import fields

from pegrna_pipeline.hpc.codon_expander import (
    CODONS_BY_AA,
    EnumerationResult,
    enumerate_with_metadata,
)


def test_enumeration_result_has_expected_fields():
    names = {f.name for f in fields(EnumerationResult)}
    assert names == {"variants", "theoretical_max", "was_capped"}


def test_theoretical_max_is_product_of_codon_counts():
    aa = "MKL"
    expected = math.prod(len(CODONS_BY_AA[a]) for a in aa)
    result = enumerate_with_metadata(aa, max_variants=10)
    assert result.theoretical_max == expected


def test_theoretical_max_for_methionine_only_is_one():
    result = enumerate_with_metadata("M", max_variants=10)
    assert result.theoretical_max == 1


def test_was_capped_true_when_more_variants_exist():
    # L has 6 codons, so "LL" has 36 combinations; easy to exceed cap=5.
    result = enumerate_with_metadata("LL", max_variants=5)
    assert result.was_capped is True
    assert len(result.variants) == 5


def test_was_capped_false_when_cap_is_high():
    result = enumerate_with_metadata("M", max_variants=10)
    assert result.was_capped is False


def test_empty_aa_returns_zero_max_no_variants():
    result = enumerate_with_metadata("", max_variants=10)
    assert result.theoretical_max == 0
    assert result.variants == []
    assert result.was_capped is False
