"""Tests for CodonVariant dataclass."""

from dataclasses import fields

import pytest

from pegrna_pipeline.hpc.codon_expander import CodonVariant


def test_codon_variant_has_expected_fields():
    names = {f.name for f in fields(CodonVariant)}
    assert names == {"dna", "c_count", "gc_content", "rank"}


def test_codon_variant_is_immutable_frozen_dataclass():
    v = CodonVariant(dna="ATG", c_count=0, gc_content=0.33, rank=0)
    with pytest.raises(Exception):
        v.dna = "CCC"


def test_codon_variant_constructs_correctly():
    v = CodonVariant(dna="GCC", c_count=2, gc_content=1.0, rank=0)
    assert v.dna == "GCC"
    assert v.c_count == 2
    assert v.gc_content == 1.0
    assert v.rank == 0


def test_codon_variants_are_equatable():
    a = CodonVariant(dna="ATG", c_count=0, gc_content=0.33, rank=1)
    b = CodonVariant(dna="ATG", c_count=0, gc_content=0.33, rank=1)
    assert a == b


def test_codon_variants_are_hashable():
    v = CodonVariant(dna="ATG", c_count=0, gc_content=0.33, rank=0)
    {v}  # should not raise
