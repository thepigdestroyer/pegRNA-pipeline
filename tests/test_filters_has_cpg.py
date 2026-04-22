"""Tests for _has_cpg.

The CpG filter is the literal substring 'CG'. Note that 'GC' does NOT count.
Per the v1.2 corrections: "GACGAT" contains 'CG' at position 2 -> True.
"""

from pegrna_pipeline.hpc.codon_expander import _has_cpg


def test_simple_cpg_present():
    assert _has_cpg("GACGAT") is True


def test_no_cpg_a_rich():
    assert _has_cpg("AAAAAA") is False


def test_no_cpg_gc_only():
    # 'GC' alone is not a CpG.
    assert _has_cpg("AGCAGC") is False


def test_no_cpg_in_gacaat():
    assert _has_cpg("GACAAT") is False


def test_cpg_at_start():
    assert _has_cpg("CGTAAA") is True


def test_cpg_at_end():
    assert _has_cpg("AAAACG") is True


def test_empty_sequence():
    assert _has_cpg("") is False
