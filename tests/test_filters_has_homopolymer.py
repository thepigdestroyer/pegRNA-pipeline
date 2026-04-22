"""Tests for _has_homopolymer with min_run parameter (regex (.)\\1{min_run-1,})."""

from pegrna_pipeline.hpc.codon_expander import _has_homopolymer


def test_no_homopolymer_below_threshold():
    assert _has_homopolymer("AAACCC", min_run=5) is False


def test_homopolymer_exactly_at_threshold():
    assert _has_homopolymer("AAAAA", min_run=5) is True


def test_homopolymer_above_threshold():
    assert _has_homopolymer("TTTTTTTT", min_run=5) is True


def test_mixed_sequence_with_embedded_run():
    assert _has_homopolymer("ACGCCCCCGAT", min_run=5) is True


def test_different_bases_same_char_counts():
    assert _has_homopolymer("GGGGG", min_run=5) is True
    assert _has_homopolymer("CCCCC", min_run=5) is True


def test_various_min_runs():
    assert _has_homopolymer("AAAA", min_run=4) is True
    assert _has_homopolymer("AAA", min_run=4) is False
    assert _has_homopolymer("AAA", min_run=3) is True


def test_empty():
    assert _has_homopolymer("", min_run=5) is False
