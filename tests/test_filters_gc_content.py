"""Tests for _gc_content returning a fraction between 0.0 and 1.0."""

import pytest

from pegrna_pipeline.hpc.codon_expander import _gc_content


def test_pure_at_is_zero():
    assert _gc_content("AAAT") == 0.0


def test_pure_gc_is_one():
    assert _gc_content("GGGC") == 1.0


def test_balanced_is_half():
    assert _gc_content("ATGC") == pytest.approx(0.5)


def test_mixed_case_is_handled():
    assert _gc_content("aAgG") == pytest.approx(0.5)


def test_empty_sequence_is_zero():
    assert _gc_content("") == 0.0
