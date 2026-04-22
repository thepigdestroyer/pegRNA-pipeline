"""Tests for :func:`compute_flap_boundaries`.

Splits an amino-acid string into three parts::

    upstream_flap (flap_length_aa)  +  bulk  +  downstream_flap (flap_length_aa)

The flaps are the regions covered by the two pegRNA RT extensions in a
Prime-Assembly (PA) construct. Default flap length is 10 aa per Plan 2
v1.1.
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.pa_workflow import compute_flap_boundaries


def test_default_flap_length_is_10_aa():
    full = "A" * 10 + "B" * 5 + "C" * 10  # 25 aa total
    up, bulk, down = compute_flap_boundaries(full)
    assert up == "A" * 10
    assert bulk == "B" * 5
    assert down == "C" * 10


def test_custom_flap_length():
    full = "A" * 5 + "B" * 10 + "C" * 5
    up, bulk, down = compute_flap_boundaries(full, flap_length_aa=5)
    assert up == "A" * 5
    assert bulk == "B" * 10
    assert down == "C" * 5


def test_minimum_length_exactly_2x_flap_gives_empty_bulk():
    full = "A" * 10 + "C" * 10  # 20 aa, exactly 2x default flap
    up, bulk, down = compute_flap_boundaries(full)
    assert up == "A" * 10
    assert bulk == ""
    assert down == "C" * 10


def test_too_short_raises():
    with pytest.raises(ValueError):
        compute_flap_boundaries("A" * 19)  # < 2 * 10


def test_too_short_with_custom_flap_raises():
    with pytest.raises(ValueError):
        compute_flap_boundaries("A" * 9, flap_length_aa=5)  # need >= 10


def test_tuple_concatenates_to_original():
    full = "MASLKDKEHVFGELRGFPKRLYCKNGGFFLRIHPDGRVDGVRE"
    up, bulk, down = compute_flap_boundaries(full, flap_length_aa=10)
    assert up + bulk + down == full


def test_empty_string_raises():
    with pytest.raises(ValueError):
        compute_flap_boundaries("")
