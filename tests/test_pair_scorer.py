"""Tests for :mod:`pegrna_pipeline.hpc.pair_scorer`.

A :class:`PegRNAPair` is a compatible upstream/downstream pegRNA pairing
for Prime Assembly. Compatibility rules (Plan 2 v1.1):

* Opposite ``Target-Strand`` values (one Fw, one Rv).
* ``|editing_position_up - editing_position_dn|`` within
  [min_distance, max_distance] — default [30, 80].
* Distance must be >= ``MIN_NICK_SEPARATION_BP`` (=23).

Joint score uses the geometric mean of the two pe2_hek293t values
(``sqrt(up_score * dn_score)``).
"""

from __future__ import annotations

import math

import pytest

from pegrna_pipeline.hpc.data_types import PegRNADesign
from pegrna_pipeline.hpc.pair_scorer import (
    MIN_NICK_SEPARATION_BP,
    PegRNAPair,
    filter_compatible_pairs,
)


def _peg(
    name: str,
    strand: str,
    position: int,
    score: float,
) -> PegRNADesign:
    return PegRNADesign(
        sequence_name=name,
        spacer="ACGT" * 5,
        pbs="A" * 13,
        rt="C" * 30,
        pbs_length=13,
        rt_length=30,
        target_strand=strand,
        editing_position=position,
        pe2_hek293t=score,
        pe2_k562=score * 0.9,
        pegrna_full="A" * 140,
    )


# ------------------------------------------------------------ constants


def test_min_nick_separation_is_23():
    assert MIN_NICK_SEPARATION_BP == 23


# --------------------------------------------------------- compatibility


def test_opposite_strand_rule():
    """Two Fw pegRNAs cannot pair."""
    up = _peg("u", "Fw", 50, 0.5)
    dn = _peg("d", "Fw", 100, 0.5)
    pairs = filter_compatible_pairs(upstream=[up], downstream=[dn], stop_codon_position=500)
    assert pairs == []


def test_cross_strand_pair_survives_if_distance_ok():
    up = _peg("u", "Fw", 50, 0.6)
    dn = _peg("d", "Rv", 100, 0.6)  # |100-50|=50, within [30,80]
    pairs = filter_compatible_pairs(upstream=[up], downstream=[dn], stop_codon_position=500)
    assert len(pairs) == 1
    assert isinstance(pairs[0], PegRNAPair)
    assert pairs[0].upstream is up
    assert pairs[0].downstream is dn


def test_distance_below_min_is_dropped():
    up = _peg("u", "Fw", 50, 0.6)
    dn = _peg("d", "Rv", 70, 0.6)  # distance 20 < 30
    pairs = filter_compatible_pairs(upstream=[up], downstream=[dn], stop_codon_position=500)
    assert pairs == []


def test_distance_above_max_is_dropped():
    up = _peg("u", "Fw", 50, 0.6)
    dn = _peg("d", "Rv", 200, 0.6)  # distance 150 > 80
    pairs = filter_compatible_pairs(upstream=[up], downstream=[dn], stop_codon_position=500)
    assert pairs == []


def test_min_nick_separation_enforced():
    """Even if min_distance would accept a pair, MIN_NICK_SEPARATION_BP does not."""
    up = _peg("u", "Fw", 50, 0.6)
    dn = _peg("d", "Rv", 60, 0.6)  # distance 10, below separation
    pairs = filter_compatible_pairs(
        upstream=[up],
        downstream=[dn],
        stop_codon_position=500,
        min_distance=5,  # caller would allow, but MIN_NICK_SEPARATION_BP blocks
        max_distance=80,
    )
    assert pairs == []


def test_custom_min_max_distance():
    up = _peg("u", "Fw", 50, 0.6)
    dn = _peg("d", "Rv", 200, 0.6)
    pairs = filter_compatible_pairs(
        upstream=[up],
        downstream=[dn],
        stop_codon_position=500,
        min_distance=100,
        max_distance=300,
    )
    assert len(pairs) == 1


# ------------------------------------------------------------- joint score


def test_joint_score_is_geometric_mean():
    up = _peg("u", "Fw", 50, 0.64)
    dn = _peg("d", "Rv", 100, 0.25)
    pairs = filter_compatible_pairs(upstream=[up], downstream=[dn], stop_codon_position=500)
    assert len(pairs) == 1
    expected = math.sqrt(0.64 * 0.25)
    assert pairs[0].joint_score == pytest.approx(expected)


def test_joint_score_zero_if_any_score_zero():
    up = _peg("u", "Fw", 50, 0.0)
    dn = _peg("d", "Rv", 100, 0.9)
    pairs = filter_compatible_pairs(upstream=[up], downstream=[dn], stop_codon_position=500)
    assert pairs[0].joint_score == 0.0


def test_joint_score_clamps_negatives_to_zero():
    """Defensive: PRIDICT sometimes emits small negative scores; we clamp at 0."""
    up = _peg("u", "Fw", 50, -0.1)
    dn = _peg("d", "Rv", 100, 0.9)
    pairs = filter_compatible_pairs(upstream=[up], downstream=[dn], stop_codon_position=500)
    assert pairs[0].joint_score == 0.0


# ------------------------------------------------------------- sorting


def test_pairs_returned_sorted_by_joint_score_descending():
    up1 = _peg("u1", "Fw", 50, 0.2)
    up2 = _peg("u2", "Fw", 50, 0.9)
    dn = _peg("d", "Rv", 100, 0.9)
    pairs = filter_compatible_pairs(upstream=[up1, up2], downstream=[dn], stop_codon_position=500)
    assert len(pairs) == 2
    assert pairs[0].joint_score >= pairs[1].joint_score
    assert pairs[0].upstream is up2  # higher-scoring pair first


def test_cartesian_product_explored():
    up1 = _peg("u1", "Fw", 50, 0.4)
    up2 = _peg("u2", "Fw", 60, 0.5)
    dn1 = _peg("d1", "Rv", 100, 0.6)
    dn2 = _peg("d2", "Rv", 110, 0.7)
    pairs = filter_compatible_pairs(
        upstream=[up1, up2], downstream=[dn1, dn2], stop_codon_position=500
    )
    assert len(pairs) == 4  # all 4 combos compatible


def test_nick_distance_preserved():
    up = _peg("u", "Fw", 50, 0.5)
    dn = _peg("d", "Rv", 95, 0.5)
    pairs = filter_compatible_pairs(upstream=[up], downstream=[dn], stop_codon_position=500)
    assert pairs[0].nick_distance == 45


def test_empty_inputs():
    assert filter_compatible_pairs(upstream=[], downstream=[], stop_codon_position=0) == []
