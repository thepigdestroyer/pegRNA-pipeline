"""Tests for :func:`select_tjpe_ngrna`.

Filters:
- ``|distance_to_edit|`` falls in [min_distance, max_distance].
- If ``pegrna_strand`` is supplied, keep only opposite-strand candidates
  (mapping {Fw: Rv, Rv: Fw}).
- Sanity filters on the spacer (hard rules, see ``guide_sanity``): no TTTT (RNA-pol III
  terminator), GC content in [0.30, 0.75]. Homopolymer runs are advisory, not filtered.

Returns the candidate with the maximum ``deepcas9_score`` that survives,
or ``None`` if no candidate survives.
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.data_types import NgRNADesign
from pegrna_pipeline.hpc.tjpe_workflow import select_tjpe_ngrna


def _ng(
    spacer: str = "AGTGCTACCTGGTACGGATC",
    distance: int = 50,
    strand: str = "Rv",
    score: float = 0.8,
) -> NgRNADesign:
    return NgRNADesign(
        spacer=spacer,
        distance_to_edit=distance,
        target_strand=strand,
        deepcas9_score=score,
        pe3b=False,
        pam_disrupt=False,
    )


# ---------------------------------------------------------- distance window


def test_distance_within_window_is_kept():
    # default window = [40, 100]
    assert select_tjpe_ngrna([_ng(distance=50)]) is not None


def test_distance_below_window_is_dropped():
    assert select_tjpe_ngrna([_ng(distance=30)]) is None


def test_distance_above_window_is_dropped():
    assert select_tjpe_ngrna([_ng(distance=150)]) is None


def test_negative_distance_respects_absolute_value():
    """|-50| = 50, within [40, 100]; keep it."""
    assert select_tjpe_ngrna([_ng(distance=-50)]) is not None
    # |-30| = 30, outside window; drop.
    assert select_tjpe_ngrna([_ng(distance=-30)]) is None


def test_custom_min_max_distance():
    cands = [_ng(distance=20), _ng(distance=30), _ng(distance=110)]
    out = select_tjpe_ngrna(cands, min_distance=20, max_distance=100)
    # distance=110 is excluded; two survive — picks highest score (tied → first).
    assert out is not None
    assert abs(out.distance_to_edit) <= 100


# ---------------------------------------------------------------- strand


def test_opposite_strand_filter_applied():
    cands = [_ng(strand="Fw", distance=50, score=0.9), _ng(strand="Rv", distance=50, score=0.1)]
    # pegrna_strand="Fw" means ngRNA must be "Rv"; the 0.1 Rv one wins.
    out = select_tjpe_ngrna(cands, pegrna_strand="Fw")
    assert out is not None
    assert out.target_strand == "Rv"
    assert out.deepcas9_score == pytest.approx(0.1)


def test_no_strand_filter_when_pegrna_strand_none():
    cands = [_ng(strand="Fw", distance=50, score=0.9), _ng(strand="Rv", distance=50, score=0.1)]
    out = select_tjpe_ngrna(cands)  # pegrna_strand=None
    assert out is not None
    assert out.deepcas9_score == pytest.approx(0.9)


# -------------------------------------------------------------- sanity filters


def test_reject_spacer_with_tttt():
    cands = [_ng(spacer="ATTTTACCTGGTACGGATCA", score=0.9)]  # has TTTT
    assert select_tjpe_ngrna(cands) is None


def test_long_homopolymer_is_not_a_hard_rejection():
    """Only TTTT (a Pol III stop signal) has a mechanism, so a run of A, C or G is reported as
    an advisory instead (guide_sanity); it used to be rejected here, on this path only, which made
    the same guide acceptable to the genome-enumerated path and not to this one."""
    cands = [_ng(spacer="AAAAAACCTGGTACGGATCA", score=0.9)]  # 6-A run
    assert select_tjpe_ngrna(cands) is not None


def test_reject_spacer_with_low_gc():
    # ~10% GC
    cands = [_ng(spacer="AAAATAATAATAAAATAAAT", score=0.9)]
    assert select_tjpe_ngrna(cands) is None


def test_reject_spacer_with_high_gc():
    # 100% GC
    cands = [_ng(spacer="GCGCGCGCGCGCGCGCGCGC", score=0.9)]
    assert select_tjpe_ngrna(cands) is None


# ---------------------------------------------------------- selection logic


def test_picks_max_deepcas9_score_among_survivors():
    cands = [
        _ng(distance=50, score=0.2),
        _ng(distance=50, score=0.9),
        _ng(distance=50, score=0.5),
    ]
    out = select_tjpe_ngrna(cands)
    assert out is not None
    assert out.deepcas9_score == pytest.approx(0.9)


def test_returns_none_when_all_filtered():
    cands = [_ng(distance=10), _ng(distance=200)]
    assert select_tjpe_ngrna(cands) is None


def test_empty_candidate_list():
    assert select_tjpe_ngrna([]) is None
