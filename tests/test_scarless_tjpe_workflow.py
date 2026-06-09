"""Tests for scarless tjPE workflow wiring (assembly, distance, ranking, gates)."""
from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.data_types import NgRNADesign
from pegrna_pipeline.hpc.tjpe_workflow import (
    assemble_tj_pegrna,
    combine_two_end_scores,
    require_b_at_s,
    select_tjpe_ngrna,
)


# ------------------------------------------------------------------ Task 5
def test_assemble_inserts_restored_segment_between_insert_and_pbs1():
    out = assemble_tj_pegrna(
        spacer="ACGTACGTACGTACGTACGT",
        insert_dna="ATGAAA",
        restored_segment="TACG",
        pbs1="CCCCGGGG",
        rc_pbs2="TTTAAAA",
    )
    assert "ATGAAA" + "TACG" + "CCCCGGGG" + "TTTAAAA" in out


def test_assemble_restored_empty_is_backwards_compatible():
    out = assemble_tj_pegrna(
        spacer="ACGTACGTACGTACGTACGT",
        insert_dna="ATGAAA",
        restored_segment="",
        pbs1="CCCCGGGG",
        rc_pbs2="TTTAAAA",
    )
    assert "ATGAAA" + "CCCCGGGG" + "TTTAAAA" in out


# ------------------------------------------------------------------ Task 6
def _ng(dist, strand="Rv", spacer="ACGTACGTACGTACGTACGT", score=50.0):
    return NgRNADesign(
        spacer=spacer,
        distance_to_edit=dist,
        target_strand=strand,
        deepcas9_score=score,
        pe3b=False,
        pam_disrupt=False,
    )


def test_select_ngrna_40bp_window():
    cands = [_ng(40), _ng(70), _ng(95)]
    chosen = select_tjpe_ngrna(cands, min_distance=36, max_distance=44, pegrna_strand="Fw")
    assert chosen is not None and abs(chosen.distance_to_edit) == 40


def test_select_ngrna_excludes_outside_window():
    cands = [_ng(70), _ng(95)]
    chosen = select_tjpe_ngrna(cands, min_distance=36, max_distance=44, pegrna_strand="Fw")
    assert chosen is None


def test_reverse_strand_requires_b_at_s():
    assert require_b_at_s(distance_to_edit=0, pegrna_strand="Rv") is True
    assert require_b_at_s(distance_to_edit=5, pegrna_strand="Rv") is False
    assert require_b_at_s(distance_to_edit=40, pegrna_strand="Fw") is True


def test_require_b_at_s_bad_strand():
    with pytest.raises(ValueError, match="pegrna_strand"):
        require_b_at_s(distance_to_edit=0, pegrna_strand="X")


# ------------------------------------------------------------------ Task 7
def test_combine_uses_minimum_bottleneck():
    assert combine_two_end_scores(insert_score=80.0, pbs2_score=30.0) == 30.0
    assert combine_two_end_scores(insert_score=10.0, pbs2_score=90.0) == 10.0
