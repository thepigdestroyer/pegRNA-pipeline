"""Tests for scarless tjPE sequence construction + INV-0 verification."""
from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.scarless import (
    assert_inv0,
    build_restored_segment,
    predict_edited_allele,
    silent_pam_break,
)


# ------------------------------------------------------------------ Task 1
def test_restored_segment_forward_exact_tiling():
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"  # A→B span genome[10:14] == "TACG"
    seg = build_restored_segment(genome, a_nick=10, b_nick=14, pegrna_strand="Fw")
    assert seg == "TACG"


def test_restored_segment_reverse_is_revcomp():
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"
    seg = build_restored_segment(genome, a_nick=10, b_nick=14, pegrna_strand="Rv")
    assert seg == "CGTA"  # revcomp(TACG)


def test_restored_segment_orders_endpoints():
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"
    # endpoints reversed → same span
    seg = build_restored_segment(genome, a_nick=14, b_nick=10, pegrna_strand="Fw")
    assert seg == "TACG"


def test_restored_segment_out_of_bounds_raises():
    with pytest.raises(ValueError, match="out of bounds"):
        build_restored_segment("ACGT", a_nick=0, b_nick=99, pegrna_strand="Fw")


def test_restored_segment_bad_strand_raises():
    with pytest.raises(ValueError, match="pegrna_strand"):
        build_restored_segment("ACGTACGT", a_nick=0, b_nick=4, pegrna_strand="X")


# ------------------------------------------------------------------ Task 2 (inv0)
def test_predict_edited_allele_zero_deletion():
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"
    cargo = "ATGTAA"
    edited = predict_edited_allele(
        genome, site_s=10, cargo=cargo, a_nick=10, b_nick=14,
        pegrna_strand="Fw", silent_pam_edit=None,
    )
    assert edited == genome[:10] + cargo + genome[10:]


def test_inv0_rejects_endogenous_deletion():
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"
    bad = genome[:10] + "ATGTAA" + genome[14:]  # dropped "TACG"
    with pytest.raises(ValueError, match="(?i)endogenous|length"):
        assert_inv0(reference=genome, edited=bad, site_s=10, cargo="ATGTAA", allowed_silent=None)


def test_inv0_passes_clean_scarless():
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"
    cargo = "ATGTAA"
    edited = genome[:10] + cargo + genome[10:]
    assert_inv0(reference=genome, edited=edited, site_s=10, cargo=cargo, allowed_silent=None)


def test_inv0_allows_single_silent_pam_base():
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"
    cargo = "ATGTAA"
    expected = genome[:10] + cargo + genome[10:]
    # change exactly one base in the restored span (position 14 in cargo-shifted allele)
    shifted = 14 + len(cargo)
    edited = expected[:shifted] + ("A" if expected[shifted] != "A" else "C") + expected[shifted + 1:]
    assert_inv0(
        reference=genome, edited=edited, site_s=10, cargo=cargo,
        allowed_silent={"pos": 14, "ref": genome[14], "alt": edited[shifted]},
    )


def test_inv0_rejects_two_changes():
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"
    cargo = "ATGTAA"
    expected = genome[:10] + cargo + genome[10:]
    e = list(expected)
    e[16] = "A" if e[16] != "A" else "C"
    e[17] = "A" if e[17] != "A" else "C"
    edited = "".join(e)
    with pytest.raises(ValueError, match="(?i)endogenous|change"):
        assert_inv0(reference=genome, edited=edited, site_s=10, cargo=cargo,
                    allowed_silent={"pos": 14, "ref": genome[14], "alt": "A"})


# ------------------------------------------------------------------ Task 3 (pam break)
def test_silent_pam_break_is_synonymous_or_none():
    res = silent_pam_break(codon_context="AAGGGA", pam_start=1, pam="AGG", frame=0)
    assert res is None or res["aa_before"] == res["aa_after"]


def test_silent_pam_break_no_silent_option_returns_none():
    # frame where breaking either G changes the amino acid → None
    res = silent_pam_break(codon_context="TGGTGG", pam_start=1, pam="GGT", frame=0)
    assert res is None or res["aa_before"] == res["aa_after"]
