"""Span-aware INV-0 (Task 2): recode_span (synonymous OK), <=2 silent PAM edits,
locked_span (non-coding byte-exact), and legacy backward-compat."""

import pytest

from pegrna_pipeline.hpc.scarless import assert_inv0


def _silent_ok(_ref, _edited):
    return True


def _silent_bad(_ref, _edited):
    return False


def test_inv0_allows_many_synonymous_changes_inside_recode_span():
    ref = "ATG" "AAA" "CCC" "GGG" "TAA"  # 15 nt
    cargo = "TTTTTT"
    # cargo inserted at S=12 (before stop); recode span [3,12) recoded synonymously
    edited = "ATG" "AAG" "CCA" "GGA" "TTTTTT" "TAA"  # 3 changes inside [3,12)
    assert_inv0(
        reference=ref, edited=edited, site_s=12, cargo=cargo,
        recode_span=(3, 12), protein_check=_silent_ok, allowed_silent=None,
    )


def test_inv0_rejects_change_outside_recode_span():
    ref = "ATGAAACCCGGGTAA"
    cargo = "TTTTTT"
    # change at pos 2 (T->C) is OUTSIDE recode span [3,12) and not an allowed silent
    edited = "ATC" "AAACCCGGG" "TTTTTT" "TAA"
    with pytest.raises(ValueError, match="outside"):
        assert_inv0(
            reference=ref, edited=edited, site_s=12, cargo=cargo,
            recode_span=(3, 12), protein_check=_silent_ok, allowed_silent=None,
        )


def test_inv0_recode_span_changes_must_be_protein_silent():
    ref = "ATGAAACCCGGGTAA"
    cargo = "TTTTTT"
    edited = "ATG" "AAG" "CCA" "GGA" "TTTTTT" "TAA"
    with pytest.raises(ValueError, match="not protein-silent"):
        assert_inv0(
            reference=ref, edited=edited, site_s=12, cargo=cargo,
            recode_span=(3, 12), protein_check=_silent_bad, allowed_silent=None,
        )


# Reference layout for span tests (S=18, cargo inserts at the very end before stop):
#   idx 0..17 : ATG AAA CCC GGG AAA CCC   (CDS, 18 nt)  then cargo, then TAA tail
_REF = "ATGAAACCCGGGAAACCC" "TAA"  # 21 nt; S=18 inserts cargo just before the stop
_CARGO = "TTT"


def _expected():
    return _REF[:18] + _CARGO + _REF[18:]


def test_inv0_allows_up_to_two_silent_pam_edits_outside_span():
    # B2-ii: two single-base silent PAM breaks OUTSIDE the recode span [9,18) are allowed
    exp = list(_expected())
    exp[1] = "C"  # change 1 (outside span, pos 1)
    exp[7] = "A"  # change 2 (outside span, pos 7)
    edited = "".join(exp)
    assert_inv0(
        reference=_REF, edited=edited, site_s=18, cargo=_CARGO,
        recode_span=(9, 18), protein_check=_silent_ok,
        allowed_silent=[{"pos": 1}, {"pos": 7}],
    )


def test_inv0_rejects_third_endogenous_change():
    exp = list(_expected())
    exp[1] = "C"
    exp[7] = "A"
    exp[4] = "G"  # third change outside span (ref[4]=A->G), only 2 allowed
    edited = "".join(exp)
    with pytest.raises(ValueError, match="outside recode_span"):
        assert_inv0(
            reference=_REF, edited=edited, site_s=18, cargo=_CARGO,
            recode_span=(9, 18), protein_check=_silent_ok,
            allowed_silent=[{"pos": 1}, {"pos": 7}],
        )


def test_inv0_locked_span_must_match_targeted_allele_exactly():
    # decision 15: recode [3,9) CDS ok; locked [9,18) non-coding must be byte-exact
    exp = list(_expected())
    exp[5] = "G"  # synonymous-ish change inside recode span [3,9) — allowed
    exp[10] = "T"  # change inside locked span [9,18) — MUST raise
    edited = "".join(exp)
    with pytest.raises(ValueError, match="non-coding locked span"):
        assert_inv0(
            reference=_REF, edited=edited, site_s=18, cargo=_CARGO,
            recode_span=(3, 9), protein_check=_silent_ok, allowed_silent=None,
            locked_span=(9, 18), locked_expected=_REF[9:18],
        )


def test_inv0_locked_span_passes_when_exact():
    exp = list(_expected())
    exp[5] = "G"  # change inside recode span [3,9) only
    edited = "".join(exp)
    assert_inv0(
        reference=_REF, edited=edited, site_s=18, cargo=_CARGO,
        recode_span=(3, 9), protein_check=_silent_ok, allowed_silent=None,
        locked_span=(9, 18), locked_expected=_REF[9:18],
    )


# ---- legacy backward-compat (recode_span=None) ----

def test_inv0_legacy_clean_insertion_passes():
    ref = "ATGAAACCCGGGTAA"
    cargo = "TTT"
    edited = "ATGAAACCC" "TTT" "GGGTAA"
    assert_inv0(reference=ref, edited=edited, site_s=9, cargo=cargo, allowed_silent=None)


def test_inv0_legacy_single_silent_dict_still_allowed():
    ref = "ATGAAACCCGGGTAA"
    cargo = "TTT"
    edited = "ATC" "AAACCC" "TTT" "GGGTAA"  # one change at pos 2
    assert_inv0(
        reference=ref, edited=edited, site_s=9, cargo=cargo,
        allowed_silent={"pos": 2},
    )


def test_inv0_legacy_rejects_change_with_no_silent():
    ref = "ATGAAACCCGGGTAA"
    cargo = "TTT"
    edited = "ATC" "AAACCC" "TTT" "GGGTAA"
    with pytest.raises(ValueError):
        assert_inv0(reference=ref, edited=edited, site_s=9, cargo=cargo, allowed_silent=None)
