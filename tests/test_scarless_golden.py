"""Golden-reference INV-0 checks: zero endogenous deletion on both strands.

Synthetic stand-ins for AAVS1/HEK3/Fah-style geometry; real coordinates from
Zheng 2023 Supplementary Table 1 can be substituted later. The invariant under
test is the make-or-break one: the predicted SENSE allele equals reference with
cargo inserted at S and nothing endogenous removed, regardless of pegRNA strand.
"""

from __future__ import annotations

from pegrna_pipeline.hpc.scarless import (
    assert_inv0,
    build_restored_segment,
    predict_edited_allele,
)


def test_golden_forward_zero_deletion():
    genome = "G" * 30 + "ATGCATGCATGC" + "C" * 30  # S at 30, A→B span = [30,42)
    cargo = "GATTACAGATTACA"
    edited = predict_edited_allele(
        genome,
        site_s=30,
        cargo=cargo,
        a_nick=30,
        b_nick=42,
        pegrna_strand="Fw",
        silent_pam_edit=None,
    )
    assert_inv0(reference=genome, edited=edited, site_s=30, cargo=cargo, allowed_silent=None)
    assert len(edited) == len(genome) + len(cargo)  # zero deletion


def test_golden_reverse_revcomp_cargo_zero_deletion():
    genome = "G" * 30 + "ATGCATGCATGC" + "C" * 30
    cargo = "GATTACAGATTACA"
    edited = predict_edited_allele(
        genome,
        site_s=30,
        cargo=cargo,
        a_nick=30,
        b_nick=42,
        pegrna_strand="Rv",
        silent_pam_edit=None,
    )
    # Sense allele identical regardless of which strand was synthesised.
    assert_inv0(reference=genome, edited=edited, site_s=30, cargo=cargo, allowed_silent=None)
    assert len(edited) == len(genome) + len(cargo)


def test_golden_restored_segment_reconstitutes_span_forward():
    genome = "G" * 30 + "ATGCATGCATGC" + "C" * 30
    seg = build_restored_segment(genome, a_nick=30, b_nick=42, pegrna_strand="Fw")
    assert seg == "ATGCATGCATGC"  # exact A→B span, in sense


def test_golden_restored_segment_revcomp_reverse():
    from Bio.Seq import Seq

    genome = "G" * 30 + "ATGCATGCATGC" + "C" * 30
    seg = build_restored_segment(genome, a_nick=30, b_nick=42, pegrna_strand="Rv")
    assert seg == str(Seq("ATGCATGCATGC").reverse_complement())


def test_golden_inv0_catches_one_base_endogenous_change_without_permission():
    genome = "G" * 30 + "ATGCATGCATGC" + "C" * 30
    cargo = "GATTACAGATTACA"
    expected = genome[:30] + cargo + genome[30:]
    # flip one endogenous base in the restored span; no silent edit declared → reject
    pos = 30 + len(cargo) + 2
    edited = expected[:pos] + ("A" if expected[pos] != "A" else "T") + expected[pos + 1 :]
    import pytest

    with pytest.raises(ValueError, match="(?i)endogenous"):
        assert_inv0(reference=genome, edited=edited, site_s=30, cargo=cargo, allowed_silent=None)
