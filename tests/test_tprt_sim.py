"""TPRT simulator: reconstruct the integrated SENSE allele from the designed
pegRNA extension and assert it equals the intended product.

The binding guarantee is sense-allele == intended product. The NEGATIVE tests
are the point: feeding the wrong strand operation (plain ``complement`` or plain
``reverse`` instead of ``reverse_complement``) MUST be rejected by the simulator,
proving the pipeline enforces the three-way distinction.
"""

from Bio.Seq import Seq

import pytest

from pegrna_pipeline.hpc.strand_orientation import orient_rtt_tail
from pegrna_pipeline.hpc.tprt_sim import TPRTMismatch, simulate_tprt


def _complement(s: str) -> str:
    return str(Seq(s).complement())


def _reverse(s: str) -> str:
    return s[::-1]


def _rc(s: str) -> str:
    return str(Seq(s).reverse_complement())


# Geometry: sense genome, S at 30, A->B arm [30,42). Cargo inserted at S.
_GENOME = "G" * 30 + "ATGCATGCATGC" + "C" * 30
_SITE_S = 30
_HI = 42
_CARGO = "GATTACAGATTACA"
_ARM = _GENOME[_SITE_S:_HI]  # "ATGCATGCATGC" (no recode for this test)


def _intended():
    return _GENOME[:_SITE_S] + _CARGO + _ARM + _GENOME[_HI:]


def test_correct_orientation_reconstructs_sense_allele_reverse():
    sense_flap = _CARGO + _ARM
    rtt_tail = orient_rtt_tail(sense_flap, pegrna_strand="Rv")  # reverse_complement
    allele = simulate_tprt(
        genome=_GENOME,
        site_s=_SITE_S,
        hi=_HI,
        rtt_tail_on_peg=rtt_tail,
        pegrna_strand="Rv",
    )
    assert allele == _intended()


def test_correct_orientation_reconstructs_sense_allele_forward():
    sense_flap = _CARGO + _ARM
    rtt_tail = orient_rtt_tail(sense_flap, pegrna_strand="Fw")  # literal
    allele = simulate_tprt(
        genome=_GENOME,
        site_s=_SITE_S,
        hi=_HI,
        rtt_tail_on_peg=rtt_tail,
        pegrna_strand="Fw",
    )
    assert allele == _intended()


def test_rejects_plain_complement_instead_of_reverse_complement():
    sense_flap = _CARGO + _ARM
    wrong = _complement(sense_flap)  # complement, NOT reverse_complement
    with pytest.raises(TPRTMismatch):
        simulate_tprt(
            genome=_GENOME,
            site_s=_SITE_S,
            hi=_HI,
            rtt_tail_on_peg=wrong,
            pegrna_strand="Rv",
            expected_allele=_intended(),
        )


def test_rejects_plain_reverse_instead_of_reverse_complement():
    sense_flap = _CARGO + _ARM
    wrong = _reverse(sense_flap)  # reverse, NOT reverse_complement
    with pytest.raises(TPRTMismatch):
        simulate_tprt(
            genome=_GENOME,
            site_s=_SITE_S,
            hi=_HI,
            rtt_tail_on_peg=wrong,
            pegrna_strand="Rv",
            expected_allele=_intended(),
        )


def test_expected_allele_guard_passes_on_correct_design():
    sense_flap = _CARGO + _ARM
    rtt_tail = orient_rtt_tail(sense_flap, pegrna_strand="Rv")
    # When expected_allele is supplied and matches, no raise and returns it.
    allele = simulate_tprt(
        genome=_GENOME,
        site_s=_SITE_S,
        hi=_HI,
        rtt_tail_on_peg=rtt_tail,
        pegrna_strand="Rv",
        expected_allele=_intended(),
    )
    assert allele == _intended()
