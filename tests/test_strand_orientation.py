"""Strand orientation for manually-designed pegRNA extension elements.

Convention (LOCKED 2026-06-19): everything is stored 5'->3'; antiparallel /
templating relationships are ``reverse_complement``. The emitted RTT-tail in the
final 5'->3' pegRNA oligo is the reverse complement of the desired sense flap.
"""

from Bio.Seq import Seq

from pegrna_pipeline.hpc.strand_orientation import orient_pbs2, orient_rtt_tail


def _rc(s: str) -> str:
    return str(Seq(s).reverse_complement())


def test_rtt_tail_reverse_strand_is_reverse_complement_of_sense_flap():
    sense_flap = "ATGC"
    out = orient_rtt_tail(sense_flap, pegrna_strand="Rv")
    assert out == _rc(sense_flap) == "GCAT"


def test_rtt_tail_forward_strand_is_literal():
    # A forward-targeting pegRNA templates the sense flap directly.
    sense_flap = "ATGCAAA"
    assert orient_rtt_tail(sense_flap, pegrna_strand="Fw") == sense_flap


def test_rtt_tail_rejects_bad_strand():
    import pytest

    with pytest.raises(ValueError, match="strand"):
        orient_rtt_tail("ACGT", pegrna_strand="X")


def test_pbs2_reverse_strand_identity_to_captured_strand():
    # PBS2-on-pegRNA must equal the captured-strand segment (complement-of-
    # complement = identity). For a reverse-targeting pegRNA the captured strand
    # is the opposite (sense top) strand window, taken literally.
    captured = "GCATGACTGA"
    assert orient_pbs2(captured, pegrna_strand="Rv") == captured


def test_pbs2_forward_strand_is_reverse_complement():
    captured = "GCATGACTGA"
    assert orient_pbs2(captured, pegrna_strand="Fw") == _rc(captured)


def test_orientations_are_inverse_consistent():
    # Re-orienting twice for the same strand round-trips the sequence.
    seq = "ATGCAATTGGCC"
    assert orient_rtt_tail(orient_rtt_tail(seq, pegrna_strand="Rv"), pegrna_strand="Rv") == seq
