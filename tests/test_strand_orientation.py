"""Strand orientation for manually-designed pegRNA extension elements.

Convention (LOCKED 2026-06-19): everything is stored 5'->3'; antiparallel /
templating relationships are ``reverse_complement``. The emitted RTT-tail in the
final 5'->3' pegRNA oligo is the reverse complement of the desired sense flap.

rcPBS2 (ADR 0014, target-strand convention): the tj-pegRNA carries ``rcPBS2 =
RC(TS)`` — the reverse complement of the target-strand window it captures. PBS2
itself lives on the synthesised flap, not the pegRNA; the old strand-dependent
``orient_pbs2`` (TS-literal for a reverse-targeting pegRNA) was a bug and is now
a deprecated alias.
"""

import pytest
from Bio.Seq import Seq

from pegrna_pipeline.hpc.strand_orientation import (
    orient_pbs2,
    orient_rcpbs2,
    orient_rtt_tail,
)


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
    with pytest.raises(ValueError, match="strand"):
        orient_rtt_tail("ACGT", pegrna_strand="X")


def test_orient_rcpbs2_is_reverse_complement_of_target_strand_window():
    # rcPBS2 carried by the tj-pegRNA is RC of the target-strand window it captures
    # (Zheng 2023; ADR 0014). Non-palindromic worked value proves it is not identity.
    ts_window = "GCATGACTGA"
    assert orient_rcpbs2(ts_window) == _rc(ts_window) == "TCAGTCATGC"


def test_orient_pbs2_is_deprecated_alias_returning_rcpbs2():
    # The old strand-dependent orient_pbs2 (TS-literal for Rv) was a bug (D4); it is
    # now a deprecated shim that warns and returns rcPBS2 = RC(TS).
    ts_window = "GCATGACTGA"
    with pytest.warns(DeprecationWarning):
        assert orient_pbs2(ts_window) == orient_rcpbs2(ts_window)
    # Backward-compatible signature: the old pegrna_strand kwarg is accepted + ignored.
    with pytest.warns(DeprecationWarning):
        assert orient_pbs2(ts_window, pegrna_strand="Rv") == orient_rcpbs2(ts_window)


def test_orientations_are_inverse_consistent():
    # Re-orienting twice for the same strand round-trips the sequence.
    seq = "ATGCAATTGGCC"
    assert orient_rtt_tail(orient_rtt_tail(seq, pegrna_strand="Rv"), pegrna_strand="Rv") == seq
