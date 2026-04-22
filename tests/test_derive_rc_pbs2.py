"""Tests for :func:`derive_rc_pbs2`.

PBS2 (the PBS of the second/downstream pegRNA) is pulled from the **top
strand** of the mRNA ending at ``ngrna_nick_position`` and then
reverse-complemented (because the pegRNA anneals to the opposite strand).

    top  -----<-nick
    ...----PBS2---|-----...   (top strand; reads 5' -> 3')
                  ^ nick

``rc_pbs2 = revcomp( top[nick - pbs2_length : nick] )``
"""

from __future__ import annotations

import pytest
from Bio.Seq import Seq

from pegrna_pipeline.hpc.tjpe_workflow import derive_rc_pbs2


def test_rc_pbs2_is_reverse_complement_of_window():
    top = "ACGTACGTACGTACGTACGTACGTACGT"
    nick = 20
    length = 5
    window = top[nick - length : nick]
    rc = str(Seq(window).reverse_complement())
    assert derive_rc_pbs2(top, ngrna_nick_position=nick, pbs2_length=length) == rc


def test_default_pbs2_length_is_17():
    top = "A" * 100
    nick = 50
    out = derive_rc_pbs2(top, ngrna_nick_position=nick)
    assert len(out) == 17
    # A reverse-complements to T, so the RC of 17 As is 17 Ts.
    assert out == "T" * 17


def test_negative_start_raises():
    """If ``nick < pbs2_length`` the window would have negative start."""
    with pytest.raises(ValueError):
        derive_rc_pbs2("ACGTACGTAC", ngrna_nick_position=5, pbs2_length=10)


def test_nick_past_sequence_raises():
    with pytest.raises(ValueError):
        derive_rc_pbs2("ACGTACGT", ngrna_nick_position=100, pbs2_length=5)


def test_output_is_acgt_only():
    top = "ACGTACGTACGTACGTACGTACGTACGT"
    out = derive_rc_pbs2(top, ngrna_nick_position=20, pbs2_length=17)
    assert set(out) <= set("ACGT")


def test_custom_length():
    top = "AAAAAAAAAAACCCCCCCCCCCC"
    # Window = top[5:10] = "AAAAA"; rc = "TTTTT"
    out = derive_rc_pbs2(top, ngrna_nick_position=10, pbs2_length=5)
    assert out == "TTTTT"
    # Window = top[0:11] = "AAAAAAAAAAA"; rc = "T"*11
    out = derive_rc_pbs2(top, ngrna_nick_position=11, pbs2_length=11)
    assert out == "T" * 11


def test_length_equal_to_nick_allowed():
    """Edge: nick == pbs2_length is valid (window starts at 0)."""
    top = "ACGTACGT"
    out = derive_rc_pbs2(top, ngrna_nick_position=4, pbs2_length=4)
    assert out == str(Seq("ACGT").reverse_complement())
