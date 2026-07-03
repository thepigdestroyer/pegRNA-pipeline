"""Stage C (C1/C3/C4): two-pegRNA PRIDICT adapter + orientation + ngRNA lifecycle.

The tj-pegRNA is assembled from TWO PRIDICT designs (ADR 0014/0016):
- pegRNA#1 (designed on the target strand) -> PBS1 + first insert end, used UNTOUCHED
  (PRIDICT already emits them as RC(TS) via PBSrevcomp/RTrevcomp).
- pegRNA#2 (designed on the OPPOSITE strand) -> rcPBS2 = RC(its PBSrevcomp),
  last_end = RC(its RTrevcomp); its spacer+scaffold BECOMES the ngRNA, bound to the design.
So every 3'-extension component ends up RC(TS).
"""

from __future__ import annotations

import dataclasses

import pytest
from Bio.Seq import Seq

from pegrna_pipeline.hpc.data_types import PegRNADesign
from pegrna_pipeline.hpc.recoded_pridict import (
    TwoPegDesign,
    TwoPegPridictFake,
    combine_two_peg,
)

SCAF = "GTTTTAGAGCTAGAAATAGCAAG"  # stand-in scaffold for tests


def _rc(s: str) -> str:
    return str(Seq(s).reverse_complement())


def _peg(name: str, spacer: str, pbs: str, rt: str, ts: str = "Rv") -> PegRNADesign:
    return PegRNADesign(
        sequence_name=name,
        spacer=spacer,
        pbs=pbs,
        rt=rt,
        pbs_length=len(pbs),
        rt_length=len(rt),
        target_strand=ts,
        editing_position=5,
        pe2_hek293t=50.0,
        pe2_k562=50.0,
        pegrna_full="",
    )


# --------------------------------------------------------------- C1: dataclass + fake
def test_twopegdesign_offtarget_risk_defaults_none():
    d = TwoPegDesign(
        pbs1="AAA",
        first_end="CCC",
        rcpbs2="GGG",
        last_end="TTT",
        ngrna="ACGT" + SCAF,
        pbs1_len=3,
        rt_len=3,
        editing_position=5,
        target_strand="Rv",
    )
    # PRIDICT2 exposes no per-spacer off-target-risk field (ADR 0016) -> None.
    assert d.offtarget_risk is None


def test_fake_emits_pridict_orientation():
    # All window args are target-strand. Real PRIDICT emits RC(TS): pegRNA#1 (on TS) ->
    # RC(TS) windows; pegRNA#2 (on OS) -> PBSrevcomp=RC(OS)=TS literally. The fake mirrors that.
    peg1, peg2 = TwoPegPridictFake().design_pair(
        first_end="ACGTACGT",
        last_end="TTTTGGGG",
        pbs1_window="CCCCAAAA",
        pbs2_window="GGGGTTTT",
        spacer1="G" * 20,
        spacer2="A" * 20,
        target_strand="Rv",
    )
    assert isinstance(peg1, PegRNADesign) and isinstance(peg2, PegRNADesign)
    assert peg1.pbs == _rc("CCCCAAAA") and peg1.rt == _rc("ACGTACGT")  # peg1 on TS -> RC(TS)
    assert peg2.pbs == "GGGGTTTT" and peg2.rt == "TTTTGGGG"  # peg2 on OS -> TS literally


def test_fake_then_combine_is_rc_of_target_strand_end_to_end():
    # The CI path (fake -> combine) must yield RC(TS) for every component starting from
    # target-strand windows -- guards a future fake regression that silently inverts a strand.
    kw = dict(
        first_end="ACGTACGTAC",
        last_end="TTTTGGGGAC",
        pbs1_window="CCCCAAAAGT",
        pbs2_window="GGGGTTTTCA",
        spacer1="G" * 20,
        spacer2="A" * 20,
        target_strand="Rv",
    )
    peg1, peg2 = TwoPegPridictFake().design_pair(**kw)
    d = combine_two_peg(peg1, peg2, scaffold=SCAF)
    assert d.pbs1 == _rc(kw["pbs1_window"])
    assert d.first_end == _rc(kw["first_end"])
    assert d.rcpbs2 == _rc(kw["pbs2_window"])
    assert d.last_end == _rc(kw["last_end"])
    assert d.rcpbs2 != kw["pbs2_window"]  # the RC is not a no-op (non-palindromic window)


# ---------------------------------------------------- C3: orientation = RC(TS) end-to-end
def test_combine_two_peg_is_rc_of_target_strand_end_to_end():
    ts_pbs1, ts_first, ts_pbs2, ts_last = "CCCCAAAAGT", "ACGTACGTAC", "GGGGTTTTCA", "TTTTGGGGAC"
    # pegRNA#1 designed on TS: PRIDICT emits RC(TS windows).
    peg1 = _peg("p1", "G" * 20, pbs=_rc(ts_pbs1), rt=_rc(ts_first))
    # pegRNA#2 designed on the OPPOSITE strand: its revcomp fields == the TS windows literally.
    peg2 = _peg("p2", "A" * 20, pbs=ts_pbs2, rt=ts_last)
    d = combine_two_peg(peg1, peg2, scaffold=SCAF)
    assert d.pbs1 == _rc(ts_pbs1)  # untouched (already RC(TS))
    assert d.first_end == _rc(ts_first)  # untouched (already RC(TS))
    assert d.rcpbs2 == _rc(ts_pbs2)  # RC of peg2's PBSrevcomp -> RC(TS)
    assert d.last_end == _rc(ts_last)  # RC of peg2's RTrevcomp -> RC(TS)
    # Every extension component is RC of its target-strand window:
    for comp, ts_win in [
        (d.pbs1, ts_pbs1),
        (d.first_end, ts_first),
        (d.rcpbs2, ts_pbs2),
        (d.last_end, ts_last),
    ]:
        assert comp == _rc(ts_win)
    # Non-tautology guard: the RC on peg2 is actually applied (would fail if combine
    # dropped the RC, since the windows are non-palindromic).
    assert d.rcpbs2 != peg2.pbs
    assert d.last_end != peg2.rt


def test_combine_pbs1_and_first_end_are_verbatim_from_peg1():
    peg1 = _peg("p1", "G" * 20, pbs="ACGTACGT", rt="TTTTCCCCAA")
    peg2 = _peg("p2", "A" * 20, pbs="GGGGAAAA", rt="CCCCGGGG")
    d = combine_two_peg(peg1, peg2, scaffold=SCAF)
    assert d.pbs1 == peg1.pbs and d.first_end == peg1.rt  # verbatim, not re-oriented
    assert d.pbs1_len == peg1.pbs_length and d.rt_len == peg1.rt_length


# ----------------------------------------------------------------- C4: ngRNA lifecycle
def test_ngrna_bound_as_peg2_spacer_plus_scaffold():
    peg1 = _peg("p1", "G" * 20, pbs="ACGT", rt="TTTT")
    peg2 = _peg("p2", "ACGTACGTACGTACGTACGT", pbs="GGGG", rt="CCCC")
    d = combine_two_peg(peg1, peg2, scaffold=SCAF)
    assert d.ngrna == peg2.spacer + SCAF  # ngRNA = pegRNA#2 spacer + scaffold, untouched


def test_ngrna_rebinds_on_new_peg2():
    # The ngRNA is bound to (never independent of) its tj-pegRNA: a new peg2 -> a new design
    # with a new ngRNA; discarding the design discards the ngRNA (they are one object).
    peg1 = _peg("p1", "G" * 20, pbs="ACGT", rt="TTTT")
    da = combine_two_peg(peg1, _peg("p2a", "A" * 20, pbs="GGGG", rt="CCCC"), scaffold=SCAF)
    db = combine_two_peg(peg1, _peg("p2b", "C" * 20, pbs="GGGG", rt="CCCC"), scaffold=SCAF)
    assert da.ngrna != db.ngrna
    assert db.ngrna == "C" * 20 + SCAF


# ---------------------------------------------------------------- misc invariants
def test_combine_propagates_peg1_edit_metadata_and_offtarget_risk():
    peg1 = _peg("p1", "G" * 20, pbs="ACGT", rt="TTTT", ts="Fw")
    peg2 = _peg("p2", "A" * 20, pbs="GGGG", rt="CCCC")
    d = combine_two_peg(peg1, peg2, scaffold=SCAF, offtarget_risk=0.7)
    assert d.editing_position == peg1.editing_position  # peg1 authoritative for edit metadata
    assert d.target_strand == peg1.target_strand == "Fw"
    assert d.offtarget_risk == 0.7  # passthrough when a caller supplies one


def test_twopegdesign_is_frozen():
    d = combine_two_peg(
        _peg("p1", "G" * 20, pbs="ACGT", rt="TTTT"),
        _peg("p2", "A" * 20, pbs="GGGG", rt="CCCC"),
        scaffold=SCAF,
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        d.rcpbs2 = "AAAA"  # type: ignore[misc]
