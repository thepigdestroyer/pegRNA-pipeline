"""Regression tests for the rcPBS2 length-selection + orientation.

The old design fixed PBS2 length to PBS1's lengths (13/15/17) and then rejected on
Tm. Because rcPBS2 is genome-determined at the nick and the 3 lengths are nested
windows there, the feasible set collapsed and whole candidate pools failed the gate
(the CTCF `no_valid_pegrna` symptom). ``_select_rcpbs2`` instead scans the length
range (10-18) per-nick to hit the Tm band, decoupled from PBS1.

Orientation (ADR 0014): the tj-pegRNA carries ``rcPBS2 = RC(TS)``. For a
reverse-targeting (``Rv``) pegRNA the target strand is the sense/top window, so the
returned rcPBS2 is the **reverse complement** of the captured sense window — NOT the
literal window (the old strand-dependent ``orient_pbs2`` bug). Tm(NN), GC and the
repetitiveness screen are reverse-complement-invariant, so the selected LENGTH is
unchanged by the orientation fix; only the emitted bases change.
"""

from __future__ import annotations

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

from pegrna_pipeline.hpc.gate_thresholds import THRESHOLDS
from pegrna_pipeline.hpc.recoded_tjpe_workflow import _is_repetitive, _select_rcpbs2

# Window (sense/top strand) where the OLD fixed lengths all miss the band: len13
# Tm=36.5 (<37), len15 Tm=45.7 (>42), len17 Tm=51.1 (>42); but a length-14 window
# lands at Tm=40.8, in [37,42]. The nick sits at the window's end.
_MIXED = "TGTTGGCCCAGTGTGAATCG"  # 20 nt, ends at the nick


def _rc(s: str) -> str:
    return str(Seq(s).reverse_complement())


def _band():
    return dict(
        tm_lo=THRESHOLDS.pbs2_tm_lo,
        tm_hi=THRESHOLDS.pbs2_tm_hi,
        gc_lo=THRESHOLDS.pbs2_gc_lo,
        gc_hi=THRESHOLDS.pbs2_gc_hi,
        len_min=THRESHOLDS.pbs2_len_min,
        len_max=THRESHOLDS.pbs2_len_max,
    )


def test_select_rcpbs2_scan_succeeds_where_fixed_lengths_all_fail():
    # Sanity: the three OLD fixed lengths really do all miss the band here.
    for length, expect_in in ((13, False), (15, False), (17, False)):
        s = _MIXED[-length:]
        in_band = THRESHOLDS.pbs2_tm_lo <= mt.Tm_NN(s) <= THRESHOLDS.pbs2_tm_hi
        assert in_band is expect_in
    # The scan finds a valid rcPBS2 (length 14) where the fixed approach returned none.
    rcpbs2, why = _select_rcpbs2(_MIXED, nick_position=len(_MIXED), pegrna_strand="Fw", **_band())
    assert rcpbs2 is not None, why
    assert len(rcpbs2) == 14  # the length the fixed (13/15/17) set could never pick
    tm = mt.Tm_NN(rcpbs2)
    assert THRESHOLDS.pbs2_tm_lo <= tm <= THRESHOLDS.pbs2_tm_hi


def test_select_rcpbs2_returns_tm_closest_to_band_centre():
    rcpbs2, _why = _select_rcpbs2(_MIXED, nick_position=len(_MIXED), pegrna_strand="Fw", **_band())
    tm = mt.Tm_NN(rcpbs2)
    center = (THRESHOLDS.pbs2_tm_lo + THRESHOLDS.pbs2_tm_hi) / 2
    # No other in-band, GC/repeat-valid length is closer to the centre than the pick
    # (length choice is RC-invariant, so compare against the raw sense windows).
    for length in range(THRESHOLDS.pbs2_len_min, THRESHOLDS.pbs2_len_max + 1):
        s = _MIXED[-length:]
        g = (s.count("G") + s.count("C")) / len(s)
        if not (THRESHOLDS.pbs2_gc_lo <= g <= THRESHOLDS.pbs2_gc_hi) or _is_repetitive(s):
            continue
        t = mt.Tm_NN(s)
        if THRESHOLDS.pbs2_tm_lo <= t <= THRESHOLDS.pbs2_tm_hi:
            assert abs(tm - center) <= abs(t - center) + 1e-9


def test_select_rcpbs2_honest_miss_returns_none_with_diagnostic():
    # An all-AT window can never reach 37 C in [10,18] nt -> per-nick miss (None),
    # not a crash and not a false pass.
    rcpbs2, why = _select_rcpbs2(
        "ATATATATATATATATATAT", nick_position=20, pegrna_strand="Fw", **_band()
    )
    assert rcpbs2 is None
    assert "no rcPBS2 length" in why and "Tm" in why


def test_select_rcpbs2_reverse_strand_is_rc_of_captured_window_not_literal():
    # ADR 0014 fix: for a reverse-targeting pegRNA (TS = sense window), the emitted
    # rcPBS2 is the REVERSE COMPLEMENT of the captured sense window, not the literal
    # window (the old orient_pbs2=identity-for-Rv bug).
    rcpbs2, why = _select_rcpbs2(_MIXED, nick_position=len(_MIXED), pegrna_strand="Rv", **_band())
    assert rcpbs2 is not None, why
    captured = _MIXED[-len(rcpbs2) :]
    assert rcpbs2 == _rc(captured)  # RC of the captured target-strand window
    assert rcpbs2 != captured  # and NOT the literal window (the fixed bug)
    assert THRESHOLDS.pbs2_tm_lo <= mt.Tm_NN(rcpbs2) <= THRESHOLDS.pbs2_tm_hi


def test_select_rcpbs2_forward_strand_is_literal_sense_window():
    # For a forward-targeting pegRNA (TS = antisense), rcPBS2 = RC(RC(sense)) = the
    # sense window literally. Dual of the Rv case; both are RC(TS).
    rcpbs2, why = _select_rcpbs2(_MIXED, nick_position=len(_MIXED), pegrna_strand="Fw", **_band())
    assert rcpbs2 is not None, why
    assert rcpbs2 == _MIXED[-len(rcpbs2) :]


def test_select_rcpbs2_rejects_window_containing_n_gap():
    # An assembly-gap N-run must not yield an N-containing rcPBS2: Tm_NN silently
    # ignores N and _is_repetitive/GC do not catch it, so the alphabet is checked
    # explicitly (per-nick miss, never a non-priming N-oligo). Real on GRCh38.
    genome = "TGTTGGCCCAGTGTGAATCGNNACGTACGT"  # NN = reference gap
    rcpbs2, why = _select_rcpbs2(genome, nick_position=24, pegrna_strand="Fw", **_band())
    assert rcpbs2 is None or "N" not in rcpbs2, f"emitted non-priming N-oligo: {rcpbs2!r}"


def test_is_repetitive_flags_homopolymer_and_low_complexity():
    assert _is_repetitive("ACGAAAAAGT") is True  # 5-run of A
    assert _is_repetitive("ATATATATAT") is True  # 2 distinct bases
    assert _is_repetitive("ACGTGACGTGACA") is False  # normal
