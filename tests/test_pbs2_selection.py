"""Regression tests for the PBS2 length-selection fix.

The old design fixed PBS2 length to PBS1's lengths (13/15/17) and then rejected on
Tm. Because PBS2 is genome-determined at the nick and the 3 lengths are nested
windows there, the feasible set collapsed and whole candidate pools failed the PBS2
gate (the CTCF `no_valid_pegrna` symptom). ``_select_pbs2`` instead scans the length
range (10-18) per-nick to hit the Tm band, decoupled from PBS1.
"""

from __future__ import annotations

from Bio.SeqUtils import MeltingTemp as mt

from pegrna_pipeline.hpc.gate_thresholds import THRESHOLDS
from pegrna_pipeline.hpc.recoded_tjpe_workflow import _is_repetitive, _select_pbs2

# Window (Fw strand -> orient_pbs2 identity) where the OLD fixed lengths all miss
# the band: len13 Tm=36.5 (<37), len15 Tm=45.7 (>42), len17 Tm=51.1 (>42); but a
# length-14 PBS2 lands at Tm=40.8, in [37,42]. The nick sits at the window's end.
_MIXED = "TGTTGGCCCAGTGTGAATCG"  # 20 nt, ends at the nick


def _band():
    return dict(
        tm_lo=THRESHOLDS.pbs2_tm_lo,
        tm_hi=THRESHOLDS.pbs2_tm_hi,
        gc_lo=THRESHOLDS.pbs2_gc_lo,
        gc_hi=THRESHOLDS.pbs2_gc_hi,
        len_min=THRESHOLDS.pbs2_len_min,
        len_max=THRESHOLDS.pbs2_len_max,
    )


def test_select_pbs2_scan_succeeds_where_fixed_lengths_all_fail():
    # Sanity: the three OLD fixed lengths really do all miss the band here.
    for length, expect_in in ((13, False), (15, False), (17, False)):
        s = _MIXED[-length:]
        in_band = THRESHOLDS.pbs2_tm_lo <= mt.Tm_NN(s) <= THRESHOLDS.pbs2_tm_hi
        assert in_band is expect_in
    # The scan finds a valid PBS2 (length 14) where the fixed approach returned none.
    pbs2, why = _select_pbs2(_MIXED, nick_position=len(_MIXED), pegrna_strand="Fw", **_band())
    assert pbs2 is not None, why
    assert len(pbs2) == 14  # the length the fixed (13/15/17) set could never pick
    tm = mt.Tm_NN(pbs2)
    assert THRESHOLDS.pbs2_tm_lo <= tm <= THRESHOLDS.pbs2_tm_hi


def test_select_pbs2_returns_tm_closest_to_band_centre():
    pbs2, _why = _select_pbs2(_MIXED, nick_position=len(_MIXED), pegrna_strand="Fw", **_band())
    tm = mt.Tm_NN(pbs2)
    center = (THRESHOLDS.pbs2_tm_lo + THRESHOLDS.pbs2_tm_hi) / 2
    # No other in-band, GC/repeat-valid length is closer to the centre than the pick.
    for length in range(THRESHOLDS.pbs2_len_min, THRESHOLDS.pbs2_len_max + 1):
        s = _MIXED[-length:]
        g = (s.count("G") + s.count("C")) / len(s)
        if not (THRESHOLDS.pbs2_gc_lo <= g <= THRESHOLDS.pbs2_gc_hi) or _is_repetitive(s):
            continue
        t = mt.Tm_NN(s)
        if THRESHOLDS.pbs2_tm_lo <= t <= THRESHOLDS.pbs2_tm_hi:
            assert abs(tm - center) <= abs(t - center) + 1e-9


def test_select_pbs2_honest_miss_returns_none_with_diagnostic():
    # An all-AT window can never reach 37 C in [10,18] nt -> per-nick miss (None),
    # not a crash and not a false pass.
    pbs2, why = _select_pbs2(
        "ATATATATATATATATATAT", nick_position=20, pegrna_strand="Fw", **_band()
    )
    assert pbs2 is None
    assert "no PBS2 length" in why and "Tm" in why


def test_select_pbs2_reverse_strand_orients_then_scans():
    # For Rv the captured window is reverse-complemented (orient_pbs2); the scan
    # still returns an in-band PBS2 when one exists.
    pbs2, why = _select_pbs2(_MIXED, nick_position=len(_MIXED), pegrna_strand="Rv", **_band())
    if pbs2 is not None:  # RC changes the sequence; assert validity when found
        assert THRESHOLDS.pbs2_tm_lo <= mt.Tm_NN(pbs2) <= THRESHOLDS.pbs2_tm_hi


def test_is_repetitive_flags_homopolymer_and_low_complexity():
    assert _is_repetitive("ACGAAAAAGT") is True  # 5-run of A
    assert _is_repetitive("ATATATATAT") is True  # 2 distinct bases
    assert _is_repetitive("ACGTGACGTGACA") is False  # normal
