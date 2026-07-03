"""rcPBS2 capture validation in the TPRT source-of-truth simulator (ADR 0014).

The tj-pegRNA carries ``rcPBS2 = RC(TS)``. The RT copies it into the newly
synthesised flap, where it appears as its reverse complement ``PBS2 = RC(rcPBS2)``.
That flap-borne PBS2 must equal the target-strand window it captures at the ngRNA
nick. A non-RC'd rcPBS2 (the old ``orient_pbs2``=TS-literal-for-Rv bug) makes
``RC(rcpbs2) != TS window`` and must fail loudly.
"""

from __future__ import annotations

import pytest
from Bio.Seq import Seq

from pegrna_pipeline.hpc.strand_orientation import orient_rcpbs2
from pegrna_pipeline.hpc.tprt_sim import TPRTMismatch, simulate_rcpbs2_capture


def _rc(s: str) -> str:
    return str(Seq(s).reverse_complement())


def test_correct_rcpbs2_reconstructs_the_capture():
    ts = "GCATGACTGA"
    rcpbs2 = orient_rcpbs2(ts)  # = RC(TS) = "TCAGTCATGC"
    pbs2_on_flap = simulate_rcpbs2_capture(rcpbs2=rcpbs2, ts_capture_window=ts)
    # PBS2 that ends up on the flap is RC(rcPBS2), which recovers the TS window.
    assert pbs2_on_flap == _rc(rcpbs2) == ts


def test_non_rcd_rcpbs2_raises_tprt_mismatch():
    # The old bug: the TS literal supplied where RC(TS) belongs.
    ts = "GCATGACTGA"
    with pytest.raises(TPRTMismatch):
        simulate_rcpbs2_capture(rcpbs2=ts, ts_capture_window=ts)


def test_expected_pbs2_mismatch_raises():
    ts = "GCATGACTGA"
    rcpbs2 = orient_rcpbs2(ts)
    with pytest.raises(TPRTMismatch):
        simulate_rcpbs2_capture(rcpbs2=rcpbs2, ts_capture_window=ts, expected_pbs2="AAAAAAAAAA")


def test_returns_pbs2_when_expected_matches():
    ts = "GCATGACTGA"
    rcpbs2 = orient_rcpbs2(ts)
    out = simulate_rcpbs2_capture(rcpbs2=rcpbs2, ts_capture_window=ts, expected_pbs2=ts)
    assert out == ts


def test_orientation_helper_and_simulator_agree_round_trip():
    # Regression fortification: orient_rcpbs2 and the capture simulator must agree for
    # any target-strand window (guards a future re-introduction of the orientation bug).
    for ts in ("ACGTTGCA", "TTTGGGCCCAAT", "GATTACAGATTACA", "CCCCGGGGAAAA"):
        rcpbs2 = orient_rcpbs2(ts)
        assert simulate_rcpbs2_capture(rcpbs2=rcpbs2, ts_capture_window=ts) == ts
