"""Task 16: scaffold-nCas9 docking gate (Boltz-2 Protocol + fake; real adapter CX3-only)."""

import pytest

from pegrna_pipeline.hpc.complex_dock import DockFake, passes_docking_gate


def test_docking_gate_blocks_interface_degradation():
    scorer = DockFake(apo_iptm=0.85, candidate_iptm=0.40)
    ok, detail = passes_docking_gate("PEGRNA", "nCas9", scorer=scorer, max_drop=0.15)
    assert ok is False
    assert "interface" in detail


def test_docking_gate_passes_when_interface_preserved():
    scorer = DockFake(apo_iptm=0.85, candidate_iptm=0.82)
    ok, _ = passes_docking_gate("PEGRNA", "nCas9", scorer=scorer, max_drop=0.15)
    assert ok is True


def test_docking_gate_boundary_drop_blocks():
    scorer = DockFake(apo_iptm=0.85, candidate_iptm=0.70)  # drop exactly 0.15
    ok, _ = passes_docking_gate("PEGRNA", "nCas9", scorer=scorer, max_drop=0.15)
    assert ok is True  # drop == max_drop is allowed; only > blocks
    scorer2 = DockFake(apo_iptm=0.85, candidate_iptm=0.69)  # drop 0.16
    ok2, _ = passes_docking_gate("PEGRNA", "nCas9", scorer=scorer2, max_drop=0.15)
    assert ok2 is False


def test_real_boltz_adapter_errors_without_gpu():
    from pegrna_pipeline.hpc.complex_dock import Boltz2DockScorer

    with pytest.raises((ImportError, RuntimeError)):
        Boltz2DockScorer().apo_interface_iptm("nCas9")
