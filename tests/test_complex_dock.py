"""Task 16: scaffold-nCas9 docking gate (Boltz-2 Protocol + fake; real adapter CX3-only)."""

import pytest

from pegrna_pipeline.hpc.complex_dock import DockFake, passes_docking_gate


def test_docking_gate_blocks_below_absolute_floor():
    # candidate ipTM below the 0.6 AF-multimer floor => blocked regardless of drop
    scorer = DockFake(apo_iptm=0.85, candidate_iptm=0.40)
    ok, detail = passes_docking_gate("PEGRNA", "nCas9", scorer=scorer)
    assert ok is False
    assert "0.6" in detail or "floor" in detail


def test_docking_gate_passes_when_interface_preserved():
    scorer = DockFake(apo_iptm=0.85, candidate_iptm=0.82)
    ok, _ = passes_docking_gate("PEGRNA", "nCas9", scorer=scorer)
    assert ok is True


def test_docking_gate_blocks_on_band_drop_even_above_floor():
    # candidate 0.64 is above the 0.6 floor but dropped 0.25 (> one 0.2 band) => block
    scorer = DockFake(apo_iptm=0.89, candidate_iptm=0.64)
    ok, detail = passes_docking_gate("PEGRNA", "nCas9", scorer=scorer)
    assert ok is False
    assert "dropped" in detail


def test_docking_gate_default_thresholds_are_literature_values():
    scorer = DockFake(apo_iptm=0.85, candidate_iptm=0.70)  # drop 0.15 (<0.2), >0.6
    ok, _ = passes_docking_gate("PEGRNA", "nCas9", scorer=scorer)
    assert ok is True  # within one confidence band and above floor


def test_real_boltz_adapter_errors_without_gpu():
    from pegrna_pipeline.hpc.complex_dock import Boltz2DockScorer

    with pytest.raises((ImportError, RuntimeError)):
        Boltz2DockScorer().apo_interface_iptm("nCas9")
