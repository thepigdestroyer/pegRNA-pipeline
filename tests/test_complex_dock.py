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


def test_boltz_yaml_offline_mode_adds_empty_msa_per_chain():
    # Offline single-sequence mode (default) must add `msa: empty` to EACH chain so
    # Boltz does not require the MSA server (the real GPU-smoke failure:
    # "Missing MSA's in input and --use_msa_server flag not set").
    from pegrna_pipeline.hpc.complex_dock import Boltz2DockScorer

    y = Boltz2DockScorer._build_complex_yaml("ACGU", "MKTAYIAKQR", single_sequence=True)
    assert y.count("msa: empty") == 2  # one per chain (rna + protein)
    assert "sequence: ACGU" in y and "sequence: MKTAYIAKQR" in y


def test_boltz_yaml_msa_server_mode_omits_empty_msa():
    from pegrna_pipeline.hpc.complex_dock import Boltz2DockScorer

    y = Boltz2DockScorer._build_complex_yaml("ACGU", "MKTAYIAKQR", single_sequence=False)
    assert "msa: empty" not in y


def test_boltz_defaults_to_offline_single_sequence():
    from pegrna_pipeline.hpc.complex_dock import Boltz2DockScorer

    assert Boltz2DockScorer().offline_single_sequence is True


def test_boltz_use_msa_server_disables_single_sequence():
    # Asking for the MSA server (online/login node) turns OFF single-sequence mode.
    from pegrna_pipeline.hpc.complex_dock import Boltz2DockScorer

    sc = Boltz2DockScorer(use_msa_server=True)
    assert sc.offline_single_sequence is False
    y = sc._build_complex_yaml("ACGU", "MKTAYIAKQR", single_sequence=sc.offline_single_sequence)
    assert "msa: empty" not in y


def test_boltz_cmd_defaults_to_no_kernels():
    # Default: pass --no_kernels so Boltz uses pure-PyTorch (no cuequivariance_torch
    # CUDA dependency). --cache is included when set.
    from pegrna_pipeline.hpc.complex_dock import Boltz2DockScorer

    sc = Boltz2DockScorer(cache_dir="/data/boltz_cache")
    cmd = sc._build_predict_cmd("in.yaml", "/out")
    assert "--no_kernels" in cmd
    assert "--cache" in cmd and "/data/boltz_cache" in cmd
    assert "--use_msa_server" not in cmd  # not unless requested


def test_boltz_use_kernels_true_omits_no_kernels():
    from pegrna_pipeline.hpc.complex_dock import Boltz2DockScorer

    cmd = Boltz2DockScorer(use_kernels=True)._build_predict_cmd("in.yaml", "/out")
    assert "--no_kernels" not in cmd
