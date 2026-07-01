"""Scorer factory: select real (HPC) vs fake (CI/dry-run) gate scorers."""

import pytest

from pegrna_pipeline.hpc.scorer_factory import GateScorers, make_scorers


def test_fake_mode_returns_working_scorers():
    s = make_scorers(mode="fake")
    assert isinstance(s, GateScorers)
    # fakes are callable/usable without any heavy deps
    assert s.splice_scorer.delta("ACGT", "ACGT") == 0.0
    ok, _ = __import__(
        "pegrna_pipeline.hpc.complex_dock", fromlist=["passes_docking_gate"]
    ).passes_docking_gate("X", "Y", scorer=s.dock_scorer)
    assert ok is True
    # off-target fake marks everything clean
    assert s.offtarget_scorer(["ACGTACGTACGTACGTACGT"]) != {}


def test_fake_mode_full_structure_clean():
    s = make_scorers(mode="fake")
    assert s.full_structure_scorer.pair_probs("ACGT") == {}


def test_fake_mode_importance_all_pass():
    s = make_scorers(mode="fake")
    m = s.importance_scorer.metrics(0)
    from pegrna_pipeline.hpc.structure_importance import residue_blocked

    assert residue_blocked(m) is False


def test_real_mode_raises_without_heavy_deps():
    # On a box without torch/tensorflow/weights, real mode must fail LOUD at
    # construction/first-use — never silently downgrade to fakes.
    with pytest.raises((ImportError, RuntimeError)):
        s = make_scorers(mode="real")
        # force first-use of a real scorer if construction was lazy
        s.full_structure_scorer.pair_probs("ACGT" * 10)


def test_unknown_mode_raises():
    with pytest.raises(ValueError, match="mode"):
        make_scorers(mode="bogus")
