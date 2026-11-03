"""A design must not be able to take the machine it runs on.

Measured on the shared 64-core deployment host: one PRIDICT2 worker ran at ~36 cores
(3,607% CPU, 355 threads, 280 CPU-minutes in five minutes) and one whole-genome off-target
scan used 69 CPU-minutes. `--cores 3` only limits PRIDICT's worker PROCESSES; PyTorch and
TensorFlow still start a thread per core in each. Every public design would do this, on a
machine other people use.
"""

import sys
from pathlib import Path

import pytest

from pegrna_pipeline.hpc.pridict_local import PridictLocalRunner

pytest.importorskip("fastapi")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "web"))

import app as webapp  # noqa: E402

THREAD_VARS = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "TF_NUM_INTRAOP_THREADS",
)


def _clear(monkeypatch):
    for v in (*THREAD_VARS, "TF_NUM_INTEROP_THREADS", "PRIDICT2_THREADS", "PEGRNA_THREADS"):
        monkeypatch.delenv(v, raising=False)


# ---- the pipeline subprocess -----------------------------------------------------------


@pytest.mark.parametrize("mode", ["fake", "api"])
def test_pipeline_subprocess_thread_pools_are_capped(monkeypatch, mode):
    _clear(monkeypatch)
    monkeypatch.setattr(webapp, "THREADS", 4)
    env = webapp._run_env(mode)
    for var in THREAD_VARS:
        assert env[var] == "4", var
    assert env["TF_NUM_INTEROP_THREADS"] == "2"


def test_an_operators_explicit_thread_setting_wins(monkeypatch):
    _clear(monkeypatch)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    env = webapp._run_env("api")
    assert env["OMP_NUM_THREADS"] == "1"


def test_the_pipeline_starts_at_lower_priority(monkeypatch):
    monkeypatch.setattr(webapp, "NICE", 10)
    monkeypatch.setattr(webapp.shutil, "which", lambda b: "/usr/bin/nice" if b == "nice" else None)
    cmd = webapp._pegrna_cli("fake")
    assert cmd[:3] == ["/usr/bin/nice", "-n", "10"]
    assert cmd[-2:] == ["-m", "pegrna_pipeline.agent.cli"]


def test_nice_can_be_disabled(monkeypatch):
    monkeypatch.setattr(webapp, "NICE", 0)
    monkeypatch.setattr(webapp.shutil, "which", lambda b: "/usr/bin/nice")
    assert "nice" not in webapp._pegrna_cli("fake")[0]


def test_a_host_without_nice_still_runs(monkeypatch):
    monkeypatch.setattr(webapp, "NICE", 10)
    monkeypatch.setattr(webapp.shutil, "which", lambda b: None)
    assert webapp._pegrna_cli("fake")[1:] == ["-m", "pegrna_pipeline.agent.cli"]


# ---- the PRIDICT2 subprocess -----------------------------------------------------------


def test_pridict_workers_get_a_small_thread_budget(monkeypatch):
    _clear(monkeypatch)
    env = PridictLocalRunner._subprocess_env()
    for var in THREAD_VARS:
        assert env[var] == "2", var
    assert env["TF_NUM_INTEROP_THREADS"] == "1"


def test_pridict_threads_are_configurable(monkeypatch):
    _clear(monkeypatch)
    monkeypatch.setenv("PRIDICT2_THREADS", "6")
    env = PridictLocalRunner._subprocess_env()
    assert env["OMP_NUM_THREADS"] == "6" and env["TF_NUM_INTRAOP_THREADS"] == "6"


def test_pridict_threads_override_the_pipeline_wide_cap(monkeypatch):
    """The parent's OMP_NUM_THREADS (4, for TensorFlow/BLAS in the main process) must not
    leak into PRIDICT, whose budget is PRIDICT2_THREADS alone."""
    _clear(monkeypatch)
    monkeypatch.setenv("OMP_NUM_THREADS", "4")
    assert PridictLocalRunner._subprocess_env()["OMP_NUM_THREADS"] == "2"
