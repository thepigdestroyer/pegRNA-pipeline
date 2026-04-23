"""Static checks for :file:`scripts/driver_pbs_wrapper.sh`."""

from __future__ import annotations

import os
from pathlib import Path


WRAPPER = Path(__file__).parent.parent / "scripts" / "driver_pbs_wrapper.sh"


def test_wrapper_exists():
    assert WRAPPER.exists()


def test_wrapper_is_executable():
    assert os.access(WRAPPER, os.X_OK)


def test_wrapper_shebang_is_bash():
    text = WRAPPER.read_text()
    assert text.startswith("#!/bin/bash")


def test_wrapper_loads_conda_and_activates_env():
    text = WRAPPER.read_text()
    assert "miniforge" in text
    assert "conda activate" in text
    assert "pegrna-pipeline" in text


def test_wrapper_runs_driver_module():
    text = WRAPPER.read_text()
    assert "python -m pegrna_pipeline.hpc.driver" in text


def test_wrapper_uses_run_dir_env_var():
    text = WRAPPER.read_text()
    assert "$RUN_DIR" in text


def test_wrapper_emits_done_marker():
    text = WRAPPER.read_text()
    assert "DRIVER_DONE" in text
