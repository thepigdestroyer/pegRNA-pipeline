"""Tests for :func:`build_pridict_pbs_body` (Plan 3 T6)."""

from __future__ import annotations

from pegrna_pipeline.hpc.pridict_submitter import build_pridict_pbs_body


def test_body_includes_module_load_and_conda():
    body = build_pridict_pbs_body(
        pridict_install_path="/apps/PRIDICT2",
        conda_env="pridict2",
        input_csv_path="/work/run/chunk_0.csv",
        output_dir="/work/run/out",
        cores=4,
    )
    assert "module load miniforge/3" in body
    assert "conda activate pridict2" in body
    assert "cd /apps/PRIDICT2" in body


def test_body_calls_pridict_batch_command():
    body = build_pridict_pbs_body(
        pridict_install_path="/apps/PRIDICT2",
        conda_env="pridict2",
        input_csv_path="/w/c.csv",
        output_dir="/w/o",
        cores=3,
    )
    assert "python pridict2_pegRNA_design.py batch" in body
    assert "--input-fname /w/c.csv" in body or "--input-fname c.csv" in body
    assert "--output-dir /w/o" in body
    assert "--nicking" in body
    assert "--cores 3" in body


def test_body_emits_success_marker():
    body = build_pridict_pbs_body(
        pridict_install_path="/apps/P",
        conda_env="pridict2",
        input_csv_path="/w/c.csv",
        output_dir="/w/o",
    )
    assert "PRIDICT2_BATCH_DONE" in body


def test_default_cores_is_three():
    body = build_pridict_pbs_body(
        pridict_install_path="/apps/P",
        conda_env="pridict2",
        input_csv_path="/w/c.csv",
        output_dir="/w/o",
    )
    assert "--cores 3" in body
