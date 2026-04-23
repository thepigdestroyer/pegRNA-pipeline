"""Tests for the laptop-side runs.db (SQLAlchemy + SQLite)."""

from __future__ import annotations

from pathlib import Path

import pytest

from pegrna_pipeline.agent.db import RunsDB


def test_record_and_get_run(tmp_path: Path):
    db = RunsDB(tmp_path / "runs.db")
    db.record_run(
        run_id="test_run_1",
        gene="CTCF",
        mode="pa",
        construct="GGGGS-FLAG",
        driver_job_id="8472193",
        remote_run_dir="/work/runs/test_run_1",
    )
    run = db.get("test_run_1")
    assert run.gene == "CTCF"
    assert run.mode == "pa"
    assert run.status == "submitted"
    assert run.driver_job_id == "8472193"


def test_get_missing_raises_keyerror(tmp_path: Path):
    db = RunsDB(tmp_path / "runs.db")
    with pytest.raises(KeyError):
        db.get("nosuch")


def test_list_recent_is_descending(tmp_path: Path):
    db = RunsDB(tmp_path / "runs.db")
    for i in range(3):
        db.record_run(run_id=f"run_{i}", gene="g", mode="tjpe", construct="c")
    runs = db.list_recent(limit=10)
    assert len(runs) == 3
    assert [r.run_id for r in runs] == ["run_2", "run_1", "run_0"]


def test_update_status_and_mark_cancelled(tmp_path: Path):
    db = RunsDB(tmp_path / "runs.db")
    db.record_run(run_id="r1", gene="g", mode="pa", construct="c")
    db.update_status("r1", "running")
    assert db.get("r1").status == "running"
    db.mark_cancelled("r1")
    assert db.get("r1").status == "cancelled"
    assert db.get("r1").completed_at is not None
