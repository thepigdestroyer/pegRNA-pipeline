"""Tests for :func:`submit_with_concurrency_cap` (Plan 3 T5)."""

from __future__ import annotations

from unittest.mock import MagicMock

from pegrna_pipeline.hpc.pbs import (
    PBSJob,
    PBSJobSpec,
    submit_with_concurrency_cap,
)


def _spec():
    return PBSJobSpec(
        name="x",
        ncpus=1,
        mem_gb=1,
        walltime_hours=1,
        log_path="/x.log",
        body="true\n",
    )


def _jobs(n):
    return [PBSJob(spec=_spec(), remote_script_path=f"/w/{i}.pbs") for i in range(n)]


def test_submits_under_cap():
    ssh = MagicMock()
    ssh.exec.return_value = ("1234.pbs\n", "", 0)
    # Active count always 0 -> no backoff
    active = MagicMock(return_value=0)
    jobs = _jobs(5)
    sleeps = []
    out = submit_with_concurrency_cap(
        ssh,
        jobs,
        max_concurrent=16,
        active_count_fn=active,
        sleep_fn=sleeps.append,
    )
    assert len(out) == 5
    assert ssh.exec.call_count == 5
    assert sleeps == []


def test_backoff_when_at_cap():
    ssh = MagicMock()
    ssh.exec.return_value = ("1234.pbs\n", "", 0)
    # First two probes hit cap, then drain.
    counts = iter([16, 16, 0, 0, 0])
    active = MagicMock(side_effect=lambda _ssh: next(counts))
    jobs = _jobs(1)
    sleeps = []
    submit_with_concurrency_cap(
        ssh,
        jobs,
        max_concurrent=16,
        backoff_seconds=120,
        active_count_fn=active,
        sleep_fn=sleeps.append,
    )
    assert sleeps == [120, 120]


def test_empty_job_list_no_calls():
    ssh = MagicMock()
    active = MagicMock(return_value=0)
    out = submit_with_concurrency_cap(
        ssh,
        [],
        max_concurrent=16,
        active_count_fn=active,
        sleep_fn=lambda _s: None,
    )
    assert out == []
    ssh.exec.assert_not_called()


def test_captures_scheduler_id_per_job():
    ssh = MagicMock()
    ssh.exec.side_effect = [
        ("100.pbs\n", "", 0),
        ("101.pbs\n", "", 0),
    ]
    active = MagicMock(return_value=0)
    jobs = _jobs(2)
    submit_with_concurrency_cap(
        ssh,
        jobs,
        max_concurrent=16,
        active_count_fn=active,
        sleep_fn=lambda _s: None,
    )
    assert jobs[0].scheduler_id == "100.pbs"
    assert jobs[1].scheduler_id == "101.pbs"
