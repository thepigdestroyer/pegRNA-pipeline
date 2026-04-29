"""Tests for :mod:`pegrna_pipeline.hpc.pbs` — render + parse + submit."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from pegrna_pipeline.hpc.pbs import (
    PBSJob,
    PBSJobSpec,
    parse_qstat_state,
    poll_until_done,
    qdel,
    qstat_one,
    qsub,
    render_pbs_script,
)


# ------------------------------------------------------------------ render
def test_render_pbs_script_basic():
    spec = PBSJobSpec(
        name="pridict_1",
        ncpus=4,
        mem_gb=16,
        walltime_hours=2,
        log_path="/work/run/pbs.log",
        body="echo hello\n",
    )
    s = render_pbs_script(spec)
    assert "#PBS -N pridict_1" in s
    assert "select=1:ncpus=4:mem=16gb" in s
    assert "walltime=02:00:00" in s
    assert "#PBS -o /work/run/pbs.log" in s
    assert "#PBS -j oe" in s
    assert s.rstrip().endswith("echo hello")


def test_render_pbs_script_with_gpu():
    spec = PBSJobSpec(
        name="gpu_job",
        ncpus=8,
        mem_gb=32,
        walltime_hours=4,
        log_path="/x/log.out",
        body="nvidia-smi\n",
        ngpus=1,
    )
    s = render_pbs_script(spec)
    assert "ncpus=8" in s and "mem=32gb" in s and "ngpus=1" in s


def test_render_pbs_script_with_queue():
    spec = PBSJobSpec(
        name="n",
        ncpus=1,
        mem_gb=1,
        walltime_hours=1,
        log_path="/x/y.log",
        body="true\n",
        queue="short",
    )
    s = render_pbs_script(spec)
    assert "#PBS -q short" in s


# ------------------------------------------------------------------ parse
def test_parse_qstat_state_running():
    out = (
        "Job id                    Name             User            Time Use S Queue\n"
        "------------------------- ---------------- --------------- -------- - -----\n"
        "1234567.pbs               pridict_1        mj              01:00:00 R workq\n"
    )
    assert parse_qstat_state(out, "1234567") == "R"


def test_parse_qstat_state_queued():
    out = (
        "Job id                    Name             User            Time Use S Queue\n"
        "------------------------- ---------------- --------------- -------- - -----\n"
        "42.pbs                    pa_driver        mj                     0 Q workq\n"
    )
    assert parse_qstat_state(out, "42") == "Q"


def test_parse_qstat_state_absent():
    assert parse_qstat_state("", "999") == ""


# ------------------------------------------------------------------ submit
def test_qsub_calls_ssh_exec_and_captures_id():
    ssh = MagicMock()
    ssh.exec.return_value = ("9999.pbs\n", "", 0)
    spec = PBSJobSpec(
        name="x",
        ncpus=1,
        mem_gb=1,
        walltime_hours=1,
        log_path="/x/log.out",
        body="true\n",
    )
    job = PBSJob(spec=spec, remote_script_path="/work/job.pbs")
    got = qsub(ssh, job)
    assert got == "9999.pbs"
    assert job.scheduler_id == "9999.pbs"
    ssh.exec.assert_called_once_with("/opt/pbs/bin/qsub /work/job.pbs")


def test_qsub_raises_on_nonzero_rc():
    ssh = MagicMock()
    ssh.exec.return_value = ("", "qsub: invalid\n", 1)
    spec = PBSJobSpec(
        name="x",
        ncpus=1,
        mem_gb=1,
        walltime_hours=1,
        log_path="/x/log.out",
        body="true\n",
    )
    job = PBSJob(spec=spec, remote_script_path="/work/bad.pbs")
    with pytest.raises(RuntimeError):
        qsub(ssh, job)


def test_qstat_one_returns_state():
    ssh = MagicMock()
    ssh.exec.return_value = ("42.pbs   name  mj  0  R  q\n", "", 0)
    spec = PBSJobSpec(
        name="n",
        ncpus=1,
        mem_gb=1,
        walltime_hours=1,
        log_path="/x",
        body="true\n",
    )
    job = PBSJob(spec=spec, remote_script_path="/a", scheduler_id="42.pbs")
    assert qstat_one(ssh, job) == "R"


def test_qdel_calls_ssh_exec():
    ssh = MagicMock()
    ssh.exec.return_value = ("", "", 0)
    spec = PBSJobSpec(
        name="n",
        ncpus=1,
        mem_gb=1,
        walltime_hours=1,
        log_path="/x",
        body="true\n",
    )
    job = PBSJob(spec=spec, remote_script_path="/a", scheduler_id="42.pbs")
    qdel(ssh, job)
    ssh.exec.assert_called_once_with("/opt/pbs/bin/qdel 42.pbs")


# ------------------------------------------------------------------ poll
def test_poll_until_done_single_job_finishes():
    ssh = MagicMock()
    spec = PBSJobSpec(
        name="n",
        ncpus=1,
        mem_gb=1,
        walltime_hours=1,
        log_path="/x",
        body="true\n",
    )
    job = PBSJob(spec=spec, remote_script_path="/a", scheduler_id="100.pbs")
    # First qstat: running; second: gone (finished).
    ssh.exec.side_effect = [
        ("100.pbs   name  mj  0  R  q\n", "", 0),
        ("", "", 0),
    ]
    sleeps = []

    def sleep_fn(secs):
        sleeps.append(secs)

    result = poll_until_done(ssh, [job], interval_seconds=1, sleep_fn=sleep_fn, max_hours=1)
    assert len(result) == 1
    assert result[0].state == "F"


def test_poll_until_done_batched_qstat():
    """The implementation must pass all job IDs in one qstat call."""
    ssh = MagicMock()
    spec = PBSJobSpec(
        name="n",
        ncpus=1,
        mem_gb=1,
        walltime_hours=1,
        log_path="/x",
        body="true\n",
    )
    jobs = [
        PBSJob(spec=spec, remote_script_path="/a", scheduler_id="100.pbs"),
        PBSJob(spec=spec, remote_script_path="/b", scheduler_id="101.pbs"),
    ]
    ssh.exec.return_value = ("", "", 0)  # both gone
    sleeps = []
    poll_until_done(ssh, jobs, interval_seconds=1, sleep_fn=sleeps.append, max_hours=1)
    cmd = ssh.exec.call_args_list[0].args[0]
    assert "100.pbs" in cmd and "101.pbs" in cmd
    assert cmd.count("qstat") == 1


def test_poll_until_done_raises_on_timeout():
    ssh = MagicMock()
    spec = PBSJobSpec(
        name="n",
        ncpus=1,
        mem_gb=1,
        walltime_hours=1,
        log_path="/x",
        body="true\n",
    )
    job = PBSJob(spec=spec, remote_script_path="/a", scheduler_id="100.pbs")
    ssh.exec.return_value = ("100.pbs name mj 0 R q", "", 0)  # never finishes
    sleeps = []
    with pytest.raises(TimeoutError):
        poll_until_done(
            ssh,
            [job],
            interval_seconds=3600,
            sleep_fn=sleeps.append,
            max_hours=0,  # immediately over budget
        )
