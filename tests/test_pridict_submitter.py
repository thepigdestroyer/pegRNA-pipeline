"""Tests for :func:`submit_pridict_batches` (Plan 3 T7).

Uses a ``FakeSSH`` that records every command and returns canned output
in a script-driven fashion. No real SSH is contacted.
"""

from __future__ import annotations

from typing import List

from pegrna_pipeline.hpc.pridict_submitter import (
    BatchResult,
    submit_pridict_batches,
)


class FakeSSH:
    """In-memory SSH stand-in for pridict submitter tests."""

    def __init__(self, marker_hits: int = 1):
        self.commands: List[str] = []
        self.uploads: List[tuple] = []
        self._qsub_ids = iter([f"{1000 + i}.pbs" for i in range(1000)])
        self._marker_hits = marker_hits

    def upload(self, local, remote):
        self.uploads.append((local, remote))

    def exec(self, cmd: str):
        self.commands.append(cmd)
        if cmd.startswith("qsub "):
            return (next(self._qsub_ids) + "\n", "", 0)
        if cmd.startswith("qstat"):
            return ("", "", 0)  # jobs already gone -> finished
        if cmd.startswith("grep "):
            return (f"{self._marker_hits}\n", "", 0)
        return ("", "", 0)


def _variants(n):
    return [{"sequence_name": f"v_{i}", "editseq": f"ACGT[AAA/{'C' * i}]TTTT"} for i in range(n)]


def test_single_chunk_happy_path():
    ssh = FakeSSH()
    result = submit_pridict_batches(
        ssh=ssh,
        variants=_variants(5),
        remote_run_dir="/work/run",
        pridict_install_path="/apps/P",
        conda_env="pridict2",
        max_per_chunk=100,
        max_concurrent=16,
    )
    assert isinstance(result, BatchResult)
    assert len(result.submitted) == 1
    assert len(result.failed) == 0
    assert result.output_dirs == ["/work/run/out_0"]
    # Verify qsub was called once
    assert sum(1 for c in ssh.commands if c.startswith("qsub ")) == 1
    # Verify the upload went to a chunk CSV
    assert any(pair[1] == "/work/run/chunk_0.csv" for pair in ssh.uploads)


def test_multiple_chunks():
    ssh = FakeSSH()
    result = submit_pridict_batches(
        ssh=ssh,
        variants=_variants(10),
        remote_run_dir="/work/run",
        pridict_install_path="/apps/P",
        conda_env="pridict2",
        max_per_chunk=3,
        max_concurrent=16,
    )
    # 10 items / 3 per chunk -> 4 chunks (3,3,3,1)
    assert len(result.output_dirs) == 4
    assert sum(1 for c in ssh.commands if c.startswith("qsub ")) == 4


def test_missing_marker_reports_failure():
    ssh = FakeSSH(marker_hits=0)
    result = submit_pridict_batches(
        ssh=ssh,
        variants=_variants(3),
        remote_run_dir="/work/run",
        pridict_install_path="/apps/P",
        conda_env="pridict2",
        max_per_chunk=100,
        max_concurrent=16,
    )
    assert len(result.submitted) == 0
    assert len(result.failed) == 1


def test_uses_ssh_upload_for_chunk_csv():
    """Verify chunks go through sftp.upload rather than a heredoc cat."""
    ssh = FakeSSH()
    submit_pridict_batches(
        ssh=ssh,
        variants=_variants(2),
        remote_run_dir="/work/run",
        pridict_install_path="/apps/P",
        conda_env="pridict2",
        max_per_chunk=100,
        max_concurrent=16,
    )
    csv_uploads = [u for u in ssh.uploads if u[1].endswith(".csv")]
    assert len(csv_uploads) == 1
