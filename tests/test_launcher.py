"""Tests for :class:`pegrna_pipeline.agent.launcher.HPCLauncher`."""

from __future__ import annotations

import json
from pathlib import Path
from typing import List

from pegrna_pipeline.agent.launcher import HPCLauncher, RunResult

CONFIG = {
    "hostname": "login.hpc",
    "username": "mj",
    "remote_runs_base": "/work/mj/pegrna-runs",
    "pridict_install_path": "/apps/PRIDICT2",
    "rscript_path": "/home/mj/pegRNA-pipeline/scripts/score_guides.R",
    "repo_dir": "/home/mj/pegRNA-pipeline",
    "conda_env": "pridict2",
    "compute_mode": "local",
    "entrez_email": "mj@test",
}


class FakeSSH:
    def __init__(self):
        self.commands: List[str] = []
        self.uploads: List[tuple] = []
        self.downloads: List[tuple] = []
        self._connected = False
        self.downloaded_files = {}

    @property
    def is_connected(self):
        return self._connected

    def connect(self):
        self._connected = True

    def disconnect(self):
        self._connected = False

    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_type, exc, tb):
        self.disconnect()

    def upload(self, local, remote):
        try:
            with open(local) as fh:
                content = fh.read()
        except Exception:
            content = ""
        self.uploads.append((local, remote, content))

    def download(self, remote, local):
        self.downloads.append((remote, local))
        # Write an empty placeholder
        Path(local).parent.mkdir(parents=True, exist_ok=True)
        Path(local).write_text(self.downloaded_files.get(remote, "stub\n"))

    def exec(self, cmd: str):
        self.commands.append(cmd)
        # Gene fetch via __main__
        if "sequence_fetcher" in cmd and "--json" in cmd:
            gene_dict = {
                "gene_name": "CTCF",
                "refseq_accession": "NM_006565.4",
                "mrna_sequence": "A" * 500,
                "cds_start": 100,
                "cds_end": 397,
                "stop_codon_start": 394,
                "stop_codon": "TAA",
                "upstream_flank": "A" * 100,
                "downstream_flank": "A" * 100,
                "flanking_bp": 100,
            }
            return (json.dumps(gene_dict) + "\n", "", 0)
        if cmd.startswith("/opt/pbs/bin/qsub "):
            return ("9999.pbs\n", "", 0)
        if cmd.startswith("/opt/pbs/bin/qstat"):
            return ("", "", 0)  # job gone -> finished
        if "DRIVER_DONE" in cmd:
            return ("1\n", "", 0)
        if cmd.startswith("mkdir") or cmd.startswith("ls"):
            return ("", "", 0)
        if cmd.startswith("find "):
            return ("/remote/run/reports/CTCF.txt\n", "", 0)
        return ("", "", 0)


def _ssh_factory(fake):
    def factory(*, hostname=None, username=None, password_provider=None, **kw):
        return fake

    return factory


def test_launch_run_happy_path(tmp_path):
    fake = FakeSSH()
    launcher = HPCLauncher(
        config=CONFIG,
        password_provider=lambda: "pw",
        ssh_factory=_ssh_factory(fake),
    )
    result = launcher.launch_run(
        gene="CTCF",
        construct="eGFP",
        mode="tjpe",
        local_output_dir=str(tmp_path / "out"),
    )
    assert isinstance(result, RunResult)
    assert result.success
    assert result.run_id
    assert result.driver_job_id == "9999.pbs"
    assert Path(result.local_reports_dir).exists()


def test_launch_run_uploads_inputs_json(tmp_path):
    fake = FakeSSH()
    launcher = HPCLauncher(
        config=CONFIG,
        password_provider=lambda: "pw",
        ssh_factory=_ssh_factory(fake),
    )
    launcher.launch_run(
        gene="CTCF",
        construct="eGFP",
        mode="tjpe",
        local_output_dir=str(tmp_path / "out"),
    )
    # inputs.json is now uploaded via ssh.upload() — check uploads, not commands.
    json_uploads = [u for u in fake.uploads if "inputs.json" in u[1]]
    assert json_uploads, "inputs.json was not uploaded"


def test_launch_run_uses_qsub_v_flag_for_run_dir(tmp_path):
    fake = FakeSSH()
    launcher = HPCLauncher(
        config=CONFIG,
        password_provider=lambda: "pw",
        ssh_factory=_ssh_factory(fake),
    )
    launcher.launch_run(
        gene="CTCF",
        construct="eGFP",
        mode="tjpe",
        local_output_dir=str(tmp_path / "out"),
    )
    qsub_cmds = [c for c in fake.commands if c.startswith("/opt/pbs/bin/qsub ")]
    assert qsub_cmds
    # Per v1.1: qsub -v RUN_DIR=...
    assert any("-v" in c and "RUN_DIR=" in c for c in qsub_cmds)


def test_poll_driver_grep_specific_log_path():
    """_poll_driver must grep <remote_run_dir>/driver.log — not a wildcard."""
    fake = FakeSSH()
    launcher = HPCLauncher(
        config=CONFIG,
        password_provider=lambda: "pw",
        ssh_factory=_ssh_factory(fake),
    )
    with fake:
        ok = launcher._poll_driver(fake, "9999.pbs", "/remote/run", interval=1)
    assert ok is True
    # grep specifically targets /remote/run/driver.log
    grep_cmds = [c for c in fake.commands if "grep" in c and "DRIVER_DONE" in c]
    assert grep_cmds
    assert any("/remote/run/driver.log" in c for c in grep_cmds)


def test_download_reports_is_public_method(tmp_path):
    """Per v1.1: download_reports (no underscore) used by 'pegrna fetch'."""
    launcher = HPCLauncher(
        config=CONFIG,
        password_provider=lambda: "pw",
        ssh_factory=_ssh_factory(FakeSSH()),
    )
    assert hasattr(launcher, "download_reports")
    # Private API should NOT exist
    assert not hasattr(launcher, "_download_reports")


def test_download_reports_pulls_files(tmp_path):
    fake = FakeSSH()
    launcher = HPCLauncher(
        config=CONFIG,
        password_provider=lambda: "pw",
        ssh_factory=_ssh_factory(fake),
    )
    dest = tmp_path / "reports"
    launcher.download_reports(
        remote_run_dir="/remote/run",
        local_reports_dir=str(dest),
    )
    assert dest.exists()
    assert len(fake.downloads) >= 1


def test_fetch_gene_record_on_login(tmp_path):
    fake = FakeSSH()
    launcher = HPCLauncher(
        config=CONFIG,
        password_provider=lambda: "pw",
        ssh_factory=_ssh_factory(fake),
    )
    gene = launcher._fetch_gene_record_on_login(fake, "CTCF")
    assert gene.gene_name == "CTCF"
    assert gene.refseq_accession == "NM_006565.4"
    # Ensure we called the login-node CLI form
    cmds = [c for c in fake.commands if "sequence_fetcher" in c]
    assert any("--json" in c for c in cmds)


def test_poll_driver_raises_when_job_disappears_without_done():
    """If qstat shows the job is gone but DRIVER_DONE was never logged, raise RuntimeError."""

    class DeadJobSSH(FakeSSH):
        def exec(self, cmd: str):
            self.commands.append(cmd)
            if cmd.startswith("/opt/pbs/bin/qstat"):
                return ("", "", 0)  # job gone, no DRIVER_DONE
            return ("0\n", "", 0)  # grep always returns 0 hits

    import pytest

    launcher = HPCLauncher(
        config=CONFIG, password_provider=lambda: "pw", ssh_factory=_ssh_factory(FakeSSH())
    )
    dead = DeadJobSSH()
    with dead:
        with pytest.raises(RuntimeError, match="DRIVER_DONE"):
            launcher._poll_driver(dead, "9999.pbs", "/remote/run", interval=1)


def test_poll_driver_raises_on_timeout():
    """TimeoutError raised when the job stays in the queue past max_hours."""

    class StuckJobSSH(FakeSSH):
        def exec(self, cmd: str):
            self.commands.append(cmd)
            if cmd.startswith("/opt/pbs/bin/qstat"):
                return ("9999.pbs  user  queue  job  R  queue\n", "", 0)
            return ("0\n", "", 0)  # grep never finds DRIVER_DONE

    import pytest

    launcher = HPCLauncher(
        config=CONFIG, password_provider=lambda: "pw", ssh_factory=_ssh_factory(FakeSSH())
    )
    stuck = StuckJobSSH()
    with stuck:
        with pytest.raises(TimeoutError):
            launcher._poll_driver(
                stuck, "9999.pbs", "/remote/run", interval=1, max_hours=0, sleep_fn=lambda _: None
            )
