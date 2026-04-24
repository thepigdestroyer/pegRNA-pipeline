"""Smoke tests for the run/status/list/cancel/fetch/setup CLI subcommands.

The `pegrna run` path itself is heavily mocked because we don't have a real
HPC connection; we just verify the command is wired and reports errors
when config is missing.
"""

from __future__ import annotations

from pathlib import Path

import yaml
from click.testing import CliRunner

from pegrna_pipeline.agent.cli import cli


def test_run_errors_when_config_missing(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("HOME", str(tmp_path))
    result = CliRunner().invoke(
        cli,
        ["run", "--gene", "CTCF", "--construct", "FLAG", "--mode", "tjpe", "--non-interactive"],
    )
    assert result.exit_code == 2
    assert "pegrna setup" in result.output


def test_status_empty_when_no_runs(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("HOME", str(tmp_path))
    result = CliRunner().invoke(cli, ["status"])
    assert result.exit_code == 0


def test_list_empty_when_no_runs(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("HOME", str(tmp_path))
    result = CliRunner().invoke(cli, ["list"])
    assert result.exit_code == 0


def test_status_specific_missing_run_errors(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("HOME", str(tmp_path))
    result = CliRunner().invoke(cli, ["status", "nosuch"])
    assert result.exit_code == 1
    assert "nosuch" in result.output


def test_setup_writes_config(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("HOME", str(tmp_path))
    answers = (
        "\n".join(
            [
                "",  # hostname (default)
                "mj1824",  # username
                "/rds/user/pegRNA-pipeline",  # repo_dir
                "",  # remote_runs_base default
                "mrinal.joshi@ukdri.ac.uk",  # entrez_email
                "/rds/user/PRIDICT2",  # pridict_install_path
                "/rds/user/pegRNA-pipeline/scripts/score_guides.R",
                "",  # conda_env default
                "local",  # compute_mode
                "",  # cell_line default
            ]
        )
        + "\n"
    )
    result = CliRunner().invoke(cli, ["setup"], input=answers)
    assert result.exit_code == 0, result.output
    cfg_path = tmp_path / ".pegrna" / "config.yaml"
    assert cfg_path.exists()
    cfg = yaml.safe_load(cfg_path.read_text())
    assert cfg["hostname"] == "login.cx3.hpc.imperial.ac.uk"
    assert cfg["username"] == "mj1824"
