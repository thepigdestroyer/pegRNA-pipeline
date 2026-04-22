"""Verify CI workflow exists and references ruff, black, and pytest."""

from pathlib import Path

import yaml

WORKFLOW = Path(__file__).resolve().parent.parent / ".github" / "workflows" / "test.yml"


def test_workflow_file_exists():
    assert WORKFLOW.is_file()


def test_workflow_runs_ruff_black_and_pytest():
    raw = WORKFLOW.read_text()
    assert "ruff" in raw
    assert "black" in raw
    assert "pytest" in raw


def test_workflow_is_valid_yaml():
    data = yaml.safe_load(WORKFLOW.read_text())
    assert isinstance(data, dict)
    # PyYAML parses the unquoted `on:` key as Python True, so check both forms.
    triggers = data.get("on") or data.get(True)
    assert triggers, "workflow must define triggers"
    assert "jobs" in data and data["jobs"]
