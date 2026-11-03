"""The api subprocess must be able to find the binaries its gates shell out to.

Most gates import a Python module, so installing them into the api environment
is enough. The off-target gate is different: it executes the `cas-offinder`
BINARY, which lives in that environment's bin directory. The pipeline subprocess
inherits the web service's PATH, which does not include it, so `shutil.which`
returns None and the gate fails loud with "requires the 'cas-offinder' binary" —
indistinguishable from it never having been installed.
"""

import os
import sys
from pathlib import Path

import pytest

pytest.importorskip("fastapi")

WEB = Path(__file__).resolve().parents[1] / "web"
sys.path.insert(0, str(WEB))

import app as webapp  # noqa: E402


def test_api_mode_puts_the_api_env_bin_on_path(monkeypatch):
    monkeypatch.setitem(webapp.MODE_PYTHON, "api", "/opt/envs/pegrna-real/bin/python")
    env = webapp._run_env("api")
    first = env["PATH"].split(os.pathsep)[0]
    assert first == "/opt/envs/pegrna-real/bin"


def test_api_mode_keeps_the_inherited_path(monkeypatch):
    monkeypatch.setitem(webapp.MODE_PYTHON, "api", "/opt/envs/pegrna-real/bin/python")
    monkeypatch.setenv("PATH", "/usr/bin:/bin")
    env = webapp._run_env("api")
    assert env["PATH"].endswith("/usr/bin:/bin")


def test_fake_mode_does_not_touch_path(monkeypatch):
    """fake mode shells out to nothing, so leave the environment alone."""
    monkeypatch.setenv("PATH", "/usr/bin:/bin")
    env = webapp._run_env("fake")
    assert env["PATH"] == "/usr/bin:/bin"
