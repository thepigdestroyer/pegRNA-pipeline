"""Tests for the FastAPI web layer (web/app.py) — deployment-hardening contract.

The web app wraps the ``pegrna design-recoded`` CLI as a subprocess. These tests
never launch the real pipeline: they monkeypatch the subprocess boundary and
exercise the HTTP contract (health, validation, concurrency cap, open access,
timeout handling, path-traversal guard, run retention).

Skipped automatically where FastAPI/starlette's TestClient isn't installed (so a
minimal CI without the web extras stays green); run for real once ``[dev]`` pulls
in fastapi + httpx.
"""

from __future__ import annotations

import importlib
import sys
import threading
from pathlib import Path

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")
from fastapi.testclient import TestClient  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_DIR = REPO_ROOT / "web"


@pytest.fixture()
def app_module(monkeypatch, tmp_path):
    """Import web/app.py fresh with a sandboxed RUNS dir and this repo as PEGRNA_REPO."""
    monkeypatch.setenv("PEGRNA_REPO", str(REPO_ROOT))
    monkeypatch.setenv("PEGRNA_RUNS_DIR", str(tmp_path / "runs"))
    monkeypatch.delenv("NVIDIA_API_KEY", raising=False)
    monkeypatch.setenv("NVIDIA_KEY_FILE", str(tmp_path / "no.key"))
    if str(WEB_DIR) not in sys.path:
        sys.path.insert(0, str(WEB_DIR))
    sys.modules.pop("app", None)
    mod = importlib.import_module("app")
    importlib.reload(mod)
    return mod


@pytest.fixture()
def client(app_module):
    return TestClient(app_module.app)


def _write_artifacts(run_dir: Path, *, report="REPORT OK\n", rc=0):
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "recoded_tjpe_report.txt").write_text(report)
    (run_dir / "recoded_tjpe_components.csv").write_text(
        "index,component,sequence,length\n0,spacer,ACGT,4\n"
    )
    (run_dir / "recoded_tjpe_gates.json").write_text('{"gate_results": []}')


def _fake_run_factory(**kw):
    """Return a subprocess.run stand-in that writes artifacts into --output-dir."""

    def _fake_run(cmd, **_):
        out = Path(cmd[cmd.index("--output-dir") + 1])
        _write_artifacts(out, **kw)

        class R:
            returncode = kw.get("rc", 0)
            stdout = "ok"
            stderr = ""

        return R()

    return _fake_run


# --------------------------------------------------------------------- health
def test_health_reports_deployment_state(client, app_module):
    r = client.get("/api/health")
    assert r.status_code == 200
    j = r.json()
    assert j["status"] == "ok"
    assert "fake" in j["modes"] and "api" in j["modes"]
    assert j["repo"] == str(REPO_ROOT)
    assert j["nvidia_key"] is False  # NVIDIA_KEY_FILE points at a missing file
    assert "api_python" in j and "concurrency" in j


def test_modes(client):
    j = client.get("/api/modes").json()
    assert j["modes"] == ["fake", "api"] and j["default"] == "fake"


# ----------------------------------------------------------------- validation
def test_design_rejects_unknown_mode(client):
    r = client.post("/api/design", json={"config": {}, "mode": "real"})
    assert r.status_code == 400


def test_design_rejects_oversized_config(client, app_module):
    monkey_big = {"blob": "x" * (app_module.MAX_CONFIG_BYTES + 1)}
    r = client.post("/api/design", json={"config": monkey_big, "mode": "fake"})
    assert r.status_code == 413


def test_design_rejects_overlong_genome(client, app_module):
    cfg = {"genome": "A" * (app_module.MAX_GENOME_LEN + 1)}
    r = client.post("/api/design", json={"config": cfg, "mode": "fake"})
    assert r.status_code == 400


# ------------------------------------------------------------------ fake run
def test_design_fake_success(client, app_module, monkeypatch):
    monkeypatch.setattr(app_module.subprocess, "run", _fake_run_factory(rc=0))
    r = client.post("/api/design", json={"config": {"genome": "ACGT"}, "mode": "fake"})
    assert r.status_code == 200
    j = r.json()
    assert j["ok"] is True and "REPORT OK" in j["report"] and j["has_csv"] is True


def test_design_timeout_is_structured_not_500(client, app_module, monkeypatch):
    import subprocess as _sp

    def _boom(cmd, **_):
        raise _sp.TimeoutExpired(cmd, 1)

    monkeypatch.setattr(app_module.subprocess, "run", _boom)
    r = client.post("/api/design", json={"config": {"genome": "ACGT"}, "mode": "fake"})
    assert r.status_code == 200  # graceful, not a crash
    j = r.json()
    assert j["ok"] is False and "timbefore" not in j["report"]
    assert "timed out" in (j.get("error", "") + j.get("stderr", "")).lower()


# --------------------------------------------------------------- concurrency
def test_concurrency_returns_429_when_full(client, app_module, monkeypatch):
    # Force a single slot and hold it so the request cannot acquire one.
    app_module._SLOTS = threading.BoundedSemaphore(1)
    assert app_module._SLOTS.acquire(blocking=False) is True
    try:
        r = client.post("/api/design", json={"config": {"genome": "ACGT"}, "mode": "fake"})
        assert r.status_code == 429
    finally:
        app_module._SLOTS.release()


# ---------------------------------------------------------------- open access
@pytest.mark.parametrize("endpoint", ["/api/design", "/api/design/stream"])
@pytest.mark.parametrize(
    "env_name,mode",
    [
        ("PEGRNA_ACCESS_TOKEN", "fake"),
        ("PEGRNA_ACCESS_TOKEN", "api"),
        ("PEGRNA_API_MODE_TOKEN", "api"),
    ],
)
def test_design_is_open_to_everyone_whatever_the_environment_holds(
    app_module, monkeypatch, endpoint, env_name, mode
):
    """The designer asks nobody for a token. These two variables were its old access controls;
    a deployment that still sets one must not lock anyone out."""
    monkeypatch.setenv(env_name, "secret")
    importlib.reload(app_module)  # the module reads its settings at import
    monkeypatch.setattr(app_module.subprocess, "run", _fake_run_factory(rc=0))

    def _not_launched(*_args, **_kwargs):
        raise OSError("the pipeline is not launched in tests")

    # The stream endpoint reports a failed launch inside the stream, so this stays a 200.
    monkeypatch.setattr(app_module.subprocess, "Popen", _not_launched)
    r = TestClient(app_module.app).post(endpoint, json={"config": {"genome": "ACGT"}, "mode": mode})
    assert r.status_code == 200, r.text


# ---------------------------------------------------------------- path guard
def test_run_file_rejects_unlisted_name(client):
    assert client.get("/api/runs/whatever/passwd").status_code == 404
    assert client.get("/api/runs/whatever/recoded_tjpe_report.txt").status_code == 404  # no run


# ---------------------------------------------------------------- het-edit pair
def _fake_run_with_het():
    """subprocess.run stand-in that writes the het-edit dual-pass artifacts the
    CLI produces: top-level primary artifacts + het_edit.json + het_edit/ pair."""

    def _fake_run(cmd, **_):
        out = Path(cmd[cmd.index("--output-dir") + 1])
        _write_artifacts(out)
        het_dir = out / "het_edit"
        het_dir.mkdir(parents=True, exist_ok=True)
        (het_dir / "recoded_tjpe_report.txt").write_text("HET REPORT\n")
        (het_dir / "recoded_tjpe_components.csv").write_text(
            "index,component,sequence,length\n0,spacer,TGCA,4\n"
        )
        (out / "het_edit.json").write_text(
            '{"requested": true, "produced": true, "note": "het-edit sense-window pair"}'
        )

        class R:
            returncode = 0
            stdout = "ok"
            stderr = ""

        return R()

    return _fake_run


def test_design_surfaces_het_edit_pair(client, app_module, monkeypatch):
    monkeypatch.setattr(app_module.subprocess, "run", _fake_run_with_het())
    r = client.post(
        "/api/design",
        json={"config": {"genome": "ACGT", "het_edit": True, "het_snps": [{}]}, "mode": "fake"},
    )
    assert r.status_code == 200
    het = r.json()["het_edit"]
    assert het is not None
    assert het["requested"] is True and het["produced"] is True
    assert "HET REPORT" in het["report"]
    assert het["components"] and het["components"][0]["sequence"] == "TGCA"


def test_design_single_pass_has_no_het_section(client, app_module, monkeypatch):
    monkeypatch.setattr(app_module.subprocess, "run", _fake_run_factory(rc=0))
    r = client.post("/api/design", json={"config": {"genome": "ACGT"}, "mode": "fake"})
    assert r.status_code == 200
    assert r.json()["het_edit"] is None


# ---------------------------------------------------------------- retention
def test_cleanup_keeps_last_n(app_module):
    runs = app_module.RUNS
    runs.mkdir(parents=True, exist_ok=True)
    for i in range(app_module.KEEP_RUNS + 5):
        (runs / f"run{i:04d}").mkdir()
    app_module._cleanup_runs()
    remaining = [p for p in runs.iterdir() if p.is_dir()]
    assert len(remaining) <= app_module.KEEP_RUNS
