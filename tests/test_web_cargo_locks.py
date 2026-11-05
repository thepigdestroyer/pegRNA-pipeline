"""The web layer exposes which parts are locked, and refuses bad lock settings with the engine's own message (ADR 0029)."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")
from fastapi.testclient import TestClient  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_DIR = REPO_ROOT / "web"
ORF = "CTGAAGGAGCTG" * 10


@pytest.fixture()
def client(monkeypatch, tmp_path):
    monkeypatch.setenv("PEGRNA_REPO", str(REPO_ROOT))
    monkeypatch.setenv("PEGRNA_RUNS_DIR", str(tmp_path / "runs"))
    monkeypatch.delenv("NVIDIA_API_KEY", raising=False)
    monkeypatch.setenv("NVIDIA_KEY_FILE", str(tmp_path / "no.key"))
    if str(WEB_DIR) not in sys.path:
        sys.path.insert(0, str(WEB_DIR))
    sys.modules.pop("app", None)
    mod = importlib.import_module("app")
    importlib.reload(mod)
    # the pipeline is never launched: a launch error is the proof the request got past validation
    monkeypatch.setattr(
        mod.subprocess, "run", lambda *a, **k: (_ for _ in ()).throw(OSError("not launched"))
    )
    return TestClient(mod.app, raise_server_exceptions=False)


def _motifs(client):
    return {
        m["name"]: m for group in client.get("/api/motifs").json()["groups"].values() for m in group
    }


def test_protein_parts_are_not_locked_and_everything_else_is(client):
    m = _motifs(client)
    assert (
        m["FLAG"]["locked"] is False
        and m["GGGGS"]["locked"] is False
        and m["P2A"]["locked"] is False
    )
    assert m["Kozak"]["locked"] is True  # a 5'UTR, not protein


def test_the_t7_elements_are_two_g_locked_and_say_so(client):
    m = _motifs(client)
    assert m["T7-G"]["dna"] == "TAATACGACTCACTATAGG" and m["rcT7-G"]["dna"] == "CCTATAGTGAGTCGTATTA"
    assert m["T7-G"]["locked"] and m["rcT7-G"]["locked"]
    assert "never recoded" in m["T7-G"]["description"] and "19 nt" in m["T7-G"]["description"]


def test_every_stop_entry_is_locked(client):
    stops = [
        m
        for g, items in client.get("/api/motifs").json()["groups"].items()
        if g == "stop"
        for m in items
    ]
    assert stops and all(m["locked"] for m in stops)


def _post(client, **cfg):
    body = {"config": {"genome": "ACGT" * 80, "cargo": ORF, **cfg}, "mode": "fake"}
    return client.post("/api/design", json=body)


@pytest.mark.parametrize(
    "cfg, text",
    [
        ({"cargo_locked": [[1, 9999]]}, "outside the cargo"),
        ({"cargo_locked": [[0, 3]]}, "outside the cargo"),
        ({"cargo_locked": [["a", 3]]}, "[first, last]"),
        ({"cargo_parts": [{"name": "x", "dna": "ACGT", "locked": False}]}, "do not add up"),
        ({"cargo_parts": [{"name": "x", "dna": "ACGN", "locked": False}]}, "ACGT"),
    ],
)
def test_bad_lock_settings_are_a_400_with_the_reason(client, cfg, text):
    r = _post(client, **cfg)
    assert r.status_code == 400 and "cargo lock settings" in r.text and text in r.text


def test_good_lock_settings_pass_validation(client):
    parts = [{"name": "ORF", "dna": ORF, "locked": False}]
    for cfg in ({"cargo_locked": [[1, 12]]}, {"cargo_parts": parts}, {}):
        assert _post(client, **cfg).status_code != 400


def test_the_changes_csv_is_served_like_the_other_run_files(client, tmp_path):
    run = tmp_path / "runs" / "r1"
    run.mkdir(parents=True)
    (run / "recoded_tjpe_changes.csv").write_text("scope,part\n")
    ok = client.get("/api/runs/r1/recoded_tjpe_changes.csv")
    assert ok.status_code == 200 and ok.text.startswith("scope,part")
    assert (
        client.get("/api/runs/r1/recoded_tjpe_secrets.csv").status_code == 404
    )  # still a whitelist
    assert (
        client.get("/api/runs/r2/recoded_tjpe_changes.csv").status_code == 404
    )  # and it must exist
