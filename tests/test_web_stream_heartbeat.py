"""A long design must not be killed by an idle-connection timeout.

`/api/design/stream` emits per-gate progress, but a single live Boltz-2 docking
call runs for minutes with no output at all. Anything between the server and the
browser that drops idle connections -- a Cloudflare Tunnel idles at ~100 s, and
most proxies are similar -- would cut the response mid-design, which the user
cannot distinguish from a crash. The stream therefore emits a heartbeat while
the pipeline is silent.
"""

import json
import sys
from pathlib import Path

import pytest

pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

WEB = Path(__file__).resolve().parents[1] / "web"
sys.path.insert(0, str(WEB))

import app as webapp  # noqa: E402


@pytest.fixture
def client(monkeypatch, tmp_path):
    monkeypatch.setattr(webapp, "RUNS", tmp_path)
    return TestClient(webapp.app)


def _slow_pipeline(seconds: float, emit_progress: bool = False):
    """A stand-in pipeline that stays silent, like a real docking call."""
    lines = []
    if emit_progress:
        lines.append('print(r\'@@P {"cand":1,"gate":"recode","info":"x","passed":true}\')')
    code = (
        "import time,sys\n"
        + ("\n".join(lines) + "\n" if lines else "")
        + "sys.stdout.flush()\n"
        + f"time.sleep({seconds})\n"
    )
    return [sys.executable, "-c", code]


def test_heartbeat_is_emitted_while_the_pipeline_is_silent(client, monkeypatch):
    monkeypatch.setattr(webapp, "STREAM_HEARTBEAT_SECONDS", 1)
    monkeypatch.setattr(webapp, "_pegrna_cli", lambda mode: _slow_pipeline(3.5))

    with client.stream(
        "POST", "/api/design/stream", json={"config": {"genome": "ACGT"}, "mode": "fake"}
    ) as resp:
        assert resp.status_code == 200
        kinds = []
        for line in resp.iter_lines():
            if not line:
                continue
            kinds.append(json.loads(line).get("type"))

    # At least one keepalive during the silence, and exactly one terminal result.
    assert "heartbeat" in kinds, kinds
    assert kinds[-1] == "result"
    assert kinds.count("result") == 1


def test_heartbeats_do_not_replace_progress_events(client, monkeypatch):
    """A gate event still arrives as a gate event, not swallowed by keepalives."""
    monkeypatch.setattr(webapp, "STREAM_HEARTBEAT_SECONDS", 1)
    monkeypatch.setattr(webapp, "_pegrna_cli", lambda mode: _slow_pipeline(2.5, True))

    with client.stream(
        "POST", "/api/design/stream", json={"config": {"genome": "ACGT"}, "mode": "fake"}
    ) as resp:
        payloads = [json.loads(line) for line in resp.iter_lines() if line]

    gates = [p for p in payloads if p.get("gate")]
    assert len(gates) == 1
    assert gates[0]["gate"] == "recode"
    assert payloads[-1]["type"] == "result"


def test_no_heartbeat_needed_for_a_fast_run(client, monkeypatch):
    """A quick design should not be padded with keepalives."""
    monkeypatch.setattr(webapp, "STREAM_HEARTBEAT_SECONDS", 30)
    monkeypatch.setattr(webapp, "_pegrna_cli", lambda mode: _slow_pipeline(0))

    with client.stream(
        "POST", "/api/design/stream", json={"config": {"genome": "ACGT"}, "mode": "fake"}
    ) as resp:
        payloads = [json.loads(line) for line in resp.iter_lines() if line]

    assert [p.get("type") for p in payloads] == ["result"]
