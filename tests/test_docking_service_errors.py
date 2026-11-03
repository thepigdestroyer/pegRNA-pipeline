"""An outage must never be reported as a verdict on a design.

Boltz-2 is a rate-limited cloud service. Previously the client retried a 429 twice
(5 s, then 10 s) and then raised an ordinary exception, which the pool walk caught
per candidate as "docking scorer error" and skipped -- so a rate limit surfaced as
"no candidate survived all gates", indistinguishable from the design genuinely being
bad. Found because a calibration burst got 21 of 26 calls answered with HTTP 429.
"""

import io
import json
import urllib.error

import pytest

from pegrna_pipeline.hpc import complex_dock_nim
from pegrna_pipeline.hpc.complex_dock import DockFake, DockingServiceError
from pegrna_pipeline.hpc.complex_dock_nim import Boltz2NIMDockScorer
from pegrna_pipeline.hpc.recoded_tjpe_workflow import NoValidPegRNA, run_recoded_tjpe
from pegrna_pipeline.hpc.scorer_factory import make_scorers

GOOD_BODY = {"pair_chains_iptm_scores": [{"0": {"0": 0.0, "1": 0.7}, "1": {"0": 0.7, "1": 0.0}}]}


class _Resp:
    def __init__(self, body):
        self._b = json.dumps(body).encode()

    def read(self):
        return self._b

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _http_error(code, headers=None, body=b""):
    return urllib.error.HTTPError("http://x", code, "err", headers or {}, io.BytesIO(body))


class _Script:
    """Replaces urlopen: yields the scripted outcomes in order, counting calls."""

    def __init__(self, *outcomes):
        self.outcomes = list(outcomes)
        self.calls = 0

    def __call__(self, req, timeout=None):
        self.calls += 1
        out = self.outcomes.pop(0)
        if isinstance(out, BaseException):
            raise out
        return _Resp(out)


@pytest.fixture
def sleeps(monkeypatch):
    waited = []
    monkeypatch.setattr(complex_dock_nim.time, "sleep", lambda s: waited.append(s))
    return waited


def _scorer(**kw):
    s = Boltz2NIMDockScorer(api_key="test-key-never-sent", **kw)
    return s


def _run(scorer, script, monkeypatch):
    monkeypatch.setattr(complex_dock_nim.urllib.request, "urlopen", script)
    return scorer._predict_iptm("ACGU" * 10, "M" * 50)


def test_a_rate_limit_is_retried_and_then_succeeds(monkeypatch, sleeps):
    script = _Script(_http_error(429), _http_error(429), GOOD_BODY)
    assert _run(_scorer(), script, monkeypatch) == 0.7
    assert script.calls == 3


def test_retry_after_is_honoured(monkeypatch, sleeps):
    script = _Script(_http_error(429, {"Retry-After": "17"}), GOOD_BODY)
    _run(_scorer(), script, monkeypatch)
    assert sleeps == [17.0]


def test_backoff_grows_and_is_capped(monkeypatch, sleeps):
    script = _Script(*[_http_error(503)] * 7)
    with pytest.raises(DockingServiceError):
        _run(_scorer(retries=6), script, monkeypatch)
    assert sleeps == sorted(sleeps)
    assert max(sleeps) <= 60.0
    assert script.calls == 7


def test_exhausted_retries_raise_a_service_error_that_says_so(monkeypatch, sleeps):
    script = _Script(*[_http_error(429)] * 3)
    with pytest.raises(DockingServiceError) as exc:
        _run(_scorer(retries=2), script, monkeypatch)
    assert "not a verdict on the design" in str(exc.value)


def test_a_bad_request_is_not_retried(monkeypatch, sleeps):
    """422/401 cannot be fixed by waiting."""
    script = _Script(_http_error(422, body=b'{"detail":"bad polymer"}'))
    with pytest.raises(DockingServiceError) as exc:
        _run(_scorer(), script, monkeypatch)
    assert script.calls == 1 and "422" in str(exc.value) and "bad polymer" in str(exc.value)


def test_a_dropped_connection_is_retried(monkeypatch, sleeps):
    script = _Script(urllib.error.URLError("connection reset"), GOOD_BODY)
    assert _run(_scorer(), script, monkeypatch) == 0.7


def test_a_non_json_body_is_retried(monkeypatch, sleeps):
    class _Garbled(_Resp):
        def read(self):
            return b"<html>502 bad gateway</html>"

    class _S(_Script):
        def __call__(self, req, timeout=None):
            self.calls += 1
            out = self.outcomes.pop(0)
            return _Garbled({}) if out == "garbled" else _Resp(out)

    assert _run(_scorer(), _S("garbled", GOOD_BODY), monkeypatch) == 0.7


def test_a_response_with_no_metric_is_a_service_error(monkeypatch, sleeps):
    with pytest.raises(DockingServiceError):
        _run(_scorer(), _Script({"structures": []}), monkeypatch)


def test_a_missing_key_is_a_service_error():
    with pytest.raises(DockingServiceError):
        Boltz2NIMDockScorer(api_key=None)._predict_iptm("ACGU", "M" * 10)


# --- the part that mattered: the walk must not swallow it -----------------------------------

import json as _json  # noqa: E402
from pathlib import Path  # noqa: E402

_CFG = _json.loads(Path("examples/ctcf_cterm/config.json").read_text())


def _walk(dock):
    fake = make_scorers(mode="fake")
    return run_recoded_tjpe(
        spacer=_CFG["spacer"],
        rtt_aa=_CFG["rtt_aa"],
        cargo=_CFG["cargo"],
        genome=_CFG["genome"],
        site_s=_CFG["site_s"],
        cds_mask=_CFG["cds_mask"],
        pbs_lengths=_CFG["pbs_lengths"],
        ngrna_distance_min=_CFG["ngrna_distance_min"],
        ngrna_distance_max=_CFG["ngrna_distance_max"],
        pegrna_strand=_CFG["pegrna_strand"],
        arm_side=_CFG["arm_side"],
        max_run=_CFG["max_run"],
        splice_scorer=fake.splice_scorer,
        full_structure_scorer=fake.full_structure_scorer,
        importance_scorer=fake.importance_scorer,
        dock_scorer=dock,
        offtarget_scorer=fake.offtarget_scorer,
        maxent_scorer=fake.maxent_scorer,
    )


class _ServiceDown(DockFake):
    def interface_iptm(self, pegrna, cas9):
        raise DockingServiceError("HTTP 429 after 7 attempts")


class _BugInTheScorer(DockFake):
    def interface_iptm(self, pegrna, cas9):
        raise ValueError("unexpected shape")


def test_an_outage_aborts_the_run_instead_of_rejecting_every_candidate():
    with pytest.raises(DockingServiceError):
        _walk(_ServiceDown(apo_iptm=0.85, candidate_iptm=0.8))


def test_ordinary_scorer_faults_are_still_recorded_per_candidate():
    """Only the service-failure type escapes; other errors keep their old handling."""
    with pytest.raises(NoValidPegRNA) as exc:
        _walk(_BugInTheScorer(apo_iptm=0.85, candidate_iptm=0.8))
    assert "docking scorer error" in str(exc.value)
