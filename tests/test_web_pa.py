"""The web API's Prime Assembly strategy (web/app.py, ``strategy: "pa"``).

Two kinds of test. The contract tests monkeypatch the subprocess boundary, as ``test_web_app.py`` does, to pin which
subcommand runs, what is refused before a job is launched, and what a failed run returns. The end-to-end tests drive
the REAL command line (fake mode) through the HTTP layer, so the API and the pipeline are tested together: the
progress stream, the result, the files, the downloads.
"""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")
from fastapi.testclient import TestClient  # noqa: E402
from pa_helpers import make_cargo, make_locus  # noqa: E402

S = 150
CARGO = make_cargo(120)
REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_DIR = REPO_ROOT / "web"


@pytest.fixture()
def app_module(monkeypatch, tmp_path):
    """web/app.py imported fresh, with a sandboxed runs directory and this repo as the pipeline."""
    monkeypatch.setenv("PEGRNA_REPO", str(REPO_ROOT))
    monkeypatch.setenv("PEGRNA_RUNS_DIR", str(tmp_path / "runs"))
    monkeypatch.delenv("NVIDIA_API_KEY", raising=False)
    monkeypatch.setenv("NVIDIA_KEY_FILE", str(tmp_path / "no.key"))
    if str(WEB_DIR) not in sys.path:
        sys.path.insert(0, str(WEB_DIR))
    sys.modules.pop("app", None)
    module = importlib.import_module("app")
    importlib.reload(module)
    return module


@pytest.fixture()
def client(app_module):
    return TestClient(app_module.app)


@pytest.fixture(scope="module")
def pa_config():
    loc = make_locus(
        [("top", 123, "AGG"), ("top", 147), ("bottom", 164), ("bottom", 190)],
        n=300,
        site_s=S,
        seed=7,
        coding_end=S,
    )
    return dict(cargo=CARGO, genome=loc.genome, site_s=S, cds_mask=list(loc.cds_mask), gene="TEST")


def _post(client, config, *, strategy="pa", mode="fake", path="/api/design"):
    return client.post(path, json={"config": config, "mode": mode, "strategy": strategy})


class _Spy:
    """A subprocess.run stand-in: records the command and writes the files a Prime Assembly run would."""

    def __init__(self, *, rc=0, files=("report", "components", "gates"), gates=None):
        self.cmds, self.rc, self.files, self.gates = [], rc, files, gates or {"gate_results": []}

    def __call__(self, cmd, **_):
        self.cmds.append(cmd)
        out = Path(cmd[cmd.index("--output-dir") + 1])
        out.mkdir(parents=True, exist_ok=True)
        name = {
            "report": ("recoded_pa_report.txt", "=== Prime Assembly design ===\n"),
            "components": (
                "recoded_pa_components.csv",
                "index,component,role,part_of,sequence,length\n",
            ),
            "gates": ("recoded_pa_gates.json", json.dumps(self.gates)),
            "fasta": ("recoded_pa_donor.fasta", ">donor\nACGT\n"),
        }
        for key in self.files:
            fname, text = name[key]
            (out / fname).write_text(text)

        class R:
            returncode = self.rc
            stdout = "ok"
            stderr = "boom" if self.rc else ""

        return R()


def _spy(app_module, monkeypatch, **kw):
    spy = _Spy(**kw)
    monkeypatch.setattr(app_module.subprocess, "run", spy)
    return spy


# ---- the strategy selects the command and the files ---------------------------------------------------------


def test_the_default_strategy_is_the_one_that_was_always_there(client, app_module, monkeypatch):
    seen = []

    def run(cmd, **_):
        seen.append(cmd)
        out = Path(cmd[cmd.index("--output-dir") + 1])
        out.mkdir(parents=True, exist_ok=True)
        (out / "recoded_tjpe_report.txt").write_text("REPORT")

        class R:
            returncode, stdout, stderr = 0, "", ""

        return R()

    monkeypatch.setattr(app_module.subprocess, "run", run)
    r = client.post("/api/design", json={"config": {"cargo": "ACGT"}, "mode": "fake"})
    assert r.status_code == 200 and r.json()["strategy"] == "tjpe" and r.json()["ok"] is True
    assert "design-recoded" in seen[0] and "design-pa" not in seen[0]


def test_strategy_pa_runs_design_pa_and_reads_the_pa_files(
    client, app_module, monkeypatch, pa_config
):
    spy = _spy(app_module, monkeypatch, files=("report", "components", "gates", "fasta"))
    j = _post(client, pa_config).json()
    assert "design-pa" in spy.cmds[0] and "design-recoded" not in spy.cmds[0]
    assert j["ok"] is True and j["strategy"] == "pa"
    assert j["report"].startswith("=== Prime Assembly design") and j["has_csv"] and j["has_fasta"]
    assert j["het_edit"] is None, "there is no het-edit pair in Prime Assembly"


def test_without_a_fasta_the_response_says_there_is_none(
    client, app_module, monkeypatch, pa_config
):
    _spy(app_module, monkeypatch, files=("report", "components", "gates"))
    assert _post(client, pa_config).json()["has_fasta"] is False


def test_the_tjpe_response_has_no_fasta_and_keeps_its_het_edit_key(client, app_module, monkeypatch):
    def run(cmd, **_):
        out = Path(cmd[cmd.index("--output-dir") + 1])
        out.mkdir(parents=True, exist_ok=True)
        (out / "recoded_tjpe_report.txt").write_text("REPORT")

        class R:
            returncode, stdout, stderr = 0, "", ""

        return R()

    monkeypatch.setattr(app_module.subprocess, "run", run)
    j = client.post("/api/design", json={"config": {}, "mode": "fake"}).json()
    assert j["has_fasta"] is False and "het_edit" in j


def test_a_failed_pa_run_returns_the_failure_record_the_cli_wrote(
    client, app_module, monkeypatch, pa_config
):
    record = {
        "gate_results": [],
        "failed": True,
        "error": "NoValidPADesign: no pair survived",
        "stand_in_gates": ["docking"],
        "scorer_mode": "fake",
        "strategy": "pa",
        "diagnostics": ["pair 101/164: pegRNA A: docking: ipTM too low"],
    }
    _spy(app_module, monkeypatch, rc=1, files=("gates",), gates=record)
    j = _post(client, pa_config).json()
    assert j["ok"] is False and j["returncode"] == 1 and j["report"] == "" and j["has_csv"] is False
    assert j["gates"]["failed"] is True and j["gates"]["diagnostics"] == record["diagnostics"]
    assert j["stderr"] == "boom"


def test_an_unreadable_gate_file_is_not_a_crash(client, app_module, monkeypatch, pa_config):
    spy = _spy(app_module, monkeypatch)

    def run(cmd, **kw):
        r = spy(cmd, **kw)
        out = Path(cmd[cmd.index("--output-dir") + 1])
        (out / "recoded_pa_gates.json").write_text("{not json")
        return r

    monkeypatch.setattr(app_module.subprocess, "run", run)
    j = _post(client, pa_config).json()
    assert j["gates"] == {} and j["ok"] is True


# ---- what is refused before a job is launched ----------------------------------------------------------------


def test_an_unknown_strategy_is_a_400(client, pa_config):
    r = _post(client, pa_config, strategy="twin")
    assert r.status_code == 400 and "strategy must be one of" in r.json()["detail"]


@pytest.mark.parametrize("missing", ["cargo", "genome", "site_s", "cds_mask"])
def test_a_pa_request_missing_what_the_pipeline_needs_is_a_400_naming_it(
    client, app_module, monkeypatch, pa_config, missing
):
    spy = _spy(app_module, monkeypatch)
    config = {k: v for k, v in pa_config.items() if k != missing}
    r = _post(client, config)
    assert r.status_code == 400 and r.json()["detail"] == f"config is missing {missing}"
    assert not spy.cmds


def test_a_cargo_too_short_to_assemble_is_refused_with_the_way_out_and_never_launched(
    client, app_module, monkeypatch, pa_config
):
    spy = _spy(app_module, monkeypatch)
    r = _post(client, dict(pa_config, cargo=CARGO[:60]))
    assert (
        r.status_code == 400
        and "design-recoded" in r.json()["detail"]
        and "60" in r.json()["detail"]
    )
    assert not spy.cmds, "a request that cannot work must not start a job"


def test_the_same_short_cargo_is_fine_for_the_other_strategy(
    client, app_module, monkeypatch, pa_config
):
    """The cargo rule is Prime Assembly's, not the API's: tj-PE takes a short insert."""
    spy = _spy(app_module, monkeypatch)
    r = _post(client, dict(pa_config, cargo=CARGO[:60]), strategy="tjpe")
    assert r.status_code == 200 and spy.cmds


@pytest.mark.parametrize(
    "over, words",
    [
        (dict(flap_length=29), "flap_length"),
        (dict(flap_length=51), "flap_length"),
        (dict(flap_length=35.0), "whole number"),
        (dict(flap_length="35"), "whole number"),
        (dict(pbs_lengths=[]), "pbs_lengths"),
        (dict(pbs_lengths=[13, "15"]), "pbs_lengths"),
        (dict(min_span=10), "at least 23"),
        (dict(min_span=80, max_span=60), "max_span"),
        (dict(max_offtargets=-1), "max_offtargets"),
        (dict(max_offtargets=0.5), "whole number"),
        (dict(max_pairs=0), "max_pairs"),
        (dict(tail="AAA"), "tail"),
        (dict(site_s=-4), "site_s"),
        (dict(cds_mask=[True, False]), "cds_mask"),
        (dict(cargo="ACGTN" * 20), "cargo"),
        (dict(cargo_locked=[[1, 9999]]), "cargo"),
    ],
)
def test_bad_pa_settings_are_a_400_in_the_engines_own_words(
    client, app_module, monkeypatch, pa_config, over, words
):
    spy = _spy(app_module, monkeypatch)
    r = _post(client, dict(pa_config, **over))
    assert r.status_code == 400 and words in r.json()["detail"] and not spy.cmds


def test_what_one_request_may_ask_of_a_shared_server_is_bounded(
    client, app_module, monkeypatch, pa_config
):
    spy = _spy(app_module, monkeypatch)
    many = app_module.PA_MAX_PAIRS + 1
    r = _post(client, dict(pa_config, max_pairs=many))
    assert r.status_code == 400 and f"at most {app_module.PA_MAX_PAIRS}" in r.json()["detail"]
    r = _post(client, dict(pa_config, pbs_lengths=[8, 9, 10, 11, 12, 13]))
    assert r.status_code == 400 and f"at most {app_module.PA_MAX_PBS_LENGTHS}" in r.json()["detail"]
    assert not spy.cmds
    ok = _post(
        client, dict(pa_config, max_pairs=app_module.PA_MAX_PAIRS, pbs_lengths=[13, 15, 17, 9, 11])
    )
    assert ok.status_code == 200 and spy.cmds


def test_unrecognised_keys_do_not_stop_a_pa_request(client, app_module, monkeypatch, pa_config):
    """tj-PE keys a client left in the config are provenance to the CLI, not an error."""
    spy = _spy(app_module, monkeypatch)
    r = _post(
        client, dict(pa_config, rtt_aa="L", pegrna_strand="Rv", arm_side="5prime", strategy="pa")
    )
    assert r.status_code == 200 and spy.cmds


def test_the_oversize_and_mode_checks_still_apply_to_pa(client, app_module, monkeypatch, pa_config):
    _spy(app_module, monkeypatch)
    assert _post(client, pa_config, mode="real").status_code == 400
    big = dict(pa_config, genome="A" * (app_module.MAX_GENOME_LEN + 1))
    assert _post(client, big).status_code == 400


def test_the_shipped_example_passes_the_apis_own_validation(client, app_module, monkeypatch):
    spy = _spy(app_module, monkeypatch)
    example = client.get("/api/example/pa")
    assert example.status_code == 200
    config = example.json()
    assert len(config["cargo"]) == 499 and config["transcript"] == "ENST00000264010"
    assert _post(client, config).status_code == 200 and spy.cmds


def test_the_example_is_a_404_when_the_file_is_missing(client, app_module, monkeypatch, tmp_path):
    monkeypatch.setattr(app_module, "PA_EXAMPLE", tmp_path / "nope.json")
    assert client.get("/api/example/pa").status_code == 404


def test_health_lists_the_strategies(client):
    j = client.get("/api/health").json()
    assert j["strategies"] == ["tjpe", "pa"]


# ---- downloads --------------------------------------------------------------------------------------------------


def test_the_pa_files_can_be_downloaded_and_nothing_else_can(
    client, app_module, monkeypatch, pa_config
):
    _spy(app_module, monkeypatch, files=("report", "components", "gates", "fasta"))
    run_id = _post(client, pa_config).json()["run_id"]
    run_dir = app_module.RUNS / run_id
    (run_dir / "recoded_pa_changes.csv").write_text("scope,part\n")
    for name in (
        "recoded_pa_report.txt",
        "recoded_pa_components.csv",
        "recoded_pa_changes.csv",
        "recoded_pa_donor.fasta",
        "config.json",
    ):
        r = client.get(f"/api/runs/{run_id}/{name}")
        assert r.status_code == 200, name
        assert r.headers["content-disposition"].endswith(f'{run_id}_{name}"')
    assert client.get(f"/api/runs/{run_id}/recoded_pa_donor.fasta").text == ">donor\nACGT\n"
    for name in ("recoded_pa_gates.json", "../config.json", "recoded_pa_secret.txt"):
        assert client.get(f"/api/runs/{run_id}/{name}").status_code == 404, name


def test_a_pa_download_that_does_not_exist_is_a_404(client, app_module, monkeypatch, pa_config):
    _spy(app_module, monkeypatch, files=("report", "components", "gates"))
    run_id = _post(client, pa_config).json()["run_id"]
    assert client.get(f"/api/runs/{run_id}/recoded_pa_donor.fasta").status_code == 404


# ---- end to end: the real command line behind the HTTP layer -------------------------------------------------


def test_a_real_pa_design_through_the_api_in_fake_mode(client, app_module, pa_config):
    r = _post(client, pa_config)
    j = r.json()
    assert r.status_code == 200 and j["ok"] is True, j["stderr"][-600:]
    assert j["strategy"] == "pa" and j["report"].startswith("=== Prime Assembly design")
    assert j["has_csv"] and j["has_fasta"] and j["components"]
    gates = j["gates"]
    assert gates["strategy"] == "pa" and gates["scorer_mode"] == "fake" and not gates.get("failed")
    pa = gates["pa"]
    assert [g["side"] for g in pa["guides"]] == ["A", "B"] and pa["cargo_length"] == len(CARGO)
    assert pa["donor"]["sequence"].count(CARGO) == 1
    fasta = client.get(f"/api/runs/{j['run_id']}/recoded_pa_donor.fasta").text
    assert fasta.count(">") == 3 and pa["donor"]["sequence"] in fasta.replace("\n", "")


def test_a_real_pa_design_streams_progress_and_then_the_result(client, app_module, pa_config):
    r = _post(client, pa_config, path="/api/design/stream")
    assert r.status_code == 200
    lines = [json.loads(line) for line in r.text.splitlines() if line.strip()]
    final = lines[-1]
    progress = [line for line in lines if "gate" in line]
    assert final["type"] == "result" and final["ok"] is True and final["strategy"] == "pa"
    assert {"importance", "structure_a", "structure_b", "docking_a", "docking_b"} <= {
        p["gate"] for p in progress
    }
    assert all(set(p) >= {"cand", "gate", "info", "passed", "metrics"} for p in progress)
    assert final["gates"]["pa"]["guides"] and final["has_fasta"]


def test_a_real_refusal_through_the_stream_endpoint_is_still_a_400_before_any_stream(
    client, pa_config
):
    r = _post(client, dict(pa_config, cargo=CARGO[:60]), path="/api/design/stream")
    assert r.status_code == 400 and "design-recoded" in r.json()["detail"]


def test_a_real_failed_pa_design_reports_why_in_the_result(client, app_module, pa_config):
    """No PAM-in pair fits the span: the CLI fails after validation, and the result carries its record."""
    r = _post(client, dict(pa_config, min_span=95, max_span=100), path="/api/design/stream")
    final = [json.loads(line) for line in r.text.splitlines() if line.strip()][-1]
    assert final["type"] == "result" and final["ok"] is False and final["returncode"] == 1
    assert final["gates"]["failed"] is True and final["gates"]["strategy"] == "pa"
    assert any("PAM-in" in d for d in final["gates"]["diagnostics"])
    assert "design-pa failed: NoValidPADesign" in final["stderr"]
    assert final["has_csv"] is False and final["report"] == ""
