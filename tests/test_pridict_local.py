"""Tests for the LOCAL (subprocess) PRIDICT2 runner (:mod:`pridict_local`).

Unlike :func:`submit_pridict_batches` (CX3 SSH/PBS), this runner shells out to a
local ``python pridict2_pegRNA_design.py batch`` on the DSI box CPU. The two
PRIDICT CLI quirks it must honour (see the 2026-07-10 wiring spec, section 3):

* ``--input-fname`` is resolved relative to ``PRIDICT2/input/`` -> pass the CSV
  *basename* only, and write the CSV into ``<PRIDICT2_DIR>/input/``.
* ``--output-dir`` is relative to the repo dir and must **pre-exist** (PRIDICT
  does not ``mkdir`` it).

Every test monkeypatches the subprocess boundary; real PRIDICT2 is never run.
"""

from __future__ import annotations

import csv
import os
import subprocess
import types
from pathlib import Path

import pytest

from pegrna_pipeline.hpc import pridict_local
from pegrna_pipeline.hpc.data_types import PegRNADesign

# ---- minimal valid PRIDICT output schemas (only the required columns) --------

_PEGRNA_HEADER = [
    "sequence_name",
    "Spacer-Sequence",
    "PBSrevcomp",
    "RTrevcomp",
    "PBSlength",
    "RTlength",
    "Target-Strand",
    "Editing_Position",
    "PRIDICT2_0_editing_Score_deep_HEK",
]
_NICKING_HEADER = [
    "Nicking-Protospacer",
    "Nicking-Position-to-edit",
    "Nicking-PAMdisrupt",
    "Target_Strand",
    "DeepCas9score",
]


def _write_pegrna_csv(path: Path, name: str) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(_PEGRNA_HEADER)
        w.writerow(
            [
                f"{name}_pegRNA1",
                "GAATTCGTGTTCCTAAGCGT",
                "TCGACGCTTAGGAACAC",
                "ACGTACGTACGTACGTACGTACGTACGTAC",
                13,
                30,
                "Fw",
                50,
                0.682,
            ]
        )


def _write_nicking_csv(path: Path, name: str) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(_NICKING_HEADER)
        w.writerow([f"{name}_ng1", "-30", "Nicking_PAM_disrupt", "Fw", "0.55"])


def _make_fake_run(calls, *, returncode: int = 0, stderr: str = "", stdout: str = ""):
    """Build a ``subprocess.run`` stand-in that records the call and (on success)
    writes the two PRIDICT output CSVs the runner will parse."""

    def fake_run(cmd, *args, **kwargs):
        cwd = Path(kwargs["cwd"])
        input_fname = cmd[cmd.index("--input-fname") + 1]
        reldir = cmd[cmd.index("--output-dir") + 1]
        calls.append(
            {
                "cmd": list(cmd),
                "cwd": kwargs.get("cwd"),
                "timeout": kwargs.get("timeout"),
                "input_fname": input_fname,
                "reldir": reldir,
                # the runner must pre-create the output dir (quirk 2) and write
                # the input CSV (quirk 1) BEFORE launching the subprocess.
                "outdir_existed": (cwd / reldir).is_dir(),
                "input_existed": (cwd / "input" / input_fname).is_file(),
            }
        )
        if returncode == 0:
            with (cwd / "input" / input_fname).open(newline="") as fh:
                names = [row["sequence_name"] for row in csv.DictReader(fh)]
            outdir = cwd / reldir
            outdir.mkdir(parents=True, exist_ok=True)
            for name in names:
                _write_pegrna_csv(outdir / f"{name}_pegRNA_Pridict_full.csv", name)
                _write_nicking_csv(outdir / f"{name}_nicking_guides.csv", name)
        return subprocess.CompletedProcess(list(cmd), returncode, stdout=stdout, stderr=stderr)

    return fake_run


def _runner(tmp_path: Path, **kw):
    kw.setdefault("pridict_dir", str(tmp_path))
    kw.setdefault("pridict_python", "/fake/envs/pridict2/bin/python")
    return pridict_local.PridictLocalRunner(**kw)


# =========================================================== command / quirks


def test_run_batch_builds_command_with_basename_and_relative_outdir(tmp_path, monkeypatch):
    monkeypatch.setattr(pridict_local, "uuid4", lambda: types.SimpleNamespace(hex="abc123"))
    calls = []
    monkeypatch.setattr(pridict_local.subprocess, "run", _make_fake_run(calls))

    runner = _runner(tmp_path, cores=5, timeout=42)
    runner.run_batch([{"sequence_name": "s1", "editseq": "AAAA(+GG)TTTT"}])

    assert len(calls) == 1
    call = calls[0]
    cmd = call["cmd"]
    # python + script + batch subcommand
    assert cmd[0] == "/fake/envs/pridict2/bin/python"
    assert "pridict2_pegRNA_design.py" in cmd
    assert cmd[cmd.index("pridict2_pegRNA_design.py") + 1] == "batch"
    # QUIRK 1: --input-fname is a bare basename (no path separator), lives in input/
    assert call["input_fname"] == "abc123.csv"
    assert os.sep not in call["input_fname"]
    assert call["input_existed"] is True
    # QUIRK 2: --output-dir is repo-relative (not absolute, not under PRIDICT2_DIR
    # as an absolute path) and was pre-created before the subprocess ran.
    assert call["reldir"] == "predictions/abc123"
    assert not os.path.isabs(call["reldir"])
    assert str(tmp_path) not in call["reldir"]
    assert call["outdir_existed"] is True
    # flags + cwd + timeout
    assert "--nicking" in cmd
    assert cmd[cmd.index("--cores") + 1] == "5"
    assert call["cwd"] == str(tmp_path)
    assert call["timeout"] == 42


def test_run_batch_uses_capture_output_and_text(tmp_path, monkeypatch):
    seen = {}

    def spy_run(cmd, *args, **kwargs):
        seen.update(kwargs)
        cwd = Path(kwargs["cwd"])
        reldir = cmd[cmd.index("--output-dir") + 1]
        outdir = cwd / reldir
        outdir.mkdir(parents=True, exist_ok=True)
        _write_pegrna_csv(outdir / "s1_pegRNA_Pridict_full.csv", "s1")
        _write_nicking_csv(outdir / "s1_nicking_guides.csv", "s1")
        return subprocess.CompletedProcess(list(cmd), 0, stdout="", stderr="")

    monkeypatch.setattr(pridict_local, "uuid4", lambda: types.SimpleNamespace(hex="u"))
    monkeypatch.setattr(pridict_local.subprocess, "run", spy_run)
    _runner(tmp_path).run_batch([{"sequence_name": "s1", "editseq": "AAAA(+GG)TTTT"}])
    assert seen.get("capture_output") is True
    assert seen.get("text") is True


# =============================================================== happy path


def test_run_batch_parses_outputs_into_pegrna_designs_and_nicking_dicts(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(pridict_local.subprocess, "run", _make_fake_run(calls))

    runner = _runner(tmp_path)
    rows = [
        {"sequence_name": "s1", "editseq": "AAAA(+GG)TTTT"},
        {"sequence_name": "s2", "editseq": "CCCC(+TT)GGGG"},
    ]
    result = runner.run_batch(rows)

    assert set(result) == {"s1", "s2"}
    for name in ("s1", "s2"):
        entry = result[name]
        assert set(entry) == {"pegrnas", "nicking"}
        pegs = entry["pegrnas"]
        assert len(pegs) == 1
        assert isinstance(pegs[0], PegRNADesign)
        assert pegs[0].spacer == "GAATTCGTGTTCCTAAGCGT"
        assert pegs[0].pbs_length == 13
        assert pegs[0].target_strand == "Fw"
        nicking = entry["nicking"]
        assert isinstance(nicking, list)
        assert isinstance(nicking[0], dict)
        assert nicking[0]["Nicking-Protospacer"] == f"{name}_ng1"
        assert nicking[0]["Target_Strand"] == "Fw"


def test_run_batch_parses_from_the_repo_relative_output_dir(tmp_path, monkeypatch):
    """The pegRNA CSV path handed to the parser is
    ``<PRIDICT2_DIR>/<reldir>/<name>_pegRNA_Pridict_full.csv``."""
    monkeypatch.setattr(pridict_local, "uuid4", lambda: types.SimpleNamespace(hex="zz"))
    calls = []
    monkeypatch.setattr(pridict_local.subprocess, "run", _make_fake_run(calls))

    seen_paths = []
    sentinel = object()

    def spy_parse(path):
        seen_paths.append(str(path))
        return [sentinel]

    monkeypatch.setattr(pridict_local, "parse_pegrna_csv", spy_parse)

    runner = _runner(tmp_path)
    result = runner.run_batch([{"sequence_name": "s1", "editseq": "AAAA(+GG)TTTT"}])

    expected = tmp_path / "predictions" / "zz" / "s1_pegRNA_Pridict_full.csv"
    assert seen_paths == [str(expected)]
    assert result["s1"]["pegrnas"] == [sentinel]


def test_run_batch_writes_input_csv_with_correct_header(tmp_path, monkeypatch):
    captured = {}

    def capture_run(cmd, *args, **kwargs):
        cwd = Path(kwargs["cwd"])
        input_fname = cmd[cmd.index("--input-fname") + 1]
        with (cwd / "input" / input_fname).open(newline="") as fh:
            captured["rows"] = list(csv.reader(fh))
        reldir = cmd[cmd.index("--output-dir") + 1]
        (cwd / reldir).mkdir(parents=True, exist_ok=True)
        _write_pegrna_csv(cwd / reldir / "s1_pegRNA_Pridict_full.csv", "s1")
        _write_nicking_csv(cwd / reldir / "s1_nicking_guides.csv", "s1")
        return subprocess.CompletedProcess(list(cmd), 0, stdout="", stderr="")

    monkeypatch.setattr(pridict_local.subprocess, "run", capture_run)
    _runner(tmp_path).run_batch([{"sequence_name": "s1", "editseq": "AAAA(+GG)TTTT"}])
    assert captured["rows"][0] == ["sequence_name", "editseq"]
    assert captured["rows"][1] == ["s1", "AAAA(+GG)TTTT"]


# ================================================================ fail loud


def test_run_batch_raises_runtimeerror_on_nonzero_exit(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(
        pridict_local.subprocess,
        "run",
        _make_fake_run(calls, returncode=1, stderr="Traceback...\nRuntimeError: boom detail"),
    )
    with pytest.raises(RuntimeError, match="boom detail"):
        _runner(tmp_path).run_batch([{"sequence_name": "s1", "editseq": "AAAA(+GG)TTTT"}])


def test_run_batch_raises_runtimeerror_on_timeout(tmp_path, monkeypatch):
    def timeout_run(cmd, *args, **kwargs):
        raise subprocess.TimeoutExpired(list(cmd), kwargs.get("timeout"), stderr=b"stalled tail")

    monkeypatch.setattr(pridict_local.subprocess, "run", timeout_run)
    with pytest.raises(RuntimeError, match="stalled tail"):
        _runner(tmp_path, timeout=7).run_batch(
            [{"sequence_name": "s1", "editseq": "AAAA(+GG)TTTT"}]
        )


# ================================================================= cleanup


def test_run_batch_removes_temp_input_csv_on_success(tmp_path, monkeypatch):
    monkeypatch.setattr(pridict_local, "uuid4", lambda: types.SimpleNamespace(hex="keepme"))
    calls = []
    monkeypatch.setattr(pridict_local.subprocess, "run", _make_fake_run(calls))
    _runner(tmp_path).run_batch([{"sequence_name": "s1", "editseq": "AAAA(+GG)TTTT"}])
    assert calls[0]["input_existed"] is True  # existed during the run
    assert not (tmp_path / "input" / "keepme.csv").exists()  # removed afterwards


def test_run_batch_removes_temp_input_csv_on_failure(tmp_path, monkeypatch):
    monkeypatch.setattr(pridict_local, "uuid4", lambda: types.SimpleNamespace(hex="failcsv"))
    calls = []
    monkeypatch.setattr(
        pridict_local.subprocess, "run", _make_fake_run(calls, returncode=2, stderr="nope")
    )
    with pytest.raises(RuntimeError):
        _runner(tmp_path).run_batch([{"sequence_name": "s1", "editseq": "AAAA(+GG)TTTT"}])
    assert not (tmp_path / "input" / "failcsv.csv").exists()


# ============================================================ config / env


def test_constructs_from_env(tmp_path, monkeypatch):
    monkeypatch.setenv("PRIDICT2_DIR", str(tmp_path))
    monkeypatch.setenv("PRIDICT2_PYTHON", "/env/python")
    monkeypatch.setenv("PRIDICT2_CORES", "12")
    runner = pridict_local.PridictLocalRunner()
    assert runner.pridict_dir == str(tmp_path)
    assert runner.pridict_python == "/env/python"
    assert runner.cores == 12


def test_cores_defaults_to_8(tmp_path, monkeypatch):
    monkeypatch.delenv("PRIDICT2_CORES", raising=False)
    runner = _runner(tmp_path)
    assert runner.cores == 8


def test_missing_config_fails_loud(monkeypatch):
    monkeypatch.delenv("PRIDICT2_DIR", raising=False)
    monkeypatch.delenv("PRIDICT2_PYTHON", raising=False)
    with pytest.raises(ValueError, match="PRIDICT2_DIR"):
        pridict_local.PridictLocalRunner(pridict_python="/x/python")


# ==================================================================== fake


def test_fake_run_batch_uses_no_subprocess(monkeypatch):
    def explode(*a, **k):  # pragma: no cover - must never be called
        raise AssertionError("PridictRunnerFake must not shell out")

    monkeypatch.setattr(pridict_local.subprocess, "run", explode)
    fake = pridict_local.PridictRunnerFake()
    result = fake.run_batch([{"sequence_name": "s1", "editseq": "AAAA(+GG)TTTT"}])
    assert set(result) == {"s1"}


def test_fake_output_is_shape_compatible_with_real_runner():
    fake = pridict_local.PridictRunnerFake()
    rows = [
        {"sequence_name": "s1", "editseq": "AAAA(+GG)TTTT"},
        {"sequence_name": "s2", "editseq": "CCCC(+TT)GGGG"},
    ]
    result = fake.run_batch(rows)
    assert set(result) == {"s1", "s2"}
    for name in ("s1", "s2"):
        entry = result[name]
        assert set(entry) == {"pegrnas", "nicking"}
        assert entry["pegrnas"] and isinstance(entry["pegrnas"][0], PegRNADesign)
        assert entry["nicking"] and isinstance(entry["nicking"][0], dict)
        # nicking dicts mirror the real parser's csv.DictReader keys.
        assert "Nicking-Protospacer" in entry["nicking"][0]
        assert "Target_Strand" in entry["nicking"][0]


def test_fake_is_deterministic():
    fake = pridict_local.PridictRunnerFake()
    rows = [{"sequence_name": "s1", "editseq": "AAAA(+GG)TTTT"}]
    first = fake.run_batch(rows)
    second = fake.run_batch(rows)
    assert first["s1"]["pegrnas"][0] == second["s1"]["pegrnas"][0]
    assert first["s1"]["nicking"] == second["s1"]["nicking"]


def test_collect_missing_pegrna_csv_raises_clear_error(tmp_path):
    """PRIDICT can exit 0 yet write no pegRNA CSV (a skipped row, e.g. an insertion
    > 40 bp). _collect must fail loud with an actionable message rather than let the
    parser raise a bare FileNotFoundError."""
    runner = pridict_local.PridictLocalRunner(
        pridict_dir=str(tmp_path), pridict_python="/usr/bin/python3"
    )
    with pytest.raises(RuntimeError, match="produced no output"):
        runner._collect(tmp_path, [{"sequence_name": "peg1", "editseq": "ACGT(+GGG)ACGT"}])
