"""FastAPI web frontend for the pegRNA (prime-editing) design pipeline.

Wraps the `pegrna design-recoded` CLI: the browser posts a design config, the
server writes it to a per-run directory, runs the pipeline, and returns the text
report + parsed CSV components. Runs the design logic on CPU. Heavy model gates
are handled by the pipeline's scorer modes (`fake` today; `api` once wired).
"""

from __future__ import annotations

import csv
import io
import json
import os
import subprocess
import sys
import time
import uuid
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

# --- paths ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
REPO = Path(os.environ.get("PEGRNA_REPO", "/mnt/data3/jmoore3/prime-editing"))
RUNS = HERE / "runs"
STATIC = HERE / "static"
EXAMPLE_CONFIG = REPO / "examples" / "ctcf_cterm" / "config.json"
RUNS.mkdir(parents=True, exist_ok=True)

# Modes offered to the browser. "fake" = design logic + stand-in scorers (light
# env). "api" = CPU-real gates (SpliceAI/MaxEntScan/ViennaRNA/AlphaFold-DB) +
# Boltz-2 docking via the NVIDIA-hosted NIM API. "real" (GPU/HPC) is not offered.
SUPPORTED_MODES = ["fake", "api"]

# Each mode runs in its own conda env: fake in the light web env, api in the
# heavy env that has TensorFlow/biotite/ViennaRNA/maxentpy installed.
CONDA = "/mnt/data3/jmoore3/miniconda3/envs"
MODE_PYTHON = {
    "fake": sys.executable,
    "api": f"{CONDA}/pegrna-real/bin/python",
}

# NVIDIA API key for the Boltz-2 NIM docking gate. Read from (in order) the
# environment or a key file, so it never has to live in the systemd unit.
NVIDIA_KEY_FILE = Path(os.environ.get("NVIDIA_KEY_FILE", "/home/jmoore3/.config/pegrna/nvidia.key"))


def _nvidia_key() -> str:
    if os.environ.get("NVIDIA_API_KEY"):
        return os.environ["NVIDIA_API_KEY"]
    if NVIDIA_KEY_FILE.exists():
        return NVIDIA_KEY_FILE.read_text().strip()
    return ""

app = FastAPI(title="pegRNA design", version="0.1.0")


class DesignRequest(BaseModel):
    config: dict[str, Any]
    mode: str = "fake"


def _pegrna_cli(mode: str) -> list[str]:
    """CLI invocation using the interpreter for the requested mode's env."""
    python = MODE_PYTHON.get(mode, sys.executable)
    return [python, "-m", "pegrna_pipeline.agent.cli"]


def _run_env(mode: str) -> dict[str, str]:
    """Environment for the pipeline subprocess (TF quiet + NVIDIA key for api)."""
    env = dict(os.environ)
    env["TF_CPP_MIN_LOG_LEVEL"] = "3"  # silence TensorFlow banner spam
    if mode == "api":
        key = _nvidia_key()
        if key:
            env["NVIDIA_API_KEY"] = key
    return env


def _parse_csv(text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    reader = csv.DictReader(io.StringIO(text))
    for row in reader:
        rows.append(dict(row))
    return rows


@app.get("/", response_class=HTMLResponse)
def index() -> HTMLResponse:
    return HTMLResponse((STATIC / "index.html").read_text())


REAL_EXAMPLE = HERE / "examples" / "ctcf_flag_real.json"


@app.get("/api/example")
def example() -> JSONResponse:
    if not EXAMPLE_CONFIG.exists():
        raise HTTPException(404, f"example config not found at {EXAMPLE_CONFIG}")
    return JSONResponse(json.loads(EXAMPLE_CONFIG.read_text()))


@app.get("/api/example/real")
def example_real() -> JSONResponse:
    """CTCF C-terminal FLAG knock-in that passes ALL real gates (incl. live
    Boltz-2 docking) in api mode."""
    if not REAL_EXAMPLE.exists():
        raise HTTPException(404, "real example config not found")
    return JSONResponse(json.loads(REAL_EXAMPLE.read_text()))


@app.get("/api/modes")
def modes() -> dict[str, Any]:
    return {"modes": SUPPORTED_MODES, "default": "fake"}


@app.get("/api/motifs")
def motifs() -> JSONResponse:
    """The motif library (epitope/affinity tags, linkers, 2A, markers, FPs, …)
    plus the T7 amplification elements, grouped by category, for the cargo builder."""
    from pegrna_pipeline.agent.motif_library import MotifLibrary
    from pegrna_pipeline.hpc.insert_elements import PROMOTERS

    groups: dict[str, list[dict[str, Any]]] = {}
    for m in MotifLibrary.load_full():
        groups.setdefault(m.category, []).append(
            {"name": m.name, "aa": m.aa, "dna": m.dna, "description": m.description}
        )
    # Frame-safe T7 promoters (from insert_elements) as a pickable category.
    for name, el in PROMOTERS.items():
        groups.setdefault("amplification", []).append(
            {"name": name, "aa": "", "dna": el.sequence,
             "description": f"{el.direction} T7 promoter for linear RNA amplification"}
        )
    for cat in groups:
        groups[cat].sort(key=lambda x: x["name"].lower())
    return JSONResponse({"groups": groups})


@app.post("/api/design")
def design(req: DesignRequest) -> JSONResponse:
    if req.mode not in SUPPORTED_MODES:
        raise HTTPException(400, f"mode must be one of {SUPPORTED_MODES}")

    run_id = time.strftime("%Y%m%d-%H%M%S") + "-" + uuid.uuid4().hex[:6]
    run_dir = RUNS / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = run_dir / "config.json"
    cfg_path.write_text(json.dumps(req.config, indent=2))

    cmd = _pegrna_cli(req.mode) + [
        "design-recoded",
        "--config",
        str(cfg_path),
        "--mode",
        req.mode,
        "--output-dir",
        str(run_dir),
    ]
    proc = subprocess.run(
        cmd,
        cwd=str(REPO),
        capture_output=True,
        text=True,
        env=_run_env(req.mode),
        timeout=int(os.environ.get("PEGRNA_TIMEOUT", "1800")),
    )
    docking_stub = req.mode == "api" and not _nvidia_key()

    report_path = run_dir / "recoded_tjpe_report.txt"
    csv_path = run_dir / "recoded_tjpe_components.csv"
    report = report_path.read_text() if report_path.exists() else ""
    components = _parse_csv(csv_path.read_text()) if csv_path.exists() else []
    gates_path = run_dir / "recoded_tjpe_gates.json"
    gates = json.loads(gates_path.read_text()) if gates_path.exists() else {}

    ok = proc.returncode == 0 and bool(report)
    return JSONResponse(
        {
            "ok": ok,
            "run_id": run_id,
            "mode": req.mode,
            "returncode": proc.returncode,
            "stdout": proc.stdout,
            "stderr": proc.stderr,
            "report": report,
            "components": components,
            "has_csv": csv_path.exists(),
            "docking_stub": docking_stub,
            "gates": gates,
        }
    )


@app.post("/api/design/stream")
def design_stream(req: DesignRequest) -> StreamingResponse:
    """Run the pipeline and stream per-gate progress as NDJSON, then a final
    `{"type":"result", ...}` line. Progress lines are `{cand,gate,info,passed}`
    emitted by the CLI's --progress hook as each gate is evaluated."""
    if req.mode not in SUPPORTED_MODES:
        raise HTTPException(400, f"mode must be one of {SUPPORTED_MODES}")
    run_id = time.strftime("%Y%m%d-%H%M%S") + "-" + uuid.uuid4().hex[:6]
    run_dir = RUNS / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "config.json").write_text(json.dumps(req.config, indent=2))
    cmd = _pegrna_cli(req.mode) + [
        "design-recoded", "--config", str(run_dir / "config.json"),
        "--mode", req.mode, "--output-dir", str(run_dir), "--progress",
    ]
    docking_stub = req.mode == "api" and not _nvidia_key()

    def gen():
        # stderr merged into stdout so a full pipe never deadlocks; non-@@P lines
        # (incl. the fail-loud NoValidPegRNA diagnostic) are collected as the log.
        proc = subprocess.Popen(
            cmd, cwd=str(REPO), env=_run_env(req.mode),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
        )
        log: list[str] = []
        for line in proc.stdout:  # type: ignore[union-attr]
            line = line.rstrip("\n")
            if line.startswith("@@P "):
                yield line[4:] + "\n"
            elif line:
                log.append(line)
        proc.wait()
        report_path = run_dir / "recoded_tjpe_report.txt"
        csv_path = run_dir / "recoded_tjpe_components.csv"
        gates_path = run_dir / "recoded_tjpe_gates.json"
        report = report_path.read_text() if report_path.exists() else ""
        components = _parse_csv(csv_path.read_text()) if csv_path.exists() else []
        gates = json.loads(gates_path.read_text()) if gates_path.exists() else {}
        final = {
            "type": "result",
            "ok": proc.returncode == 0 and bool(report),
            "run_id": run_id, "mode": req.mode, "returncode": proc.returncode,
            "stderr": "\n".join(log[-60:]), "report": report, "components": components,
            "has_csv": csv_path.exists(), "docking_stub": docking_stub, "gates": gates,
        }
        yield json.dumps(final) + "\n"

    return StreamingResponse(gen(), media_type="application/x-ndjson")


@app.get("/api/runs/{run_id}/{name}")
def run_file(run_id: str, name: str) -> FileResponse:
    # Guard against path traversal: only allow the two known artifacts.
    allowed = {"recoded_tjpe_report.txt", "recoded_tjpe_components.csv", "config.json"}
    if name not in allowed:
        raise HTTPException(404, "not found")
    path = RUNS / run_id / name
    if not path.exists():
        raise HTTPException(404, "not found")
    return FileResponse(str(path), filename=f"{run_id}_{name}")


app.mount("/static", StaticFiles(directory=str(STATIC)), name="static")
