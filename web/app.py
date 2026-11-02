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
import shutil
import subprocess
import sys
import threading
import time
import uuid
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

# --- paths ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent


def _default_repo() -> Path:
    """Repo root: ``PEGRNA_REPO`` if set, else the parent of ``web/`` when it
    looks like the checkout, else the historical deployment path."""
    env = os.environ.get("PEGRNA_REPO")
    if env:
        return Path(env)
    if (HERE.parent / "pyproject.toml").exists():
        return HERE.parent
    return Path("/mnt/data3/jmoore3/prime-editing")


REPO = _default_repo()
STATIC = HERE / "static"
RUNS = Path(os.environ.get("PEGRNA_RUNS_DIR", str(HERE / "runs")))
# Auto-spacer example (no manual `spacer`): matches the page's "PRIDICT designs the
# guide" flow. The 503 bp FLAG-P2A-ZeoR cassette is realistic and passes in fake mode
# (its large RTT is structure-gate-rejected in api mode -- see config_autospacer_flag
# for a minimal FLAG preset that clears the real gates).
EXAMPLE_CONFIG = REPO / "examples" / "ctcf_cterm" / "config_autospacer.json"
RUNS.mkdir(parents=True, exist_ok=True)

# Modes offered to the browser. "fake" = design logic + stand-in scorers (light
# env). "api" = CPU-real gates (SpliceAI/MaxEntScan/ViennaRNA/AlphaFold-DB) +
# Boltz-2 docking via the NVIDIA-hosted NIM API. "real" (GPU/HPC) is not offered.
SUPPORTED_MODES = ["fake", "api"]

# Each mode runs in its own conda env. The env root and both interpreter paths
# are env-overridable so this deploys on any host (not just the original box).
CONDA_ENVS = os.environ.get("PEGRNA_CONDA_ENVS", "/mnt/data3/jmoore3/miniconda3/envs")


def _api_python() -> str:
    """Interpreter for `api` mode.

    Prefers an explicit ``PEGRNA_API_PYTHON``, then the conventional
    ``pegrna-real`` conda env, and finally the interpreter already running this
    app. Without the fallback a deployment that installed the api stack into the
    current environment (pip/uv rather than that specific conda layout) would
    point at a path that does not exist and every api run would die with a bare
    "failed to launch" — the behaviour seen on any host but the original box.
    """
    explicit = os.environ.get("PEGRNA_API_PYTHON")
    if explicit:
        return explicit
    conda = Path(f"{CONDA_ENVS}/pegrna-real/bin/python")
    if conda.exists():
        return str(conda)
    return sys.executable


MODE_PYTHON = {
    "fake": os.environ.get("PEGRNA_WEB_PYTHON", sys.executable),
    "api": _api_python(),
}

# NVIDIA API key for the Boltz-2 NIM docking gate. Read from (in order) the
# environment or a key file, so it never has to live in the systemd unit.
NVIDIA_KEY_FILE = Path(
    os.environ.get("NVIDIA_KEY_FILE", str(Path.home() / ".config" / "pegrna" / "nvidia.key"))
).expanduser()

# --- deployment / abuse guards (all env-tunable) ----------------------------
# A public deployment must bound resource use: Boltz-2 NIM calls take minutes and
# cost credits, so cap concurrent runs and reject (429) rather than queue-flood.
MAX_CONCURRENCY = max(1, int(os.environ.get("PEGRNA_MAX_CONCURRENCY", "2")))
MAX_CONFIG_BYTES = int(os.environ.get("PEGRNA_MAX_CONFIG_BYTES", "200000"))  # 200 KB
MAX_GENOME_LEN = int(os.environ.get("PEGRNA_MAX_GENOME_LEN", "50000"))
KEEP_RUNS = int(os.environ.get("PEGRNA_KEEP_RUNS", "200"))
RUN_TIMEOUT = int(os.environ.get("PEGRNA_TIMEOUT", "1800"))
# Optional shared-secret gate. Empty => open (anyone can design); set it to
# require an `X-Access-Token` header (or `?token=`) on the design endpoints.
ACCESS_TOKEN = os.environ.get("PEGRNA_ACCESS_TOKEN", "")

_SLOTS = threading.BoundedSemaphore(MAX_CONCURRENCY)


def _nvidia_key() -> str:
    if os.environ.get("NVIDIA_API_KEY"):
        return os.environ["NVIDIA_API_KEY"]
    if NVIDIA_KEY_FILE.exists():
        return NVIDIA_KEY_FILE.read_text().strip()
    return ""


def _check_token(request: Request) -> None:
    """401 unless the configured access token matches (no-op when unset)."""
    if not ACCESS_TOKEN:
        return
    supplied = request.headers.get("X-Access-Token") or request.query_params.get("token") or ""
    if supplied != ACCESS_TOKEN:
        raise HTTPException(401, "missing or invalid access token")


def _validate(req: "DesignRequest") -> None:
    """Reject bad mode (400), oversized config (413), or an overlong genome (400)."""
    if req.mode not in SUPPORTED_MODES:
        raise HTTPException(400, f"mode must be one of {SUPPORTED_MODES}")
    if not isinstance(req.config, dict):
        raise HTTPException(400, "config must be a JSON object")
    if len(json.dumps(req.config)) > MAX_CONFIG_BYTES:
        raise HTTPException(413, f"config exceeds {MAX_CONFIG_BYTES} bytes")
    genome = req.config.get("genome")
    if isinstance(genome, str) and len(genome) > MAX_GENOME_LEN:
        raise HTTPException(400, f"genome exceeds {MAX_GENOME_LEN} nt (api mode is CPU-bounded)")


class _Slot:
    """Non-blocking concurrency slot; raises 429 when the pool is full. Reads the
    module-global ``_SLOTS`` at call time so it stays correct if reconfigured."""

    def __enter__(self) -> "_Slot":
        if not _SLOTS.acquire(blocking=False):
            raise HTTPException(429, f"server busy (max {MAX_CONCURRENCY} concurrent designs)")
        return self

    def __exit__(self, *exc: object) -> bool:
        _SLOTS.release()
        return False


def _cleanup_runs() -> None:
    """Keep only the newest ``KEEP_RUNS`` run directories (best-effort)."""
    try:
        dirs = sorted(
            (p for p in RUNS.iterdir() if p.is_dir()),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        for old in dirs[KEEP_RUNS:]:
            shutil.rmtree(old, ignore_errors=True)
    except FileNotFoundError:
        pass


app = FastAPI(title="pegRNA design", version="0.2.0")


@app.get("/api/health")
def health() -> dict[str, Any]:
    """Deployment health + capability probe (for systemd, monitoring, and the
    Tailscale-Funnel smoke test)."""
    api_py = MODE_PYTHON["api"]
    return {
        "status": "ok",
        "version": app.version,
        "modes": SUPPORTED_MODES,
        "default_mode": "fake",
        "repo": str(REPO),
        "repo_ok": (REPO / "pyproject.toml").exists(),
        "runs_dir": str(RUNS),
        "api_python": api_py,
        "api_python_ok": Path(api_py).exists(),
        "nvidia_key": bool(_nvidia_key()),
        "access_token_required": bool(ACCESS_TOKEN),
        "concurrency": {"limit": MAX_CONCURRENCY},
    }


class DesignRequest(BaseModel):
    config: dict[str, Any]
    mode: str = "fake"


def _pegrna_cli(mode: str) -> list[str]:
    """CLI invocation using the interpreter for the requested mode's env."""
    python = MODE_PYTHON.get(mode, sys.executable)
    return [python, "-m", "pegrna_pipeline.agent.cli"]


def _run_env(mode: str) -> dict[str, str]:
    """Environment for the pipeline subprocess (TF quiet + NVIDIA key and the
    PRIDICT2 location for api)."""
    env = dict(os.environ)
    env["TF_CPP_MIN_LOG_LEVEL"] = "3"  # silence TensorFlow banner spam
    if mode == "api":
        key = _nvidia_key()
        if key:
            env["NVIDIA_API_KEY"] = key
        # PRIDICT2 auto-spacer design (local CPU). Pass the designer's location
        # through to the pipeline when configured (mirrors NVIDIA_API_KEY above).
        for var in ("PRIDICT2_DIR", "PRIDICT2_PYTHON", "PRIDICT2_CORES"):
            val = os.environ.get(var)
            if val:
                env[var] = val
    return env


def _parse_csv(text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    reader = csv.DictReader(io.StringIO(text))
    for row in reader:
        rows.append(dict(row))
    return rows


def _read_het_edit(run_dir: Path) -> dict[str, Any] | None:
    """Het-edit dual-pass artifacts (Stage 3c), or None for a plain single pass.

    When the config requested the het-edit pair, the CLI writes ``het_edit.json``
    (requested / produced / note) at the run root and, if a het pair was designed,
    its report + components under ``het_edit/``. A single-pass run has neither, so
    this returns None and the response carries no het section."""
    meta_path = run_dir / "het_edit.json"
    if not meta_path.exists():
        return None
    try:
        meta = json.loads(meta_path.read_text())
    except (ValueError, OSError):
        return None
    het_dir = run_dir / "het_edit"
    report_path = het_dir / "recoded_tjpe_report.txt"
    csv_path = het_dir / "recoded_tjpe_components.csv"
    return {
        "requested": bool(meta.get("requested")),
        "produced": bool(meta.get("produced")),
        "note": meta.get("note", ""),
        "report": report_path.read_text() if report_path.exists() else "",
        "components": _parse_csv(csv_path.read_text()) if csv_path.exists() else [],
    }


@app.get("/", response_class=HTMLResponse)
def index() -> HTMLResponse:
    return HTMLResponse((STATIC / "index.html").read_text())


# Auto-spacer minimal-FLAG preset (24 bp cargo, no manual spacer): verified live on the
# DSI box to run the full api cascade through the Boltz-2 docking gate. Docking pass is
# stochastic (Boltz-2 is a diffusion model), so a given api run may or may not yield a
# final design; importance / structure / splice pass reliably.
REAL_EXAMPLE = REPO / "examples" / "ctcf_cterm" / "config_autospacer_flag.json"


@app.get("/api/example")
def example() -> JSONResponse:
    if not EXAMPLE_CONFIG.exists():
        raise HTTPException(404, f"example config not found at {EXAMPLE_CONFIG}")
    return JSONResponse(json.loads(EXAMPLE_CONFIG.read_text()))


@app.get("/api/example/real")
def example_real() -> JSONResponse:
    """CTCF C-terminal FLAG knock-in (auto-spacer) that runs the full api cascade
    through the live Boltz-2 docking gate. Docking pass is stochastic."""
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
            {
                "name": name,
                "aa": "",
                "dna": el.sequence,
                "description": f"{el.direction} T7 promoter for linear RNA amplification",
            }
        )
    for cat in groups:
        groups[cat].sort(key=lambda x: x["name"].lower())
    return JSONResponse({"groups": groups})


@app.post("/api/design")
def design(req: DesignRequest, request: Request) -> JSONResponse:
    _check_token(request)
    _validate(req)
    with _Slot():
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
        docking_stub = req.mode == "api" and not _nvidia_key()
        error = ""
        returncode, stdout, stderr = -1, "", ""
        try:
            proc = subprocess.run(
                cmd,
                cwd=str(REPO),
                capture_output=True,
                text=True,
                env=_run_env(req.mode),
                timeout=RUN_TIMEOUT,
            )
            returncode, stdout, stderr = proc.returncode, proc.stdout, proc.stderr
        except subprocess.TimeoutExpired:
            error = f"pipeline timed out after {RUN_TIMEOUT}s"
        except Exception as exc:  # noqa: BLE001 - surface as structured error, not a 500
            error = f"pipeline failed to launch: {exc}"

        report_path = run_dir / "recoded_tjpe_report.txt"
        csv_path = run_dir / "recoded_tjpe_components.csv"
        report = report_path.read_text() if report_path.exists() else ""
        components = _parse_csv(csv_path.read_text()) if csv_path.exists() else []
        gates_path = run_dir / "recoded_tjpe_gates.json"
        gates = json.loads(gates_path.read_text()) if gates_path.exists() else {}
        het_edit = _read_het_edit(run_dir)
        _cleanup_runs()

        return JSONResponse(
            {
                "ok": returncode == 0 and bool(report) and not error,
                "run_id": run_id,
                "mode": req.mode,
                "returncode": returncode,
                "stdout": stdout,
                "stderr": stderr,
                "error": error,
                "report": report,
                "components": components,
                "has_csv": csv_path.exists(),
                "docking_stub": docking_stub,
                "gates": gates,
                "het_edit": het_edit,
            }
        )


@app.post("/api/design/stream")
def design_stream(req: DesignRequest, request: Request) -> StreamingResponse:
    """Run the pipeline and stream per-gate progress as NDJSON, then a final
    `{"type":"result", ...}` line. Progress lines are `{cand,gate,info,passed}`
    emitted by the CLI's --progress hook as each gate is evaluated."""
    _check_token(request)
    _validate(req)
    # Acquire a slot up front so an over-capacity request is rejected (429)
    # before any work starts; the generator releases it in a finally.
    if not _SLOTS.acquire(blocking=False):
        raise HTTPException(429, f"server busy (max {MAX_CONCURRENCY} concurrent designs)")
    run_id = time.strftime("%Y%m%d-%H%M%S") + "-" + uuid.uuid4().hex[:6]
    run_dir = RUNS / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "config.json").write_text(json.dumps(req.config, indent=2))
    cmd = _pegrna_cli(req.mode) + [
        "design-recoded",
        "--config",
        str(run_dir / "config.json"),
        "--mode",
        req.mode,
        "--output-dir",
        str(run_dir),
        "--progress",
    ]
    docking_stub = req.mode == "api" and not _nvidia_key()

    def gen():
        error = ""
        returncode = -1
        log: list[str] = []
        try:
            # stderr merged into stdout so a full pipe never deadlocks; non-@@P lines
            # (incl. the fail-loud NoValidPegRNA diagnostic) are collected as the log.
            proc = subprocess.Popen(
                cmd,
                cwd=str(REPO),
                env=_run_env(req.mode),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            for line in proc.stdout:  # type: ignore[union-attr]
                line = line.rstrip("\n")
                if line.startswith("@@P "):
                    yield line[4:] + "\n"
                elif line:
                    log.append(line)
            proc.wait()
            returncode = proc.returncode
        except Exception as exc:  # noqa: BLE001 - surface as a result line, not a broken stream
            error = f"pipeline failed to launch: {exc}"
        finally:
            _SLOTS.release()
        report_path = run_dir / "recoded_tjpe_report.txt"
        csv_path = run_dir / "recoded_tjpe_components.csv"
        gates_path = run_dir / "recoded_tjpe_gates.json"
        report = report_path.read_text() if report_path.exists() else ""
        components = _parse_csv(csv_path.read_text()) if csv_path.exists() else []
        gates = json.loads(gates_path.read_text()) if gates_path.exists() else {}
        het_edit = _read_het_edit(run_dir)
        _cleanup_runs()
        final = {
            "type": "result",
            "ok": returncode == 0 and bool(report) and not error,
            "run_id": run_id,
            "mode": req.mode,
            "returncode": returncode,
            "error": error,
            "stderr": "\n".join(log[-60:]),
            "report": report,
            "components": components,
            "has_csv": csv_path.exists(),
            "docking_stub": docking_stub,
            "gates": gates,
            "het_edit": het_edit,
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
