"""FastAPI web frontend for the pegRNA (prime-editing) design pipeline.

Wraps the `pegrna design-recoded` CLI (strategy "tjpe": the insert rides inside the pegRNA's
extension) and the `pegrna design-pa` CLI (strategy "pa", Prime Assembly: two pegRNAs write flaps
that a separate donor DNA matches): the browser posts a design config, the server writes it to a
per-run directory, runs the pipeline, and returns the text report + parsed CSV components. Runs the
design logic on CPU. Heavy model gates are handled by the pipeline's scorer modes (`fake` and `api`).
"""

from __future__ import annotations

import csv
import io
import json
import os
import queue
import shutil
import subprocess
import sys
import threading
import time
import uuid
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
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
# Prime Assembly example: the CTCF C-terminal cassette (499 nt, two-G T7 element) as a separate donor.
PA_EXAMPLE = REPO / "examples" / "ctcf_cterm" / "config_pa.json"
RUNS.mkdir(parents=True, exist_ok=True)

# Modes offered to the browser. "fake" = design logic + stand-in scorers (light
# env). "api" = CPU-real gates (SpliceAI/MaxEntScan/ViennaRNA/AlphaFold-DB) +
# Boltz-2 docking via the NVIDIA-hosted NIM API. "real" (GPU/HPC) is not offered.
SUPPORTED_MODES = ["fake", "api"]

# The two ways the insert can be delivered. ``command`` is the CLI subcommand; ``prefix`` starts the name of every
# file it writes into the run directory (report, components, changes, gates, and for "pa" the FASTA).
STRATEGIES = {
    "tjpe": {"command": "design-recoded", "prefix": "recoded_tjpe"},
    "pa": {"command": "design-pa", "prefix": "recoded_pa"},
}
# Prime Assembly can need a paid docking call per pegRNA, and the walk docks each pegRNA of each pair it tries, so a
# shared server bounds how much one request may ask for.
PA_MAX_PAIRS = int(os.environ.get("PEGRNA_PA_MAX_PAIRS", "60"))
PA_MAX_PBS_LENGTHS = 5

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
# A design must not be able to take the machine it runs on. The deployment host is a shared
# 64-core server: one unconstrained design ran PRIDICT2 at ~36 cores and a cas-offinder scan
# for 69 CPU-minutes. Threads are capped and the pipeline starts at a lower scheduling
# priority, so other people's jobs win any contention.
THREADS = max(1, int(os.environ.get("PEGRNA_THREADS", "4")))
NICE = max(0, int(os.environ.get("PEGRNA_NICE", "10")))
RUN_TIMEOUT = int(os.environ.get("PEGRNA_TIMEOUT", "1800"))
# How long the design stream may stay silent before emitting a keepalive. Must
# stay comfortably under the shortest idle timeout between here and the browser
# (a Cloudflare Tunnel is ~100 s).
STREAM_HEARTBEAT_SECONDS = int(os.environ.get("PEGRNA_HEARTBEAT_SECONDS", "15"))
# There is no access token: the designer is open to everyone, in both modes. What protects the
# shared box and the NVIDIA credits is the cap on simultaneous designs (429 when full), the
# thread and priority limits, and the run timeout, all above.
_SLOTS = threading.BoundedSemaphore(MAX_CONCURRENCY)


def _nvidia_key() -> str:
    if os.environ.get("NVIDIA_API_KEY"):
        return os.environ["NVIDIA_API_KEY"]
    if NVIDIA_KEY_FILE.exists():
        return NVIDIA_KEY_FILE.read_text().strip()
    return ""


def _validate_pa(config: dict[str, Any]) -> None:
    """Refuse a Prime Assembly request that cannot work (400), in the engine's own words, before a job is launched.

    It uses the engine's validator, so the page and the pipeline cannot disagree about what is valid, and adds only
    the bounds a shared server needs. A cargo shorter than two flaps is refused here with the way out (the other
    strategy), which is much kinder than a failed job.
    """
    import inspect

    from pegrna_pipeline.hpc.recoded_pa_workflow import NoValidPADesign, validate_pa_request

    missing = [k for k in ("cargo", "genome", "site_s", "cds_mask") if k not in config]
    if missing:
        raise HTTPException(400, f"config is missing {', '.join(missing)}")
    pbs = config.get("pbs_lengths")
    if isinstance(pbs, list) and len(pbs) > PA_MAX_PBS_LENGTHS:
        raise HTTPException(
            400,
            f"pbs_lengths may name at most {PA_MAX_PBS_LENGTHS} lengths (each one is folded and docked)",
        )
    pairs = config.get("max_pairs")
    if isinstance(pairs, int) and not isinstance(pairs, bool) and pairs > PA_MAX_PAIRS:
        raise HTTPException(
            400,
            f"max_pairs may be at most {PA_MAX_PAIRS} on this server (each pair can need docking)",
        )
    accepted = inspect.signature(validate_pa_request).parameters
    try:
        validate_pa_request(**{k: v for k, v in config.items() if k in accepted})
    except NoValidPADesign as exc:
        raise HTTPException(400, str(exc)) from exc
    except (ValueError, TypeError, AttributeError) as exc:
        raise HTTPException(400, f"Prime Assembly settings: {exc}") from exc


def _validate(req: "DesignRequest") -> None:
    """Reject bad strategy or mode (400), oversized config (413), or an overlong genome (400)."""
    if req.strategy not in STRATEGIES:
        raise HTTPException(400, f"strategy must be one of {list(STRATEGIES)}")
    if req.mode not in SUPPORTED_MODES:
        raise HTTPException(400, f"mode must be one of {SUPPORTED_MODES}")
    if not isinstance(req.config, dict):
        raise HTTPException(400, "config must be a JSON object")
    if len(json.dumps(req.config)) > MAX_CONFIG_BYTES:
        raise HTTPException(413, f"config exceeds {MAX_CONFIG_BYTES} bytes")
    genome = req.config.get("genome")
    if isinstance(genome, str) and len(genome) > MAX_GENOME_LEN:
        raise HTTPException(400, f"genome exceeds {MAX_GENOME_LEN} nt (api mode is CPU-bounded)")
    # The lock inputs are checked by the engine's own parser, so the page and the pipeline cannot disagree
    # about what is valid.
    cargo = req.config.get("cargo")
    if isinstance(cargo, str) and (req.config.get("cargo_parts") or req.config.get("cargo_locked")):
        from pegrna_pipeline.hpc.cargo_spec import build_cargo_spec

        try:
            build_cargo_spec(
                cargo,
                parts=req.config.get("cargo_parts"),
                locked=req.config.get("cargo_locked"),
            )
        except (ValueError, TypeError, AttributeError) as exc:
            raise HTTPException(400, f"cargo lock settings: {exc}") from exc
    if req.strategy == "pa":
        _validate_pa(req.config)


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


app = FastAPI(title="pegRNA design", version="1.0.0")

# Browser origins allowed to call this API cross-origin, comma-separated.
# The UI can be hosted away from the pipeline (e.g. the static app on Vercel,
# the models on a lab box). A real `api` design outruns a serverless function's
# execution cap, so the browser has to talk to this backend directly rather than
# being proxied -- which makes it a cross-origin call, hence CORS.
# Empty => no cross-origin access (same-origin only), which is the default.
_origins = [o.strip() for o in os.environ.get("PEGRNA_CORS_ORIGINS", "").split(",") if o.strip()]
if _origins:
    app.add_middleware(
        CORSMiddleware,
        allow_origins=_origins,
        allow_credentials=False,
        allow_methods=["GET", "POST", "OPTIONS"],
        allow_headers=["Content-Type"],
    )


@app.get("/api/health")
def health() -> dict[str, Any]:
    """Deployment health + capability probe (for systemd, monitoring, and the
    Tailscale-Funnel smoke test)."""
    api_py = MODE_PYTHON["api"]
    return {
        "status": "ok",
        "version": app.version,
        "modes": SUPPORTED_MODES,
        "strategies": list(STRATEGIES),
        "default_mode": "fake",
        "repo": str(REPO),
        "repo_ok": (REPO / "pyproject.toml").exists(),
        "runs_dir": str(RUNS),
        "api_python": api_py,
        "api_python_ok": Path(api_py).exists(),
        "nvidia_key": bool(_nvidia_key()),
        "concurrency": {"limit": MAX_CONCURRENCY},
    }


class DesignRequest(BaseModel):
    config: dict[str, Any]
    mode: str = "fake"
    #: "tjpe": the insert rides in the pegRNA's extension. "pa": Prime Assembly, a separate donor DNA.
    strategy: str = "tjpe"


def _pegrna_cli(mode: str) -> list[str]:
    """CLI invocation using the interpreter for the requested mode's env.

    Prefixed with ``nice`` when available. The ``nice`` binary rather than
    ``preexec_fn=os.nice``: this runs from a threaded server, where preexec_fn is documented
    as unsafe.
    """
    python = MODE_PYTHON.get(mode, sys.executable)
    cmd = [python, "-m", "pegrna_pipeline.agent.cli"]
    nice = shutil.which("nice")
    if NICE and nice:
        cmd = [nice, "-n", str(NICE), *cmd]
    return cmd


def _run_env(mode: str) -> dict[str, str]:
    """Environment for the pipeline subprocess (TF quiet + NVIDIA key and the
    PRIDICT2 location for api)."""
    env = dict(os.environ)
    env["TF_CPP_MIN_LOG_LEVEL"] = "3"  # silence TensorFlow banner spam
    # Cap the thread pools of everything the pipeline imports (TensorFlow for SpliceAI,
    # BLAS under numpy/scikit-learn, ViennaRNA). Left alone each starts one thread per core.
    for var in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "TF_NUM_INTRAOP_THREADS",
    ):
        env.setdefault(var, str(THREADS))
    env.setdefault("TF_NUM_INTEROP_THREADS", "2")
    if mode == "api":
        # Put the api env's bin directory on PATH. Gates that shell out to a
        # BINARY rather than import a module -- cas-offinder for off-target --
        # are installed into that env, but the subprocess otherwise inherits the
        # service's PATH and cannot see them. The gate then fails loud with
        # "requires the 'cas-offinder' binary", which reads like a missing
        # install even though it is right there.
        api_bin = str(Path(MODE_PYTHON["api"]).parent)
        env["PATH"] = api_bin + os.pathsep + env.get("PATH", "")
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


def _collect(run_dir: Path, strategy: str) -> dict[str, Any]:
    """What a finished run left in its directory, for the response: the report, the components, the gate record and
    which downloads exist. The file names depend on the strategy (``STRATEGIES``)."""
    prefix = STRATEGIES[strategy]["prefix"]
    report_path = run_dir / f"{prefix}_report.txt"
    csv_path = run_dir / f"{prefix}_components.csv"
    gates_path = run_dir / f"{prefix}_gates.json"
    gates: dict[str, Any] = {}
    if gates_path.exists():
        try:
            gates = json.loads(gates_path.read_text())
        except ValueError:
            gates = {}
    return {
        "report": report_path.read_text() if report_path.exists() else "",
        "components": _parse_csv(csv_path.read_text()) if csv_path.exists() else [],
        "has_csv": csv_path.exists(),
        "has_fasta": (run_dir / f"{prefix}_donor.fasta").exists(),
        "gates": gates,
        # only a tj-PE run can leave a het-edit pair behind; for any other run this is None
        "het_edit": _read_het_edit(run_dir),
    }


@app.get("/", response_class=HTMLResponse)
def index() -> HTMLResponse:
    return HTMLResponse((STATIC / "index.html").read_text())


# Auto-spacer minimal-FLAG preset (24 bp cargo, no manual spacer): runs the full api
# cascade through the Boltz-2 docking gate. Results are reproducible: the Boltz-2 service
# returns identical scores for identical requests, and the docking gate averages several
# diffusion samples and compares against a neutral RNA of the same length.
REAL_EXAMPLE = REPO / "examples" / "ctcf_cterm" / "config_autospacer_flag.json"


@app.get("/api/example")
def example() -> JSONResponse:
    if not EXAMPLE_CONFIG.exists():
        raise HTTPException(404, f"example config not found at {EXAMPLE_CONFIG}")
    return JSONResponse(json.loads(EXAMPLE_CONFIG.read_text()))


@app.get("/api/example/real")
def example_real() -> JSONResponse:
    """CTCF C-terminal FLAG knock-in (auto-spacer) that runs the full api cascade
    through the live Boltz-2 docking gate."""
    if not REAL_EXAMPLE.exists():
        raise HTTPException(404, "real example config not found")
    return JSONResponse(json.loads(REAL_EXAMPLE.read_text()))


@app.get("/api/example/pa")
def example_pa() -> JSONResponse:
    """The CTCF C-terminal cassette as a Prime Assembly donor: a 499 nt insert, far too long for the pegRNA's own
    extension, with the transcript the SpliceAI cross-check needs."""
    if not PA_EXAMPLE.exists():
        raise HTTPException(404, "Prime Assembly example config not found")
    return JSONResponse(json.loads(PA_EXAMPLE.read_text()))


@app.get("/api/modes")
def modes() -> dict[str, Any]:
    return {"modes": SUPPORTED_MODES, "default": "fake"}


@app.get("/api/motifs")
def motifs() -> JSONResponse:
    """The motif library (epitope/affinity tags, linkers, 2A, markers, FPs, …)
    plus the T7 amplification elements, grouped by category, for the cargo builder."""
    from pegrna_pipeline.agent.motif_library import MotifLibrary
    from pegrna_pipeline.hpc.cargo_spec import library_part_is_coding
    from pegrna_pipeline.hpc.insert_elements import PROMOTERS

    # ``locked`` is the one rule the recoder uses: a part that is not protein is never changed (ADR 0029).
    groups: dict[str, list[dict[str, Any]]] = {}
    for m in MotifLibrary.load_full():
        groups.setdefault(m.category, []).append(
            {
                "name": m.name,
                "aa": m.aa,
                "dna": m.dna,
                "description": m.description,
                "locked": not library_part_is_coding(aa=m.aa, category=m.category),
            }
        )
    # The T7 promoters (from insert_elements) as a pickable category. Never recoded.
    for name, el in PROMOTERS.items():
        groups.setdefault("amplification", []).append(
            {
                "name": name,
                "aa": "",
                "dna": el.sequence,
                "description": (
                    f"{el.direction} T7 promoter for linear RNA amplification, "
                    f"{len(el.sequence)} nt (two G's after TATA): never recoded"
                ),
                "locked": True,
            }
        )
    for cat in groups:
        groups[cat].sort(key=lambda x: x["name"].lower())
    return JSONResponse({"groups": groups})


@app.post("/api/design")
def design(req: DesignRequest) -> JSONResponse:
    _validate(req)
    with _Slot():
        run_id = time.strftime("%Y%m%d-%H%M%S") + "-" + uuid.uuid4().hex[:6]
        run_dir = RUNS / run_id
        run_dir.mkdir(parents=True, exist_ok=True)
        cfg_path = run_dir / "config.json"
        cfg_path.write_text(json.dumps(req.config, indent=2))

        cmd = _pegrna_cli(req.mode) + [
            STRATEGIES[req.strategy]["command"],
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

        found = _collect(run_dir, req.strategy)
        _cleanup_runs()

        return JSONResponse(
            {
                "ok": returncode == 0 and bool(found["report"]) and not error,
                "run_id": run_id,
                "mode": req.mode,
                "strategy": req.strategy,
                "returncode": returncode,
                "stdout": stdout,
                "stderr": stderr,
                "error": error,
                "docking_stub": docking_stub,
                **found,
            }
        )


@app.post("/api/design/stream")
def design_stream(req: DesignRequest) -> StreamingResponse:
    """Run the pipeline and stream per-gate progress as NDJSON, then a final
    `{"type":"result", ...}` line. Progress lines are `{cand,gate,info,passed}`
    emitted by the CLI's --progress hook as each gate is evaluated."""
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
        STRATEGIES[req.strategy]["command"],
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

            # Pump the subprocess on a worker thread so the generator can emit a
            # heartbeat while the pipeline is silent. A single live Boltz-2 call
            # runs for minutes with no output, and anything between the browser
            # and here -- a Cloudflare Tunnel idles at ~100 s, and most proxies
            # and load balancers are similar -- will drop a connection that sends
            # nothing. Without this a real api design dies mid-run, which looks
            # to the user like the pipeline crashed.
            lines: "queue.Queue[str | None]" = queue.Queue()

            def _pump() -> None:
                try:
                    for raw in proc.stdout:  # type: ignore[union-attr]
                        lines.put(raw.rstrip("\n"))
                finally:
                    lines.put(None)  # sentinel: stdout closed

            threading.Thread(target=_pump, daemon=True).start()

            while True:
                try:
                    line = lines.get(timeout=STREAM_HEARTBEAT_SECONDS)
                except queue.Empty:
                    # Ignored by the client (no `gate`), but it keeps the socket warm.
                    yield json.dumps({"type": "heartbeat"}) + "\n"
                    continue
                if line is None:
                    break
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
        found = _collect(run_dir, req.strategy)
        _cleanup_runs()
        final = {
            "type": "result",
            "ok": returncode == 0 and bool(found["report"]) and not error,
            "run_id": run_id,
            "mode": req.mode,
            "strategy": req.strategy,
            "returncode": returncode,
            "error": error,
            "stderr": "\n".join(log[-60:]),
            "docking_stub": docking_stub,
            **found,
        }
        yield json.dumps(final) + "\n"

    return StreamingResponse(gen(), media_type="application/x-ndjson")


@app.get("/api/runs/{run_id}/{name}")
def run_file(run_id: str, name: str) -> FileResponse:
    # Guard against path traversal: only the known artifacts, for either strategy.
    allowed = {
        "recoded_tjpe_report.txt",
        "recoded_tjpe_components.csv",
        "recoded_tjpe_changes.csv",
        "recoded_pa_report.txt",
        "recoded_pa_components.csv",
        "recoded_pa_changes.csv",
        "recoded_pa_donor.fasta",
        "config.json",
    }
    if name not in allowed:
        raise HTTPException(404, "not found")
    path = RUNS / run_id / name
    if not path.exists():
        raise HTTPException(404, "not found")
    return FileResponse(str(path), filename=f"{run_id}_{name}")


app.mount("/static", StaticFiles(directory=str(STATIC)), name="static")
