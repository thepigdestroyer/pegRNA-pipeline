#!/usr/bin/env bash
# Build the two conda envs the web app needs, on a GPU-less host (e.g. the DSI box).
#
#   pegrna-web   light: the CLI package + FastAPI/uvicorn. Serves `fake` mode.
#   pegrna-real  heavy: adds the CPU/API model stack for `api` mode
#                (ViennaRNA, SpliceAI+TensorFlow-CPU, biotite, dnachisel, maxentpy).
#                NO torch/boltz — docking is the NVIDIA NIM API, not local GPU.
#
# Idempotent and fail-loud. Bootstraps micromamba under $PEGRNA_ROOT if no
# conda/mamba/micromamba is found. Everything lives under $PEGRNA_ROOT so nothing
# touches system paths or the (small) home quota.
#
#   PEGRNA_ROOT   storage root      (default /mnt/data2/mj1824/pegrna)
#   PEGRNA_REPO   repo checkout     (default $PEGRNA_ROOT/prime-editing)
#   PY_VER        python version    (default 3.11)
#
# Usage:  bash web/deploy/setup_envs.sh
set -euo pipefail

PEGRNA_ROOT="${PEGRNA_ROOT:-/mnt/data2/mj1824/pegrna}"
PEGRNA_REPO="${PEGRNA_REPO:-$PEGRNA_ROOT/prime-editing}"
PY_VER="${PY_VER:-3.11}"
MAMBA_ROOT="$PEGRNA_ROOT/micromamba"

log() { echo "[setup_envs $(date -u +%H:%M:%S)] $*"; }

mkdir -p "$PEGRNA_ROOT"

# --- resolve a conda-family front-end (prefer existing; else bootstrap mm) ----
MM=""
for c in micromamba mamba conda; do command -v "$c" >/dev/null 2>&1 && { MM="$c"; break; }; done
if [ -z "$MM" ]; then
  log "no conda/mamba/micromamba found — bootstrapping micromamba into $MAMBA_ROOT"
  mkdir -p "$MAMBA_ROOT/bin"
  arch="$(uname -m)"; case "$arch" in x86_64) MA=linux-64;; aarch64) MA=linux-aarch64;; *) MA=linux-64;; esac
  curl -Ls "https://micro.mamba.pm/api/micromamba/${MA}/latest" | tar -xvj -C "$MAMBA_ROOT" bin/micromamba
  MM="$MAMBA_ROOT/bin/micromamba"
fi
export MAMBA_ROOT_PREFIX="$MAMBA_ROOT"
log "using front-end: $MM  (root prefix $MAMBA_ROOT_PREFIX)"

# micromamba/mamba/conda share the create/run verbs used below.
create() { "$MM" create -y -p "$1" -c conda-forge -c bioconda "python=$PY_VER" "${@:2}"; }
run()    { "$MM" run -p "$1" "${@:2}"; }

WEB_ENV="$PEGRNA_ROOT/envs/pegrna-web"
REAL_ENV="$PEGRNA_ROOT/envs/pegrna-real"

# --- pegrna-web (fake mode) --------------------------------------------------
if [ ! -x "$WEB_ENV/bin/python" ]; then
  log "creating pegrna-web"
  create "$WEB_ENV"
fi
log "pegrna-web: install CLI (editable) + web layer"
run "$WEB_ENV" python -m pip install -q --upgrade pip
run "$WEB_ENV" python -m pip install -q -e "$PEGRNA_REPO"
run "$WEB_ENV" python -m pip install -q -r "$PEGRNA_REPO/web/requirements.txt"

# --- pegrna-real (api mode) --------------------------------------------------
if [ ! -x "$REAL_ENV/bin/python" ]; then
  log "creating pegrna-real (+ ViennaRNA from conda-forge)"
  create "$REAL_ENV" viennarna
fi
log "pegrna-real: install CLI + api model stack (CPU) + web layer"
run "$REAL_ENV" python -m pip install -q --upgrade pip
# TensorFlow-CPU first so spliceai doesn't drag in the CUDA build.
run "$REAL_ENV" python -m pip install -q "tensorflow-cpu"
# CLI + the api-relevant extras: splice (SpliceAI), importance (biotite),
# scarless (ViennaRNA python bindings), recode (dnachisel). NOT complex/fullstructure.
run "$REAL_ENV" python -m pip install -q -e "$PEGRNA_REPO[splice,importance,scarless,recode]"
run "$REAL_ENV" python -m pip install -q -r "$PEGRNA_REPO/web/requirements.txt"

# maxentpy (MaxEntScan) — vendored; NOT pip-installable on py>=3.11 because its
# pre-generated Cython _hashseq.c includes the private longintrepr.h. The pure-
# Python maxent.py (score5/score3) is all the pipeline uses, so install WITHOUT
# the C extension: strip ext_modules from setup.py, then --no-build-isolation.
MAXSRC="$PEGRNA_REPO/web/maxentpy_src"
if ! run "$REAL_ENV" python -c "import maxentpy.maxent" >/dev/null 2>&1; then
  log "pegrna-real: install vendored maxentpy (pure-Python, no C ext)"
  TMP="$(mktemp -d)"; cp -r "$MAXSRC/." "$TMP/"
  sed -i 's/ext_modules=ext_modules,//' "$TMP/setup.py" 2>/dev/null || true
  run "$REAL_ENV" python -m pip install -q --no-build-isolation "$TMP"
  rm -rf "$TMP"
fi

log "verify api-mode imports resolve in pegrna-real"
run "$REAL_ENV" python - <<'PY'
import importlib, sys
mods = ["pegrna_pipeline", "RNA", "biotite", "spliceai", "maxentpy.maxent", "fastapi"]
bad = []
for m in mods:
    try:
        importlib.import_module(m)
    except Exception as e:  # noqa: BLE001
        bad.append(f"{m}: {e}")
print("api-mode import check:", "OK" if not bad else "MISSING -> " + "; ".join(bad))
sys.exit(1 if bad else 0)
PY

log "done. envs:"
log "  pegrna-web  -> $WEB_ENV/bin/python"
log "  pegrna-real -> $REAL_ENV/bin/python  (set PEGRNA_API_PYTHON to this)"
