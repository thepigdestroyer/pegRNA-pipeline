#!/usr/bin/env bash
# Build the two conda envs the web app needs, on a GPU-less host (e.g. the DSI box).
#
#   pegrna-web   light: the CLI package + FastAPI/uvicorn. Serves `fake` mode.
#   pegrna-real  heavy: adds the CPU/API model stack for `api` mode
#                (ViennaRNA, SpliceAI+TensorFlow-CPU, biotite, dnachisel, maxentpy).
#                NO torch/boltz — docking is the NVIDIA NIM API, not local GPU.
#   pridict2     PRIDICT2 pegRNA designer for `api`-mode auto-spacer design.
#                Cloned from upstream + its own env (python 3.10, CPU torch==2.0.1).
#                Torch lives ONLY here; pegrna-real stays torch-free.
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
# TensorFlow-CPU first so spliceai doesn't drag in the CUDA build. setuptools<81
# restores pkg_resources (removed in setuptools>=81), which spliceai imports.
run "$REAL_ENV" python -m pip install -q "tensorflow-cpu" "setuptools<81"
# CLI + the api-relevant extras: splice (SpliceAI), importance (biotite),
# scarless (ViennaRNA python bindings), recode (dnachisel). NOT complex/fullstructure.
run "$REAL_ENV" python -m pip install -q -e "$PEGRNA_REPO[splice,importance,scarless,recode]"
run "$REAL_ENV" python -m pip install -q -r "$PEGRNA_REPO/web/requirements.txt"

# maxentpy (MaxEntScan) — vendored, installed without its C extension. The logic
# lives in scripts/install_maxentpy.sh so this and a plain local install cannot
# drift apart; it used to exist only here, which is why api mode failed at the
# cryptic-splice gate on every host except this one.
log "pegrna-real: install vendored maxentpy (pure-Python, no C ext)"
PYTHON="$REAL_ENV/bin/python" bash "$PEGRNA_REPO/scripts/install_maxentpy.sh"

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

# --- pridict2 (api mode: auto-spacer pegRNA design, CPU) ---------------------
# PRIDICT2 designs the pegRNA (spacer/PBS/RTT) at the insert site. Cloned WITH its
# trained_models (~95M). Env from the repo's own pridict2_repo.yml (python 3.10),
# then CPU torch==2.0.1 (the build PRIDICT2 officially supports).
PRIDICT_DIR="$PEGRNA_ROOT/PRIDICT2"
PRIDICT_ENV="$PEGRNA_ROOT/envs/pridict2"

if [ ! -d "$PRIDICT_DIR/.git" ]; then
  log "cloning PRIDICT2 into $PRIDICT_DIR (includes trained_models, ~95M)"
  git clone https://github.com/uzh-dqbm-cmi/PRIDICT2.git "$PRIDICT_DIR"
else
  log "PRIDICT2 already present at $PRIDICT_DIR"
fi

if [ ! -x "$PRIDICT_ENV/bin/python" ]; then
  log "creating pridict2 env from $PRIDICT_DIR/pridict2_repo.yml (python 3.10)"
  "$MM" create -y -p "$PRIDICT_ENV" -f "$PRIDICT_DIR/pridict2_repo.yml"
fi
log "pridict2: install CPU torch==2.0.1"
run "$PRIDICT_ENV" python -m pip install -q "torch==2.0.1" \
  --index-url https://download.pytorch.org/whl/cpu

log "verify pridict2 env: torch (CPU), Bio, pandas import"
run "$PRIDICT_ENV" python -c \
  "import torch,Bio,pandas; print(torch.__version__, torch.cuda.is_available())"

log "done. envs:"
log "  pegrna-web  -> $WEB_ENV/bin/python"
log "  pegrna-real -> $REAL_ENV/bin/python  (set PEGRNA_API_PYTHON to this)"
log "  pridict2    -> $PRIDICT_ENV/bin/python  (set PRIDICT2_PYTHON to this)"
log ""
log "The frontend is served from web/static/, which is committed pre-built —"
log "no Node is needed here. Rebuild only from web/frontend/ (see its README)."
