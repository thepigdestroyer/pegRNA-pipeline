#!/usr/bin/env bash
# One-command HPC install for the recoded TJ-PE pipeline.
#
#   git clone <repo> && cd pegRNA-pipeline
#   export PEGRNA_DATA=$EPHEMERAL/pegrna-data      # where weights/genome are staged
#   bash scripts/hpc/install.sh
#
# Run this ON A CX3 LOGIN NODE (it needs internet to build the env + download
# weights). It is idempotent: re-running skips an existing env and already-staged
# weights. After it finishes, submit a GPU run:
#
#   qsub -v PEGRNA_DATA=$PEGRNA_DATA,CONFIG=examples/ctcf_cterm/config.json \
#        scripts/hpc/run_recoded_pipeline.pbs
#
# Steps:
#   1. create/refresh the conda env from scripts/hpc/environment.yml
#   2. pip-install the package + all real-gate extras (editable)
#   3. stage model weights + reference genome (scripts/hpc/stage_models.sh)
#   4. run the model smoke check (CPU-capable gates on the login node)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

: "${PEGRNA_DATA:?set PEGRNA_DATA to a staging dir, e.g. export PEGRNA_DATA=\$EPHEMERAL/pegrna-data}"
ENV_NAME="${PEGRNA_ENV:-pegrna-recoded}"
CONDA_ROOT="${CONDA_ROOT:-$HOME/miniforge3}"
log() { echo "[$(date -u +%FT%TZ)] $*"; }

# --- locate a conda/mamba ---
if [ ! -x "$CONDA_ROOT/bin/conda" ] && ! command -v conda >/dev/null 2>&1; then
  log "ERROR: no conda found. Install miniforge to \$CONDA_ROOT ($CONDA_ROOT) first:"
  log "  wget -qO- https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh | bash -b -p $CONDA_ROOT"
  exit 1
fi
CONDA_BIN="$CONDA_ROOT/bin/conda"; [ -x "$CONDA_BIN" ] || CONDA_BIN="$(command -v conda)"
MAMBA_BIN="$CONDA_ROOT/bin/mamba"; [ -x "$MAMBA_BIN" ] || MAMBA_BIN="$CONDA_BIN"
# shellcheck disable=SC1091
source "$("$CONDA_BIN" info --base)/etc/profile.d/conda.sh"

# --- 1. env ---
if "$CONDA_BIN" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  log "conda env '$ENV_NAME' already exists (skipping create; use 'conda env remove -n $ENV_NAME' to rebuild)"
else
  log "Creating conda env '$ENV_NAME' from scripts/hpc/environment.yml"
  "$MAMBA_BIN" env create -n "$ENV_NAME" -f scripts/hpc/environment.yml
fi
conda activate "$ENV_NAME"
PIP="$CONDA_ROOT/envs/$ENV_NAME/bin/pip"; [ -x "$PIP" ] || PIP="pip"

# --- 2. package + all real-gate extras (editable) ---
log "Installing pegrna_pipeline + real-gate extras (recode,splice,importance,fullstructure,complex)"
"$PIP" install -e '.[recode,splice,importance,fullstructure,complex]'

# --- 3. stage weights + genome (idempotent) ---
log "Staging model weights + reference genome to \$PEGRNA_DATA=$PEGRNA_DATA"
PEGRNA_PIP="$PIP" bash scripts/hpc/stage_models.sh

# --- 4. smoke check (login-node CPU-capable gates: splice, importance, structure/RNAformer) ---
log "Model smoke check (login-node gates). Boltz/RiNALMo need a GPU node -> run the GPU smoke:"
log "  qsub -v PEGRNA_DATA=$PEGRNA_DATA scripts/hpc/gpu_smoke.pbs"
export RNAFORMER_CKPT="${RNAFORMER_CKPT:-$PEGRNA_DATA/rnaformer_weights/RNAformer_32M_state_dict_biophysical.pth}"
export PEGRNA_GENOME="${PEGRNA_GENOME:-$PEGRNA_DATA/GRCh38.fa}"
python scripts/hpc/smoke_models.py --data "$PEGRNA_DATA" --only splice importance structure || \
  log "WARN: a login-node gate did not pass; inspect output above before a real run"

log "Install complete. Next:"
log "  # GPU smoke (all 5 gates, incl. Boltz + RiNALMo):"
log "  qsub -v PEGRNA_DATA=$PEGRNA_DATA scripts/hpc/gpu_smoke.pbs"
log "  # end-to-end CTCF design (real gates):"
log "  qsub -v PEGRNA_DATA=$PEGRNA_DATA,CONFIG=examples/ctcf_cterm/config.json scripts/hpc/run_recoded_pipeline.pbs"
