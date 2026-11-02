#!/usr/bin/env bash
# Install PRIDICT2 + its own Python environment for the `api`-mode auto-spacer path.
#
# PRIDICT2 designs the pegRNA#1 spacer (and the nicking guides the ngRNA is chosen
# from) locally on CPU. It needs its OWN environment: it pins torch, and its
# bundled DeepCas9 is TF1 graph-mode code that only runs under Keras 2. Keeping it
# separate is why the main pegrna env stays torch-free.
#
# Two version constraints are load-bearing and were both found the hard way:
#   * numpy < 2        -- torch wheels below 2.3 are built against the numpy 1 ABI,
#                         so numpy 2 breaks `import torch` with "_ARRAY_API not found".
#   * tf-keras         -- DeepCas9 passes a symbolic tf.bool placeholder as
#                         `training=` to tf.keras.layers.Dropout. Keras 3 rejects
#                         it; PRIDICT then exits 0 with no pegRNA CSV, which the
#                         pipeline can only report as "outside its design domain".
#                         PridictLocalRunner sets TF_USE_LEGACY_KERAS=1 for the
#                         subprocess, which needs this package present.
#
# Usage:
#   scripts/install_pridict2.sh [target-dir]      # default: ../PRIDICT2
#
# Then export (or put in web/deploy/pegrna.env):
#   PRIDICT2_DIR=<target-dir>
#   PRIDICT2_PYTHON=<target-dir>/.venv/bin/python
#   PRIDICT2_CORES=3            # PRIDICT caps at 3 to bound memory
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="${1:-$(dirname "$REPO_ROOT")/PRIDICT2}"
PRIDICT_REPO="https://github.com/uzh-dqbm-cmi/PRIDICT2.git"

log() { printf '\033[1;34m==>\033[0m %s\n' "$*"; }

if ! command -v uv >/dev/null 2>&1; then
  echo "error: this script uses 'uv' to build the environment (https://docs.astral.sh/uv/)" >&2
  exit 1
fi

if [ ! -d "$TARGET/.git" ]; then
  log "cloning PRIDICT2 into $TARGET (~127 MB, includes trained weights)"
  # Retry: the weights make this a large fetch and GitHub drops it intermittently.
  for attempt in 1 2 3 4 5; do
    rm -rf "$TARGET"
    if git -c http.postBuffer=524288000 -c http.lowSpeedLimit=0 \
           clone --depth 1 "$PRIDICT_REPO" "$TARGET"; then
      break
    fi
    echo "  clone attempt $attempt failed; retrying..." >&2
    sleep 10
  done
  [ -d "$TARGET/.git" ] || { echo "error: could not clone PRIDICT2" >&2; exit 1; }
else
  log "PRIDICT2 already present at $TARGET"
fi

log "building the pridict2 environment (python 3.10)"
uv venv --python 3.10 "$TARGET/.venv"
uv pip install --python "$TARGET/.venv/bin/python" \
  "torch==2.2.2" \
  "numpy<2" \
  "tensorflow>=2.13" \
  "tf-keras" \
  "biopython==1.81" \
  joblib matplotlib "pandas>=2.0" prettytable primer3-py scikit-learn scipy seaborn tqdm

log "verifying the PRIDICT2 CLI loads"
TF_CPP_MIN_LOG_LEVEL=3 TF_USE_LEGACY_KERAS=1 \
  "$TARGET/.venv/bin/python" "$TARGET/pridict2_pegRNA_design.py" --help >/dev/null

cat <<EOF

PRIDICT2 ready. Export these (or add them to web/deploy/pegrna.env):

  export PRIDICT2_DIR="$TARGET"
  export PRIDICT2_PYTHON="$TARGET/.venv/bin/python"
  export PRIDICT2_CORES=3

Smoke-test the auto-spacer path:

  pegrna design-recoded --config examples/ctcf_cterm/config_autospacer_flag.json \\
      --mode api --output-dir /tmp/pegrna-autospacer
EOF
