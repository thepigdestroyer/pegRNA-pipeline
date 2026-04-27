#!/bin/bash
# driver_pbs_wrapper.sh — PBS wrapper around the pegRNA-pipeline driver.
#
# Submit with:
#     qsub -v RUN_DIR=/path/to/run driver_pbs_wrapper.sh
#
# The script:
#   1. Loads miniforge/3 and activates the pegrna-pipeline conda env.
#   2. cds to the pegRNA-pipeline git checkout.
#   3. Runs ``python -m pegrna_pipeline.hpc.driver "$RUN_DIR"``.
#   4. On success, echoes a DRIVER_DONE marker line (the launcher greps
#      for this in the job log to decide whether to download reports).
#
# Configure the PBS resources above this comment at submit-time (the
# launcher renders them dynamically via render_pbs_script).
set -euo pipefail

if [[ -z "${RUN_DIR:-}" ]]; then
    echo "ERROR: RUN_DIR environment variable must be set (use qsub -v RUN_DIR=...)" >&2
    exit 2
fi

module load miniforge/3 || true
# Imperial CX3's miniforge module only ships a bootstrap; the real conda
# install lives in ~/miniforge3 (created by setup_hpc.sh).
# shellcheck disable=SC1091
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate pegrna-pipeline

REPO_DIR="${PEGRNA_REPO_DIR:-$HOME/pegRNA-pipeline}"
cd "$REPO_DIR"

python -m pegrna_pipeline.hpc.driver "$RUN_DIR"

echo "DRIVER_DONE"
