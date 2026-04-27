#!/usr/bin/env bash
# setup_hpc.sh — bootstrap pegRNA-pipeline on Imperial CX3.
#
# Run from a CX3 login node:
#   bash setup_hpc.sh <GITHUB_PAT>
#
# Takes ~30-60 min (conda + Bioconductor installs).

set -euo pipefail

GITHUB_PAT="${1:?ERROR: pass your GitHub PAT as the first argument}"

REPO_URL="https://${GITHUB_PAT}@github.com/thepigdestroyer/pegRNA-pipeline.git"
REPO_DIR="${HOME}/pegRNA-pipeline"
PRIDICT_DIR="${HOME}/PRIDICT2"
RUNS_DIR="${EPHEMERAL}/pegrna-runs"

echo "╔══════════════════════════════════════════════════════╗"
echo "║      pegRNA-pipeline HPC bootstrap (Imperial CX3)   ║"
echo "╚══════════════════════════════════════════════════════╝"
echo "HOME:       $HOME"
echo "EPHEMERAL:  $EPHEMERAL"
echo ""

# ── 0. conda init ──────────────────────────────────────────────────────────────
# Imperial CX3 ships only a `miniforge-setup` bootstrapper; the per-user
# install lives in $HOME/miniforge3 and must be run once before `conda` is
# available on PATH.
module load miniforge/3
if [[ ! -x "$HOME/miniforge3/bin/conda" ]]; then
    echo "      installing miniforge into \$HOME/miniforge3 (one-time, ~3 min)..."
    miniforge-setup
fi
# shellcheck disable=SC1091
source "$HOME/miniforge3/etc/profile.d/conda.sh"

# ── 1. Clone / update pipeline repo ───────────────────────────────────────────
echo "[1/6] pegRNA-pipeline source..."
if [[ -d "$REPO_DIR/.git" ]]; then
    git -C "$REPO_DIR" remote set-url origin "$REPO_URL"
    git -C "$REPO_DIR" pull --ff-only
    echo "      updated"
else
    git clone "$REPO_URL" "$REPO_DIR"
    echo "      cloned"
fi

# ── 2. pegrna-pipeline env (driver + codon expander + reports) ─────────────────
echo "[2/6] pegrna-pipeline conda env..."
if conda env list | grep -qE "^pegrna-pipeline "; then
    echo "      already exists — reinstalling package"
else
    conda env create -f "$REPO_DIR/envs/pipeline.yml"
fi
conda activate pegrna-pipeline
pip install -e "$REPO_DIR" -q
conda deactivate
echo "      done"

# ── 3. Clone PRIDICT2 (public repo) ───────────────────────────────────────────
echo "[3/6] PRIDICT2 source..."
if [[ -d "$PRIDICT_DIR/.git" ]]; then
    git -C "$PRIDICT_DIR" pull --ff-only
    echo "      updated"
else
    git clone https://github.com/uzh-dqbm-cmi/PRIDICT2.git "$PRIDICT_DIR"
    echo "      cloned"
fi

# ── 4. pridict2 conda env ─────────────────────────────────────────────────────
echo "[4/6] pridict2 conda env (~15 min)..."
if conda env list | grep -qE "^pridict2 "; then
    echo "      already exists, skipping"
else
    conda env create -n pridict2 -f "$PRIDICT_DIR/pridict2_repo.yml"
    conda activate pridict2
    # PyTorch CPU build must be installed separately per PRIDICT2 README
    pip install torch==2.0.1 --index-url https://download.pytorch.org/whl/cpu -q
    conda deactivate
fi
echo "      done"

# ── 5. pegrna-ngrna env (R + crisprScore via Bioconductor) ────────────────────
echo "[5/6] pegrna-ngrna conda env + Bioconductor (~30 min)..."
if conda env list | grep -qE "^pegrna-ngrna "; then
    echo "      already exists — checking Bioconductor packages"
else
    conda env create -f "$REPO_DIR/envs/offtarget.yml"
fi
conda activate pegrna-ngrna
Rscript - <<'REOF'
suppressMessages(library(BiocManager))
pkgs <- c("crisprScore", "crisprDesign", "Biostrings", "BSgenome.Hsapiens.UCSC.hg38")
needed <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(needed) > 0) {
    message("Installing: ", paste(needed, collapse = ", "))
    BiocManager::install(needed, ask = FALSE, update = FALSE)
} else {
    message("All Bioconductor packages already present")
}
REOF
conda deactivate
echo "      done"

# ── 6. Ephemeral runs directory ────────────────────────────────────────────────
echo "[6/6] Ephemeral runs directory..."
mkdir -p "$RUNS_DIR"
echo "      $RUNS_DIR"

# ── Summary ────────────────────────────────────────────────────────────────────
echo ""
echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║  Setup complete. Use these values when running 'pegrna setup'       ║"
echo "║  on your laptop (or they are already written to ~/.pegrna/config):  ║"
echo "╠══════════════════════════════════════════════════════════════════════╣"
printf "║  %-22s %s\n" "hostname:"             "login.cx3.hpc.imperial.ac.uk                  ║"
printf "║  %-22s %s\n" "username:"             "mj1824                                        ║"
printf "║  %-22s %s\n" "repo_dir:"             "${REPO_DIR}                     ║"
printf "║  %-22s %s\n" "remote_runs_base:"     "${RUNS_DIR}         ║"
printf "║  %-22s %s\n" "entrez_email:"         "mrinal.joshi@ukdri.ac.uk                      ║"
printf "║  %-22s %s\n" "pridict_install_path:" "${PRIDICT_DIR}                   ║"
printf "║  %-22s %s\n" "rscript_path:"         "${REPO_DIR}/scripts/score_guides.R  ║"
printf "║  %-22s %s\n" "conda_env:"            "pridict2                                      ║"
printf "║  %-22s %s\n" "compute_mode:"         "local                                         ║"
printf "║  %-22s %s\n" "cell_line:"            "hek293t                                       ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
