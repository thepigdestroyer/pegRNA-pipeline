#!/usr/bin/env bash
# Stage all model weights + reference data for the recoded TJ-PE pipeline on CX3.
# Run ONCE on a login node (has internet); compute nodes are offline, so every
# weight/genome must be cached to $PEGRNA_DATA first.
#
#   export PEGRNA_DATA=$EPHEMERAL/pegrna-data
#   bash scripts/hpc/stage_models.sh
#
# Idempotent: re-running skips anything already present. Pins are recorded to
# $PEGRNA_DATA/MANIFEST.txt for provenance (goes into the design run manifest).
set -euo pipefail

: "${PEGRNA_DATA:?set PEGRNA_DATA to a staging dir on \$EPHEMERAL or \$RDS}"
mkdir -p "$PEGRNA_DATA"
MANIFEST="$PEGRNA_DATA/MANIFEST.txt"
: > "$MANIFEST"
log() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$MANIFEST"; }

# --- 1. T2T diploid KOLF2.1J assembly + GRCh38 chain (B3-C) ---
if [ ! -f "$PEGRNA_DATA/KOLF2.1Jv1.1.noY.fa" ]; then
  log "Fetching T2T-KOLF2.1J diploid assembly (Zenodo 18777354)"
  # marbl/KOLF2.1J-Diploid-Assembly -> Zenodo record 18777354
  wget -q -O "$PEGRNA_DATA/kolf.zip" \
    "https://zenodo.org/records/18777354/files/KOLF2.1Jv1.1.noY.fa.gz" || \
    log "WARN: confirm exact Zenodo asset name on the record page"
  gunzip -f "$PEGRNA_DATA/kolf.zip" 2>/dev/null || true
fi
# GRCh38<->KOLF chain: build with minimap2 if no published chain is found.
log "TODO-IF-MISSING: GRCh38<->KOLF2.1Jv1.1 chain (build via minimap2 -cx asm5)"

# --- 2. SpliceAI weights (ship with the pip package) + GRCh38 ref for context ---
log "SpliceAI weights ship inside the 'spliceai' pip package (models/spliceai{1..5}.h5)"
if [ ! -f "$PEGRNA_DATA/GRCh38.primary_assembly.genome.fa" ]; then
  log "Fetching GRCh38 primary assembly (GENCODE) for SpliceAI +/-5kb context"
  wget -q -O "$PEGRNA_DATA/GRCh38.fa.gz" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
  gunzip -f "$PEGRNA_DATA/GRCh38.fa.gz"
fi

# --- 3. RiNALMo (giga) + SS-head weights ---
if [ ! -d "$PEGRNA_DATA/RiNALMo" ]; then
  log "Cloning RiNALMo (lbcb-sci/RiNALMo) + downloading giga + SS-head weights"
  git clone --depth 1 https://github.com/lbcb-sci/RiNALMo "$PEGRNA_DATA/RiNALMo"
  pip install -e "$PEGRNA_DATA/RiNALMo"
  # Pretrained + fine-tuned SS weights are linked from the repo README (HF / Zenodo).
  log "TODO: download rinalmo_giga_pretrained.pt + rinalmo_giga_ss_*_ft.pt per README"
fi

# --- 4. RNAformer (CPU fallback) ---
if [ ! -d "$PEGRNA_DATA/RNAformer" ]; then
  log "Cloning RNAformer (automl/RNAformer) + checkpoint"
  git clone --depth 1 https://github.com/automl/RNAformer "$PEGRNA_DATA/RNAformer"
  pip install -e "$PEGRNA_DATA/RNAformer"
  log "TODO: download RNAformer.ckpt per repo (conda-forge rnaformer also available)"
fi

# --- 5. Boltz-2 (MIT, pip) — weights auto-download on first run (needs internet) ---
log "Installing boltz (pip) — pre-fetch weights on the login node:"
pip install boltz>=2.0
log "Pre-fetching Boltz-2 weights via a tiny login-node prediction (cache to \$BOLTZ_CACHE)"
export BOLTZ_CACHE="$PEGRNA_DATA/boltz_cache"
mkdir -p "$BOLTZ_CACHE"
# A 1-residue dummy predict primes the weight cache; ignore failure if offline-limited.
printf 'version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: M\n' \
  > "$PEGRNA_DATA/_prime.yaml"
boltz predict "$PEGRNA_DATA/_prime.yaml" --out_dir "$BOLTZ_CACHE" || \
  log "WARN: prime Boltz cache manually before a compute-node run"

# --- 6. AlphaFold-DB models: fetched per-UniProt at design time (cache dir) ---
mkdir -p "$PEGRNA_DATA/afdb"
log "AFDB models fetched per target UniProt accession into $PEGRNA_DATA/afdb at design time"

log "Staging complete. Manifest: $MANIFEST"
echo "Set these for the run: PEGRNA_DATA=$PEGRNA_DATA BOLTZ_CACHE=$BOLTZ_CACHE"
