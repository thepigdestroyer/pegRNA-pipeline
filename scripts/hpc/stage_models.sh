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
# Use the pegrna-recoded env's pip, NOT whatever pip is first on PATH (the GPU
# smoke run revealed RiNALMo/RNAformer landing in the wrong site-packages).
PIP="${PEGRNA_PIP:-$HOME/miniforge3/envs/pegrna-recoded/bin/pip}"
pip() { "$PIP" "$@"; }
# MUST run on a LOGIN NODE (internet). CX3 compute nodes are offline, so every
# download below has to happen here; the GPU smoke job then runs fully offline.
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
# Staged as GRCh38.fa (the name smoke_models.py / CasOffinderScorer expect via
# PEGRNA_GENOME). Guard on the SAME filename we create, else it re-downloads 3 GB.
if [ ! -f "$PEGRNA_DATA/GRCh38.fa" ]; then
  log "Fetching GRCh38 primary assembly (GENCODE) for SpliceAI +/-5kb context + off-target"
  wget -q -O "$PEGRNA_DATA/GRCh38.fa.gz" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
  gunzip -f "$PEGRNA_DATA/GRCh38.fa.gz"
fi

# --- 3. RiNALMo (giga) + SS-head weights (Zenodo 15043668, resolved) ---
if [ ! -d "$PEGRNA_DATA/RiNALMo" ]; then
  log "Cloning RiNALMo (lbcb-sci/RiNALMo) + installing"
  git clone --depth 1 https://github.com/lbcb-sci/RiNALMo "$PEGRNA_DATA/RiNALMo"
  pip install -e "$PEGRNA_DATA/RiNALMo" flash-attn --no-build-isolation || \
    pip install -e "$PEGRNA_DATA/RiNALMo"
fi
mkdir -p "$PEGRNA_DATA/rinalmo_weights"
ZEN_R=https://zenodo.org/records/15043668/files
for w in rinalmo_giga_pretrained.pt rinalmo_giga_ss_archiveII-23s_ft.pt \
         rinalmo_giga_ss_archiveII-srp_ft.pt; do
  [ -f "$PEGRNA_DATA/rinalmo_weights/$w" ] || { log "Fetching $w"; \
    wget -q -O "$PEGRNA_DATA/rinalmo_weights/$w" "$ZEN_R/$w?download=1"; }
done
log "RiNALMo SS checkpoint default: rinalmo_giga_ss_archiveII-23s_ft.pt (broad family; \
validate choice for engineered pegRNA against RNAformer + ViennaRNA)"

# --- 4. RNAformer (CPU-capable fallback; Freiburg artifacts, resolved) ---
if [ ! -d "$PEGRNA_DATA/RNAformer" ]; then
  log "Cloning RNAformer (automl/RNAformer) + installing"
  git clone --depth 1 https://github.com/automl/RNAformer "$PEGRNA_DATA/RNAformer"
  pip install -e "$PEGRNA_DATA/RNAformer"
fi
mkdir -p "$PEGRNA_DATA/rnaformer_weights"
RF=https://ml.informatik.uni-freiburg.de/research-artifacts/RNAformer/models
# BOTH the state_dict AND its matching config .yml are required: RiboFormer is
# built from the config's dims and load_state_dict is strict=True, so a mismatched
# / missing config fails (CX3 smoke: 'NoneType has no attribute items' from the
# template default_config.yaml). Names pair as _state_dict_X.pth <-> _config_X.yml.
for w in RNAformer_32M_state_dict_biophysical.pth \
         RNAformer_32M_config_biophysical.yml \
         RNAformer_32M_state_dict_inter_family_finetuned.pth \
         RNAformer_32M_config_inter_family_finetuned.yml; do
  [ -f "$PEGRNA_DATA/rnaformer_weights/$w" ] || { log "Fetching $w"; \
    wget -q -O "$PEGRNA_DATA/rnaformer_weights/$w" "$RF/$w"; }
done

# --- 5. Boltz-2 (MIT, pip) — weights auto-download on first run (needs internet) ---
log "Installing boltz (pip) — pre-fetch weights on the login node:"
pip install boltz>=2.0
log "Pre-fetching Boltz-2 weights via a tiny login-node prediction (cache to \$BOLTZ_CACHE)"
export BOLTZ_CACHE="$PEGRNA_DATA/boltz_cache"
mkdir -p "$BOLTZ_CACHE"
# A 1-residue dummy predict primes the weight cache; ignore failure if offline-limited.
# msa: empty => single-sequence mode, so priming does not need the MSA server.
printf 'version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: M\n      msa: empty\n' \
  > "$PEGRNA_DATA/_prime.yaml"
boltz predict "$PEGRNA_DATA/_prime.yaml" --out_dir "$BOLTZ_CACHE" --cache "$BOLTZ_CACHE" --no_kernels || \
  log "WARN: prime Boltz cache manually before a compute-node run"

# --- 6. AlphaFold-DB models: pre-fetch on the login node (compute is offline) ---
mkdir -p "$PEGRNA_DATA/afdb"
# Smoke test protein (calmodulin) + any design targets in $PEGRNA_AFDB_ACCS.
for acc in P0DP23 ${PEGRNA_AFDB_ACCS:-}; do
  if ! ls "$PEGRNA_DATA/afdb/AF-$acc-F1-model_v"*.cif >/dev/null 2>&1; then
    url=$(curl -s "https://alphafold.ebi.ac.uk/api/prediction/$acc" \
      | python3 -c "import sys,json;d=json.load(sys.stdin);print(next((m['cifUrl'] for m in d if '-F1-' in m['cifUrl']), d[0]['cifUrl']))" 2>/dev/null)
    [ -n "$url" ] && { log "Fetching AFDB $acc"; wget -q -P "$PEGRNA_DATA/afdb" "$url"; }
  fi
done
log "AFDB models for the smoke protein + \$PEGRNA_AFDB_ACCS staged to $PEGRNA_DATA/afdb"

log "Staging complete. Manifest: $MANIFEST"
echo "Set these for the run: PEGRNA_DATA=$PEGRNA_DATA BOLTZ_CACHE=$BOLTZ_CACHE"
