#!/usr/bin/env bash
# Stage a GRCh38 genome FASTA so the off-target gate runs for real.
#
# Without a FASTA, _api_bundle falls back to a clean-pass stand-in: every guide
# is reported as having zero off-targets. That verdict is now declared in the
# report (`stand_in_gates`), but on a production deployment you want the real
# thing -- an unchecked guide is exactly what the gate exists to catch.
#
# Pairs with the `cas-offinder` binary, installed into the pegrna-real env by
# web/deploy/setup_envs.sh.
#
# Downloads the Ensembl GRCh38 primary assembly (~900 MB gzipped, ~3.1 GB
# unpacked), matching the GRCh38 coordinate frame the rest of the gates use.
#
# Usage:
#   scripts/install_genome.sh [target-dir]     # default: $PEGRNA_ROOT/genomes
#
# Then set in web/deploy/pegrna.env:
#   PEGRNA_GENOME=<target-dir>/GRCh38.primary_assembly.genome.fa
set -euo pipefail

PEGRNA_ROOT="${PEGRNA_ROOT:-/mnt/data2/mj1824/pegrna}"
TARGET="${1:-$PEGRNA_ROOT/genomes}"
URL="https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
FASTA="$TARGET/GRCh38.primary_assembly.genome.fa"

log() { echo "[install_genome $(date -u +%H:%M:%S)] $*"; }

mkdir -p "$TARGET"

if [ -s "$FASTA" ]; then
  log "already present: $FASTA ($(du -h "$FASTA" | cut -f1))"
else
  log "downloading GRCh38 primary assembly (~900 MB gz) -> $TARGET"
  # --continue so an interrupted transfer resumes instead of restarting.
  curl -L --continue-at - --retry 5 --retry-delay 10 -o "$FASTA.gz" "$URL"
  log "decompressing (~3.1 GB)"
  gunzip -c "$FASTA.gz" > "$FASTA"
  rm -f "$FASTA.gz"
fi

log "verifying it looks like a FASTA"
head -c 1 "$FASTA" | grep -q '>' || { echo "error: $FASTA does not start with '>'" >&2; exit 1; }
contigs=$(grep -c '^>' "$FASTA" || true)
log "OK — $contigs sequences"

if command -v cas-offinder >/dev/null 2>&1; then
  log "cas-offinder found: $(command -v cas-offinder)"
else
  log "WARN: cas-offinder is not on PATH here. It lives in the pegrna-real env;"
  log "      the web app invokes it from there, so this warning is fine locally."
fi

cat <<EOF

Genome staged. Point the pipeline at it:

  PEGRNA_GENOME=$FASTA

Add that line to ~/.config/pegrna/pegrna.env and restart:

  systemctl --user restart pegrna-web

The off-target gate then runs cas-offinder for real, and 'offtarget' disappears
from stand_in_gates in the report and the UI.
EOF
