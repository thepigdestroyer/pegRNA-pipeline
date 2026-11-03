#!/usr/bin/env bash
# Stage a GRCh38 genome for the off-target gate so cas-offinder runs for real.
#
# Without a genome, _api_bundle falls back to a clean-pass stand-in: every guide
# is reported as having zero off-targets. That is declared in the report now
# (`stand_in_gates`), but on a production deployment you want the real check --
# an unscreened guide is exactly what this gate exists to catch.
#
# Uses UCSC's .2bit (~950 MB) rather than an unpacked FASTA (~3.1 GB).
# cas-offinder >= 2.4 reads .2bit directly, and the gate only COUNTS hits per
# mismatch tier -- it never reports coordinates -- so UCSC's `chr1` naming
# versus Ensembl's `1` is irrelevant here. On a shared box where the deployment
# lives on the system disk, 2.2 GB of headroom matters more than the convention.
#
# Usage:
#   scripts/install_genome.sh [target-dir]      # default: $PEGRNA_ROOT/genomes
set -euo pipefail

PEGRNA_ROOT="${PEGRNA_ROOT:-$HOME/pegrna}"
TARGET="${1:-$PEGRNA_ROOT/genomes}"
URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit"
GENOME="$TARGET/hg38.2bit"
MIN_BYTES=800000000   # a truncated download must not be accepted silently

log() { echo "[install_genome $(date -u +%H:%M:%S)] $*"; }

mkdir -p "$TARGET"

# Refuse to fill a disk that other people depend on.
avail_kb="$(df -Pk "$TARGET" | tail -1 | awk '{print $4}')"
if [ "$avail_kb" -lt 2000000 ]; then
  echo "error: only $((avail_kb / 1024)) MB free at $TARGET; need ~2 GB." >&2
  echo "       Free space or pass a target on a larger volume." >&2
  exit 1
fi

if [ -s "$GENOME" ] && [ "$(stat -c %s "$GENOME" 2>/dev/null || stat -f %z "$GENOME")" -ge "$MIN_BYTES" ]; then
  log "already present: $GENOME ($(du -h "$GENOME" | cut -f1))"
else
  log "downloading hg38.2bit (~950 MB) -> $TARGET"
  curl -L --fail --continue-at - --retry 5 --retry-delay 10 -o "$GENOME" "$URL"
fi

size="$(stat -c %s "$GENOME" 2>/dev/null || stat -f %z "$GENOME")"
if [ "$size" -lt "$MIN_BYTES" ]; then
  echo "error: $GENOME is only $((size / 1024 / 1024)) MB — download looks truncated." >&2
  exit 1
fi
# .2bit starts with the magic number 0x1A412743.
if [ "$(head -c 4 "$GENOME" | od -An -tx1 | tr -d ' \n')" != "4327411a" ]; then
  echo "error: $GENOME does not have the .2bit magic number." >&2
  exit 1
fi
log "OK — $(du -h "$GENOME" | cut -f1), valid .2bit"

cat <<EOF

Genome staged. Point the pipeline at it:

  PEGRNA_GENOME=$GENOME

Add that to ~/.config/pegrna/pegrna.env and restart:

  systemctl --user restart pegrna-web

'offtarget' then disappears from stand_in_gates in the report and the UI.
EOF
