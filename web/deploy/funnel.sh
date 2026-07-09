#!/usr/bin/env bash
# Expose the loopback web app publicly via Tailscale Funnel — WITHOUT opening any
# inbound port on the host (the tunnel is outbound to Tailscale's relays). This is
# the security-safe hosting path for a shared DSI compute box.
#
# Runs tailscaled in USERSPACE mode (no root, no TUN device) under $PEGRNA_ROOT,
# so it works as an unprivileged user. Funnel must be permitted for this node in
# the tailnet ACLs (admin console: Access controls -> nodeAttrs "funnel"), and
# HTTPS/MagicDNS enabled for the tailnet.
#
#   PEGRNA_ROOT   storage root  (default /mnt/data2/mj1824/pegrna)
#   PEGRNA_PORT   loopback port (default 8092)
#
# Usage:
#   bash web/deploy/funnel.sh up      # start tailscaled(userspace) + funnel :443 -> app
#   bash web/deploy/funnel.sh status  # show funnel URL + state
#   bash web/deploy/funnel.sh down     # stop funnel
set -euo pipefail

PEGRNA_ROOT="${PEGRNA_ROOT:-/mnt/data2/mj1824/pegrna}"
PEGRNA_PORT="${PEGRNA_PORT:-8092}"
TS_DIR="$PEGRNA_ROOT/tailscale"
TS_SOCK="$TS_DIR/tailscaled.sock"
TS_STATE="$TS_DIR/tailscaled.state"
mkdir -p "$TS_DIR"
log() { echo "[funnel $(date -u +%H:%M:%S)] $*"; }

# Resolve tailscale/tailscaled (system install preferred; else fetch static bins).
TS="$(command -v tailscale || true)"
TSD="$(command -v tailscaled || true)"
if [ -z "$TS" ] || [ -z "$TSD" ]; then
  log "tailscale not on PATH — fetching static binaries into $TS_DIR"
  arch="$(uname -m)"; case "$arch" in x86_64) A=amd64;; aarch64) A=arm64;; *) A=amd64;; esac
  ver="1.78.1"
  curl -Ls "https://pkgs.tailscale.com/stable/tailscale_${ver}_${A}.tgz" | tar -xz -C "$TS_DIR" --strip-components=1
  TS="$TS_DIR/tailscale"; TSD="$TS_DIR/tailscaled"
fi
TSC="$TS --socket=$TS_SOCK"

cmd="${1:-up}"
case "$cmd" in
  up)
    if ! pgrep -f "tailscaled.*$TS_SOCK" >/dev/null 2>&1; then
      log "starting tailscaled (userspace-networking)"
      nohup "$TSD" --tun=userspace-networking --socket="$TS_SOCK" \
        --statedir="$TS_DIR" >"$TS_DIR/tailscaled.log" 2>&1 &
      sleep 3
    fi
    log "bringing the node up (a browser login URL prints if not yet authed)"
    $TSC up --hostname=pegrna-dsi || true
    log "enabling Funnel: https -> 127.0.0.1:$PEGRNA_PORT"
    $TSC funnel --bg "$PEGRNA_PORT"
    $TSC funnel status || true
    log "if you see a public https URL above, the site is live."
    ;;
  status) $TSC status || true; $TSC funnel status || true ;;
  down)   $TSC funnel reset || true; log "funnel stopped" ;;
  *) echo "usage: $0 {up|status|down}"; exit 2 ;;
esac
