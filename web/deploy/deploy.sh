#!/usr/bin/env bash
# One-command deploy of the pegRNA web app on a GPU-less host (the DSI box).
# Orchestrates: clone/pull -> build envs -> install systemd --user unit + env
# file -> (re)start the service -> health check. Public exposure is a separate,
# explicit step: `web/deploy/funnel.sh up`.
#
#   PEGRNA_ROOT  storage root (default /mnt/data2/mj1824/pegrna)
#   PEGRNA_PORT  loopback port (default 8092)
#   REPO_URL     git remote   (default the neurogenomics repo)
#
# Usage:  bash web/deploy/deploy.sh
set -euo pipefail

PEGRNA_ROOT="${PEGRNA_ROOT:-/mnt/data2/mj1824/pegrna}"
PEGRNA_REPO="${PEGRNA_REPO:-$PEGRNA_ROOT/prime-editing}"
PEGRNA_PORT="${PEGRNA_PORT:-8092}"
REPO_URL="${REPO_URL:-https://github.com/neurogenomics/prime-editing.git}"
export PEGRNA_ROOT PEGRNA_REPO PEGRNA_PORT
log() { echo "[deploy $(date -u +%H:%M:%S)] $*"; }

mkdir -p "$PEGRNA_ROOT"

# 1) checkout ----------------------------------------------------------------
if [ -d "$PEGRNA_REPO/.git" ]; then
  log "updating $PEGRNA_REPO"; git -C "$PEGRNA_REPO" pull --ff-only
else
  log "cloning $REPO_URL -> $PEGRNA_REPO"; git clone "$REPO_URL" "$PEGRNA_REPO"
fi

# 2) envs --------------------------------------------------------------------
log "building conda envs (setup_envs.sh)"
bash "$PEGRNA_REPO/web/deploy/setup_envs.sh"

# 3) config + systemd --user -------------------------------------------------
mkdir -p "$HOME/.config/systemd/user" "$HOME/.config/pegrna"
if [ ! -f "$HOME/.config/pegrna/pegrna.env" ]; then
  log "installing pegrna.env (edit it if your paths differ)"
  cp "$PEGRNA_REPO/web/deploy/pegrna.env" "$HOME/.config/pegrna/pegrna.env"
fi
cp "$PEGRNA_REPO/web/deploy/pegrna-web.service" "$HOME/.config/systemd/user/"
if [ ! -f "$HOME/.config/pegrna/nvidia.key" ]; then
  log "NOTE: no NVIDIA key at ~/.config/pegrna/nvidia.key — api mode will stub docking."
  log "      Add it with: printf %s 'nvapi-...' > ~/.config/pegrna/nvidia.key && chmod 600 ~/.config/pegrna/nvidia.key"
fi

loginctl enable-linger "$USER" 2>/dev/null || log "WARN: could not enable linger (service won't survive logout)"
systemctl --user daemon-reload
systemctl --user enable --now pegrna-web
sleep 3

# 4) health ------------------------------------------------------------------
log "health check"
if curl -sf "http://127.0.0.1:$PEGRNA_PORT/api/health" >/dev/null; then
  curl -s "http://127.0.0.1:$PEGRNA_PORT/api/health"
  echo
  log "app is UP on 127.0.0.1:$PEGRNA_PORT"
  log "expose publicly (no open ports):  bash $PEGRNA_REPO/web/deploy/funnel.sh up"
else
  log "health check FAILED — logs:"
  systemctl --user status pegrna-web --no-pager -l | tail -20 || true
  exit 1
fi
