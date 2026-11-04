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
# PEGRNA_SKIP_PULL=1 uses whatever is already in $PEGRNA_REPO. That covers the
# case where the tree was pushed to the host directly (rsync/scp) rather than
# cloned -- a private repo would otherwise need a PAT or deploy key stored on a
# shared box purely so deploy.sh can re-fetch what is already there.
if [ "${PEGRNA_SKIP_PULL:-0}" = "1" ]; then
  if [ ! -d "$PEGRNA_REPO" ]; then
    echo "error: PEGRNA_SKIP_PULL=1 but $PEGRNA_REPO does not exist" >&2
    exit 1
  fi
  log "PEGRNA_SKIP_PULL=1 — using the tree already at $PEGRNA_REPO"
elif [ -d "$PEGRNA_REPO/.git" ]; then
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
  # Mode 600 from the start: the access token is appended to this file later, and the host is
  # shared with group-readable home directories.
  install -m 600 "$PEGRNA_REPO/web/deploy/pegrna.env" "$HOME/.config/pegrna/pegrna.env"
fi
# Tighten an older install too: nothing in this directory is for anyone but this user.
chmod 700 "$HOME/.config/pegrna"
find "$HOME/.config/pegrna" -type f -exec chmod 600 {} +
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
