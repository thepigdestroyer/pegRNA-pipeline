#!/usr/bin/env bash
# Re-point the public site at the backend's CURRENT tunnel URL, in one command.
#
# A Cloudflare quick tunnel gets a new random hostname every time it restarts,
# and the frontend bakes that hostname in at build time. So a tunnel restart
# silently breaks the site: it loads, then every API call fails. This reads the
# live URL off the backend host, updates CORS there, rebuilds and redeploys.
#
# Run it whenever the site loads but cannot reach the backend.
#
#   scripts/refresh-backend.sh            # preview deployment
#   scripts/refresh-backend.sh --prod     # production
#
#   SSH_HOST   ssh alias/host for the backend   (default: dsi)
#   SITE       public site origin, for CORS     (default: https://pegrna-designer.vercel.app)
#
# The permanent fix is a NAMED Cloudflare tunnel on a domain you own, which
# keeps one hostname forever: see web/deploy/GO-LIVE.md.
set -euo pipefail

SSH_HOST="${SSH_HOST:-dsi}"
SITE="${SITE:-https://pegrna-designer.vercel.app}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

log() { printf '\033[1;34m==>\033[0m %s\n' "$*"; }

log "reading the live tunnel URL from $SSH_HOST"
# tunnel-run.sh publishes whichever provider actually came up (cloudflared, or
# localhost.run when Cloudflare is rate-limiting), so read the published URL
# rather than grepping for one provider's hostname pattern.
URL="$(ssh "$SSH_HOST" 'cat ~/pegrna/tunnel-url.txt 2>/dev/null' | tr -d '[:space:]')"

if [ -z "$URL" ]; then
  echo "error: no tunnel URL published on $SSH_HOST (~/pegrna/tunnel-url.txt is empty)." >&2
  echo "       Check the tunnel:  ssh $SSH_HOST systemctl --user status pegrna-tunnel" >&2
  echo "       Recent attempts:   ssh $SSH_HOST tail -20 ~/pegrna/cloudflared.log" >&2
  exit 1
fi
log "tunnel: $URL"

log "checking the backend answers on it"
if ! curl -sf --max-time 30 "$URL/api/health" >/dev/null; then
  echo "error: $URL/api/health did not respond — not deploying a broken site." >&2
  exit 1
fi

log "allow-listing $SITE for CORS on the backend"
ssh "$SSH_HOST" "
  set -e
  sed -i '/^PEGRNA_CORS_ORIGINS=/d' ~/.config/pegrna/pegrna.env
  echo 'PEGRNA_CORS_ORIGINS=$SITE' >> ~/.config/pegrna/pegrna.env
  systemctl --user restart pegrna-web
"
sleep 4
ssh "$SSH_HOST" "systemctl --user is-active pegrna-web"

log "rebuilding + redeploying the site"
PEGRNA_API_BASE="$URL" bash "$REPO_ROOT/web/vercel/deploy.sh" "${1:-}"

log "verifying the public site can reach the backend"
if curl -sf --max-time 30 -H "Origin: $SITE" "$URL/api/health" | grep -q '"status"'; then
  log "OK — $SITE is pointed at $URL"
else
  echo "warning: the backend did not answer a CORS preflight-style probe." >&2
fi
