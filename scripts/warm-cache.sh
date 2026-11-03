#!/usr/bin/env bash
# Run the shipped example once so its off-target results are cached.
#
# The first design for a set of guides runs a whole-genome cas-offinder scan (~69 CPU-minutes,
# a few minutes of wall time at the default thread cap). Results are cached per guide, so
# every later run of the same guides is fast. Do this after a deploy, or after replacing the
# genome (which invalidates the cache), so the first visitor to try the demo does not pay.
#
#   scripts/warm-cache.sh            # on the backend host
#   SSH_HOST=dsi scripts/warm-cache.sh   # from a laptop
#
# Needs the api-mode token if the deployment is token-gated: PEGRNA_API_MODE_TOKEN in the
# environment, or in ~/.config/pegrna/pegrna.env on the backend host.
set -euo pipefail

PORT="${PEGRNA_PORT:-8093}"
REPO="${PEGRNA_REPO:-$HOME/pegrna/prime-editing}"

run() {
  local token="${PEGRNA_API_MODE_TOKEN:-}"
  if [ -z "$token" ] && [ -f "$HOME/.config/pegrna/pegrna.env" ]; then
    token="$(grep '^PEGRNA_API_MODE_TOKEN=' "$HOME/.config/pegrna/pegrna.env" | cut -d= -f2 || true)"
  fi
  local started; started="$(date +%s)"
  echo "running the example design (api mode) to warm the off-target cache..."
  curl -sS --max-time 3000 -X POST "http://127.0.0.1:$PORT/api/design" \
    -H 'Content-Type: application/json' \
    ${token:+-H "X-Access-Token: $token"} \
    -d "{\"mode\":\"api\",\"config\":$(cat "$REPO/examples/ctcf_cterm/config_autospacer_flag.json")}" \
    | python3 -c 'import json,sys; r=json.load(sys.stdin); print("ok:", r["ok"], "| stand-in gates:", r["gates"].get("stand_in_gates"))'
  echo "done in $(( $(date +%s) - started )) s. Re-run it: it should now be much faster."
}

if [ -n "${SSH_HOST:-}" ]; then
  ssh "$SSH_HOST" "PEGRNA_PORT=$PORT PEGRNA_REPO=$REPO bash -s" < "$0"
else
  run
fi
