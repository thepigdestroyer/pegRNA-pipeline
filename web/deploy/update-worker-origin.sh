#!/usr/bin/env bash
# Point the Cloudflare Worker at the backend's CURRENT tunnel URL.
#
# The public hostname (https://pegrna.ukdri.org) never changes; the tunnel
# behind it does. This re-uploads the Worker with the new origin, which takes a
# second or two and needs no redeploy of the site. Run by the backend host
# itself whenever its tunnel rotates, so the site never goes dark.
#
# Credentials: a Cloudflare API token with Workers Scripts:Edit, read from
# $CF_TOKEN_FILE (default ~/.config/pegrna/cloudflare.token, mode 600).
#
#   web/deploy/update-worker-origin.sh [origin-url]
#
# With no argument it reads $PEGRNA_ROOT/tunnel-url.txt.
set -euo pipefail

PEGRNA_ROOT="${PEGRNA_ROOT:-$HOME/pegrna}"
CF_TOKEN_FILE="${CF_TOKEN_FILE:-$HOME/.config/pegrna/cloudflare.token}"
CF_ACCOUNT="${CF_ACCOUNT:-03adf89089ebc97cadf187a3de136e09}"
WORKER_NAME="${WORKER_NAME:-pegrna-proxy}"
SCRIPT_SRC="${SCRIPT_SRC:-$PEGRNA_ROOT/prime-editing/web/deploy/cloudflare-worker.js}"

ORIGIN="${1:-$(cat "$PEGRNA_ROOT/tunnel-url.txt" 2>/dev/null || true)}"
ORIGIN="$(printf '%s' "$ORIGIN" | tr -d '[:space:]')"

if [ -z "$ORIGIN" ]; then
  echo "error: no origin given and $PEGRNA_ROOT/tunnel-url.txt is empty" >&2
  exit 1
fi
if [ ! -s "$CF_TOKEN_FILE" ]; then
  echo "error: no Cloudflare token at $CF_TOKEN_FILE" >&2
  exit 1
fi
if [ ! -f "$SCRIPT_SRC" ]; then
  echo "error: worker source not found at $SCRIPT_SRC" >&2
  exit 1
fi

TOKEN="$(cat "$CF_TOKEN_FILE")"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# Substitute via python so a URL containing sed metacharacters cannot corrupt it.
python3 - "$SCRIPT_SRC" "$TMP/worker.js" "$ORIGIN" <<'PY'
import sys
src, dst, origin = sys.argv[1], sys.argv[2], sys.argv[3]
text = open(src).read().replace("__ORIGIN__", origin)
open(dst, "w").write(text)
PY

resp="$(curl -s --max-time 60 -X PUT \
  -H "Authorization: Bearer $TOKEN" \
  -F "metadata={\"main_module\":\"worker.js\",\"compatibility_date\":\"2024-11-01\"};type=application/json" \
  -F "worker.js=@$TMP/worker.js;type=application/javascript+module" \
  "https://api.cloudflare.com/client/v4/accounts/$CF_ACCOUNT/workers/scripts/$WORKER_NAME")"

if printf '%s' "$resp" | python3 -c "import json,sys; sys.exit(0 if json.load(sys.stdin).get('success') else 1)"; then
  echo "worker '$WORKER_NAME' now proxies to $ORIGIN"
else
  echo "error: worker upload failed:" >&2
  printf '%s\n' "$resp" | python3 -m json.tool 2>/dev/null | head -20 >&2
  exit 1
fi
