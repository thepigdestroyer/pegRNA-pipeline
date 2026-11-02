#!/usr/bin/env bash
# Publish the pegRNA UI to Vercel.
#
# Vercel hosts the frontend ONLY. The pipeline cannot run there: it needs
# PRIDICT2 (torch), SpliceAI (TensorFlow), ViennaRNA, cas-offinder and a ~3 GB
# genome FASTA, and a single design takes minutes.
#
# The browser therefore calls the pipeline backend DIRECTLY, not through a
# Vercel proxy: a serverless function is capped at 300 s on this plan, and a real
# api-mode design with live Boltz-2 docking runs longer than that, so a proxied
# design would always time out. Direct calls have no such ceiling.
#
# That makes it a cross-origin request, so the backend must allow this origin:
#   PEGRNA_CORS_ORIGINS=https://<your-vercel-domain>   (in ~/.config/pegrna/pegrna.env)
#
# Required:
#   PEGRNA_API_BASE   public URL of the pipeline backend (the tunnel hostname)
#
#   PEGRNA_API_BASE=https://xxxx.trycloudflare.com ./deploy.sh --prod
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRONTEND="$HERE/../frontend"

: "${PEGRNA_API_BASE:?set PEGRNA_API_BASE to the public URL of the pipeline backend}"

echo "==> backend: $PEGRNA_API_BASE"
if ! curl -sf --max-time 20 "$PEGRNA_API_BASE/api/health" >/dev/null; then
  echo "WARNING: $PEGRNA_API_BASE/api/health is not responding." >&2
  echo "         Deploying anyway, but the site will not be able to design." >&2
fi

echo "==> building the frontend against that backend"
cd "$FRONTEND"
[ -d node_modules ] || npm ci
# VITE_BASE=/ because Vercel serves these files from the root, whereas the
# default build targets FastAPI's /static mount.
VITE_BASE=/ VITE_API_BASE="$PEGRNA_API_BASE" npm run build

echo "==> staging into public/"
rm -rf "$HERE/public"
mkdir -p "$HERE/public"
cp -R "$HERE/../static/." "$HERE/public/"

echo "==> deploying to Vercel"
cd "$HERE"
if [ "${1:-}" = "--prod" ]; then
  vercel deploy --prod --yes
else
  vercel deploy --yes
fi

cat <<EOF

Done. Two things have to agree or the site cannot call the backend:

  1. This build points at : $PEGRNA_API_BASE
  2. The backend must allow the deployment's origin:

       # in ~/.config/pegrna/pegrna.env on the backend host
       PEGRNA_CORS_ORIGINS=https://<your-vercel-domain>
       systemctl --user restart pegrna-web

Note the frontend is baked with the backend URL at BUILD time, so if the tunnel
hostname changes you must re-run this script.
EOF
