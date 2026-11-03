#!/usr/bin/env bash
# Publish the pegRNA UI to Vercel.
#
# Vercel hosts the frontend ONLY. The pipeline cannot run there: it needs
# PRIDICT2 (torch), SpliceAI (TensorFlow), ViennaRNA, cas-offinder and an
# 835 MB genome, and a single design takes ~15 minutes.
#
# The browser calls the backend DIRECTLY rather than through a Vercel proxy:
# serverless functions are capped at 300 s on this plan and a real api-mode
# design runs far longer, so a proxied design would always time out.
#
# Two non-obvious things this handles, both of which cost an afternoon:
#
#   * Vercel reads git metadata from the working directory and attributes the
#     deployment to the last commit's author. If that author lacks deploy rights
#     on the project, the build is Blocked with a message about commit authors
#     that has nothing to do with the code. Staging outside the repo makes it a
#     plain upload with no git metadata at all.
#
#   * Vite's base must be "/" here. The default build targets FastAPI's /static
#     mount; with the wrong base the bundle 404s and the page renders blank --
#     and because the SPA fallback used to serve index.html for missing assets,
#     it did so with HTTP 200 and no console error.
#
# Usage:
#   PEGRNA_API_BASE=https://host ./deploy.sh [--prod]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRONTEND="$HERE/../frontend"

: "${PEGRNA_API_BASE:?set PEGRNA_API_BASE to the public URL of the pipeline backend}"
PEGRNA_API_BASE="${PEGRNA_API_BASE%/}"

echo "==> backend: $PEGRNA_API_BASE"
if curl -sf --max-time 25 "$PEGRNA_API_BASE/api/health" >/dev/null; then
  echo "    health check OK"
else
  echo "    WARNING: $PEGRNA_API_BASE/api/health did not respond." >&2
  echo "             Deploying anyway, but the site will not be able to design." >&2
fi

# Stage OUTSIDE the repository (see the git-metadata note above).
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
mkdir -p "$STAGE/public"

echo "==> building the frontend against that backend"
cd "$FRONTEND"
[ -d node_modules ] || npm ci
# VITE_BASE=/ because Vercel serves these files from the root, whereas the default build
# targets FastAPI's /static mount. VITE_OUT_DIR sends the output straight to the staging
# directory: building into ../static would overwrite the COMMITTED FastAPI-served bundle with
# one that has this deployment's backend URL baked in (which has happened).
VITE_OUT_DIR="$STAGE/public" VITE_BASE=/ VITE_API_BASE="$PEGRNA_API_BASE" \
  sh -c 'npx tsc -b && npx vite build'

echo "==> staging into $STAGE"
cp "$HERE/vercel.json" "$STAGE/"
# Carry the project link across so this deploys to the same Vercel project.
[ -d "$HERE/.vercel" ] && cp -R "$HERE/.vercel" "$STAGE/"

echo "==> deploying"
cd "$STAGE"
if [ "${1:-}" = "--prod" ]; then
  vercel deploy --prod --yes
else
  vercel deploy --yes
fi

# Persist the project link created on first deploy.
[ -d "$STAGE/.vercel" ] && cp -R "$STAGE/.vercel" "$HERE/"

cat <<EOF

Deployed against: $PEGRNA_API_BASE

The backend must allow this site's origin, or the browser will block every call:

    # ~/.config/pegrna/pegrna.env on the backend host
    PEGRNA_CORS_ORIGINS=https://pegrna-designer.vercel.app
    systemctl --user restart pegrna-web

The backend URL is baked in at BUILD time. If the tunnel hostname changes, run
scripts/refresh-backend.sh (which does all of this for you).
EOF
