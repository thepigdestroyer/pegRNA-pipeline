#!/usr/bin/env bash
# Install the vendored MaxEntScan (maxentpy) into the active environment.
#
# maxentpy is the PRIMARY cryptic-splice gate, but it is not pip-installable on
# Python >= 3.11: its pre-generated Cython `_hashseq.c` includes the private
# CPython header `longintrepr.h`. The pipeline only uses the pure-Python
# `maxent.py` scorer (score5/score3), so this installs the vendored copy with the
# C extension stripped out.
#
# Previously this lived only inside web/deploy/setup_envs.sh, which meant api
# mode failed at the cryptic-splice gate on every host except the DSI box.
#
# Usage:  scripts/install_maxentpy.sh          (installs into the active env)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$REPO_ROOT/web/maxentpy_src"
PYTHON="${PYTHON:-python}"

if [ ! -d "$SRC" ]; then
  echo "error: vendored maxentpy source not found at $SRC" >&2
  exit 1
fi

if "$PYTHON" -c "import maxentpy.maxent" >/dev/null 2>&1; then
  echo "maxentpy already importable — nothing to do."
  exit 0
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp -r "$SRC/." "$TMP/"

# Strip the C extension from setup.py (portable across GNU and BSD sed).
"$PYTHON" - "$TMP/setup.py" <<'PY'
import re, sys
path = sys.argv[1]
with open(path) as fh:
    text = fh.read()
with open(path, "w") as fh:
    fh.write(re.sub(r"ext_modules\s*=\s*ext_modules\s*,", "", text))
PY

"$PYTHON" -m pip install --no-build-isolation "$TMP"
"$PYTHON" -c "import maxentpy.maxent; print('maxentpy installed OK (pure-Python)')"
