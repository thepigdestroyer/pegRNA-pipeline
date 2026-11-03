#!/usr/bin/env bash
# Run the Cloudflare tunnel and publish its current public URL to a stable file.
#
# Two problems this solves:
#
#   * A quick tunnel's hostname is random and only ever appears once, in the
#     startup log. Writing it to a known path means tooling can find the current
#     URL without grepping an ever-growing log.
#
#   * cloudflared is chatty. Appending to a file forever fills the disk on a
#     shared box eventually, so the verbose stream goes to journald (which
#     systemd rotates) and only the URL is persisted.
#
# Reads PEGRNA_PORT (default 8093) and PEGRNA_ROOT (default ~/pegrna).
set -euo pipefail

PEGRNA_ROOT="${PEGRNA_ROOT:-$HOME/pegrna}"
PEGRNA_PORT="${PEGRNA_PORT:-8093}"
URL_FILE="$PEGRNA_ROOT/tunnel-url.txt"
LOG_FILE="$PEGRNA_ROOT/cloudflared.log"
BIN="$PEGRNA_ROOT/bin/cloudflared"
MAX_LOG_BYTES=$((5 * 1024 * 1024))

# The wrapper owns its own log rather than using StandardOutput=journal: the
# user journal is not guaranteed to be available on a shared host (it is not on
# the DSI box, where it makes the unit fail to start outright). Truncating on
# each start bounds growth without needing logrotate.
# NB: written as `if`, not `[ ... ] && ...`. Under `set -e` the latter exits the
# script whenever the condition is false, which silently turned "log is small"
# into a service that failed to start.
if [ -f "$LOG_FILE" ]; then
  size="$(stat -c %s "$LOG_FILE" 2>/dev/null || echo 0)"
  if [ "$size" -gt "$MAX_LOG_BYTES" ]; then
    : > "$LOG_FILE"
  fi
fi

: > "$URL_FILE"

# Line-buffer so the URL is captured as soon as cloudflared prints it, not when
# some pipe buffer happens to flush.
stdbuf -oL -eL "$BIN" tunnel --no-autoupdate --url "http://127.0.0.1:$PEGRNA_PORT" 2>&1 |
  while IFS= read -r line; do
    printf '%s\n' "$line" >> "$LOG_FILE"
    case "$line" in
      *trycloudflare.com*)
        url="$(printf '%s' "$line" | grep -oE 'https://[a-z0-9-]+\.trycloudflare\.com' | head -1)"
        if [ -n "$url" ] && [ "$url" != "$(cat "$URL_FILE" 2>/dev/null)" ]; then
          printf '%s' "$url" > "$URL_FILE"
          printf 'pegrna: tunnel URL -> %s\n' "$url" >> "$LOG_FILE"
        fi
        ;;
    esac
  done
