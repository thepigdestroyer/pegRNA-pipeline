#!/usr/bin/env bash
# Expose the loopback web app publicly, without opening an inbound port.
#
# DSI's guidance was explicit that ports must not be exposed on a shared compute
# box, so every provider here works the same way: an OUTBOUND connection to a
# relay, which then serves the public hostname. Nothing listens externally on
# this machine.
#
# Providers are tried in order and the first that yields a URL wins:
#
#   1. cloudflared quick tunnel  — best behaved, but Cloudflare RATE-LIMITS
#                                  account-less provisioning (HTTP 429). Once
#                                  tripped it can stay tripped for a long while,
#                                  which would otherwise leave the site dark.
#   2. localhost.run (ssh -R)    — no account, no extra software; used only when
#                                  cloudflared cannot provision.
#
# Both give a RANDOM hostname that changes whenever the tunnel restarts. That is
# the one genuinely non-production property of this setup; the permanent fix is
# a named Cloudflare tunnel on a domain you own (see GO-LIVE.md). The live URL is
# always written to $PEGRNA_ROOT/tunnel-url.txt so tooling never has to guess.
set -uo pipefail

PEGRNA_ROOT="${PEGRNA_ROOT:-$HOME/pegrna}"
PEGRNA_PORT="${PEGRNA_PORT:-8093}"
URL_FILE="$PEGRNA_ROOT/tunnel-url.txt"
LOG_FILE="$PEGRNA_ROOT/cloudflared.log"
BIN="$PEGRNA_ROOT/bin/cloudflared"
MAX_LOG_BYTES=$((5 * 1024 * 1024))

# Own the log rather than using StandardOutput=journal: the user journal is not
# available on every host (it is not here, and the unit then fails to start at
# all). Truncate on start so it cannot grow without bound.
# NB: written as `if`, not `[ cond ] && action` — under `set -e` the latter
# exits whenever the condition is false.
if [ -f "$LOG_FILE" ]; then
  size="$(stat -c %s "$LOG_FILE" 2>/dev/null || echo 0)"
  if [ "$size" -gt "$MAX_LOG_BYTES" ]; then : > "$LOG_FILE"; fi
fi

log() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG_FILE"; }
publish() {
  printf '%s' "$1" > "$URL_FILE"
  log "pegrna: tunnel URL -> $1"
}

: > "$URL_FILE"

# --- provider 1: cloudflared -------------------------------------------------
try_cloudflared() {
  [ -x "$BIN" ] || { log "cloudflared: binary missing at $BIN"; return 1; }
  log "cloudflared: requesting a quick tunnel"

  stdbuf -oL -eL "$BIN" tunnel --no-autoupdate --url "http://127.0.0.1:$PEGRNA_PORT" 2>&1 |
    while IFS= read -r line; do
      printf '%s\n' "$line" >> "$LOG_FILE"
      case "$line" in
        *trycloudflare.com*)
          url="$(printf '%s' "$line" | grep -oE 'https://[a-z0-9-]+\.trycloudflare\.com' | head -1)"
          if [ -n "$url" ] && [ "$url" != "$(cat "$URL_FILE" 2>/dev/null)" ]; then
            publish "$url"
          fi
          ;;
      esac
    done
  return 0
}

# While on the fallback, keep checking whether Cloudflare will serve us again.
# The fallback rotates its hostname roughly every 20 minutes, which breaks the
# site each time; a cloudflared tunnel keeps one hostname for the life of the
# process. So as soon as the rate limit lifts we want to be back on it. On
# success this exits the whole script so systemd restarts us into provider 1.
watch_for_cloudflared() {
  local probe
  while sleep 900; do
    [ -x "$BIN" ] || continue
    probe="$(mktemp)"
    timeout 45 "$BIN" tunnel --no-autoupdate --url "http://127.0.0.1:$PEGRNA_PORT" \
      > "$probe" 2>&1 || true
    if grep -q 'trycloudflare.com' "$probe" && ! grep -q '429' "$probe"; then
      log "cloudflared: available again — restarting onto it for a stable hostname"
      rm -f "$probe"
      kill -TERM -$$ 2>/dev/null
      exit 0
    fi
    rm -f "$probe"
    log "cloudflared: still unavailable; staying on the fallback"
  done
}

# --- provider 2: localhost.run ----------------------------------------------
try_localhost_run() {
  log "localhost.run: opening reverse tunnel"
  ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=~/.ssh/known_hosts \
      -o ServerAliveInterval=30 -o ServerAliveCountMax=3 -o ExitOnForwardFailure=yes \
      -R "80:127.0.0.1:$PEGRNA_PORT" nokey@localhost.run 2>&1 |
    while IFS= read -r line; do
      printf '%s\n' "$line" >> "$LOG_FILE"
      case "$line" in
        *lhr.life*|*localhost.run*)
          url="$(printf '%s' "$line" | grep -oE 'https://[a-z0-9-]+\.lhr\.life' | head -1)"
          if [ -n "$url" ] && [ "$url" != "$(cat "$URL_FILE" 2>/dev/null)" ]; then
            publish "$url"
          fi
          ;;
      esac
    done
  return 0
}

# A provider "succeeds" only if it actually published a URL; cloudflared exits
# almost immediately on a 429, so a short probe distinguishes the two cases.
run_provider() {
  local fn="$1" name="$2" waited=0
  "$fn" &
  local pid=$!
  while [ "$waited" -lt 40 ]; do
    if [ -s "$URL_FILE" ]; then
      log "$name: up at $(cat "$URL_FILE")"
      wait "$pid"
      return 0
    fi
    if ! kill -0 "$pid" 2>/dev/null; then
      log "$name: exited without providing a URL"
      return 1
    fi
    sleep 2
    waited=$((waited + 2))
  done
  log "$name: no URL after ${waited}s — giving up on this provider"
  kill "$pid" 2>/dev/null
  return 1
}

run_provider try_cloudflared "cloudflared" && exit 0

log "falling back to localhost.run"
watch_for_cloudflared &
WATCHDOG=$!
trap 'kill "$WATCHDOG" 2>/dev/null' EXIT
run_provider try_localhost_run "localhost.run" && exit 0

log "no tunnel provider succeeded"
exit 1
