#!/usr/bin/env bash
# Production audit for the live deployment.
#
# Asserts the things that must be true and names any that are not. Run it after
# any deployment change, and whenever the site misbehaves -- most outages here
# are one of these checks failing (usually the tunnel hostname rotating).
#
#   scripts/audit-production.sh
#
# Requires an ssh alias `dsi` for the backend host.
SITE="https://pegrna-designer.vercel.app"
pass=0; fail=0
ok()   { printf '  \033[32mPASS\033[0m  %s\n' "$1"; pass=$((pass+1)); }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }
chk()  { if eval "$2" >/dev/null 2>&1; then ok "$1"; else bad "$1"; fi; }

URL=$(ssh dsi 'cat ~/pegrna/tunnel-url.txt' 2>/dev/null | tr -d '[:space:]')

echo "BACKEND (DSI)"
chk "web service active"            "ssh dsi 'systemctl --user is-active pegrna-web | grep -q active'"
chk "tunnel service active"         "ssh dsi 'systemctl --user is-active pegrna-tunnel | grep -q active'"
chk "both enabled at boot"          "ssh dsi 'systemctl --user is-enabled pegrna-web pegrna-tunnel | grep -qc enabled'"
chk "lingering on (survives logout)" "ssh dsi 'loginctl show-user mj1824 | grep -q Linger=yes'"
chk "health endpoint responds"      "ssh dsi 'curl -sf --max-time 10 localhost:8093/api/health'"
chk "NVIDIA key present"            "ssh dsi 'curl -s localhost:8093/api/health | grep -q \"\\\"nvidia_key\\\":true\"'"
chk "api python resolves"           "ssh dsi 'curl -s localhost:8093/api/health | grep -q \"\\\"api_python_ok\\\":true\"'"
chk "genome staged for off-target"  "ssh dsi 'test -s ~/pegrna/genomes/hg38.2bit'"
chk "cas-offinder installed"        "ssh dsi 'test -x ~/pegrna/envs/pegrna-real/bin/cas-offinder'"
chk "PRIDICT2 env present"          "ssh dsi 'test -x ~/pegrna/envs/pridict2/bin/python'"
# The host is shared and its home directories are group-readable. pegrna.env holds the api-mode
# token (it gates paid model runs), so it is a secret like the NVIDIA key and the tunnel token.
chk "secrets private (config dir 700, every file in it 600)" "ssh dsi 'test \$(stat -c %a ~/.config/pegrna) = 700 && ! find ~/.config/pegrna -type f ! -perm 600 | grep -q .'"
chk "jmoore3 instance untouched"    "ssh dsi 'curl -s --max-time 10 localhost:8092/api/health | grep -q jmoore3'"
chk "disk under 90%"                "ssh dsi '[ \$(df --output=pcent / | tail -1 | tr -dc 0-9) -lt 90 ]'"
# The host is shared; a design must not be able to take it (see GO-LIVE.md, "Resource limits").
chk "resource caps configured"      "ssh dsi 'grep -q ^PEGRNA_THREADS= ~/.config/pegrna/pegrna.env && grep -q ^PEGRNA_NICE= ~/.config/pegrna/pegrna.env && grep -q ^PRIDICT2_THREADS= ~/.config/pegrna/pegrna.env && grep -q ^PEGRNA_CASOFFINDER_THREADS= ~/.config/pegrna/pegrna.env'"
chk "off-target cache enabled"      "ssh dsi 'grep -q ^PEGRNA_CACHE_DIR= ~/.config/pegrna/pegrna.env'"
chk "host not saturated (load < 48)" "ssh dsi '[ \$(cut -d. -f1 /proc/loadavg) -lt 48 ]'"

# A quick tunnel works but its hostname is random and changes on restart, so it is reported, not
# failed: the fix needs a Cloudflare token (GO-LIVE.md, "Setting up the named tunnel").
if ssh dsi 'test -s ~/.config/pegrna/tunnel.token' 2>/dev/null; then
  ok "named tunnel: permanent hostname ($URL)"
else
  printf '  \033[33mNOTE\033[0m  quick tunnel: the hostname is random and changes on restart (scripts/setup-named-tunnel.py)\n'
fi

echo
echo "ZERO STAND-IN GATES"
STUBS=$(ssh dsi 'set -a; . ~/.config/pegrna/pegrna.env; set +a
  export NVIDIA_API_KEY=$(cat ~/.config/pegrna/nvidia.key)
  ~/pegrna/envs/pegrna-real/bin/python -c "
import os
from pegrna_pipeline.hpc.scorer_factory import make_scorers
s=make_scorers(mode=\"api\", uniprot_acc=\"P49711\", genome=os.environ.get(\"PEGRNA_GENOME\"), nvidia_api_key=os.environ.get(\"NVIDIA_API_KEY\"))
print(len(s.stub_gates))"' 2>/dev/null | tr -d '[:space:]')
if [ "$STUBS" = "0" ]; then ok "every gate is a real model (stub_gates empty)"; else bad "stub_gates has $STUBS entries"; fi

echo
echo "ACCESS CONTROL"
chk "api mode requires a token"     "[ \$(curl -s -o /dev/null -w '%{http_code}' --max-time 30 -X POST $URL/api/design -H 'Content-Type: application/json' -d '{\"mode\":\"api\",\"config\":{}}') = 401 ]"
chk "fake mode stays open"          "[ \$(curl -s -o /dev/null -w '%{http_code}' --max-time 60 -X POST $URL/api/design -H 'Content-Type: application/json' -d '{\"mode\":\"fake\",\"config\":{\"genome\":\"ACGT\",\"site_s\":2,\"cargo\":\"ACG\",\"rtt_aa\":\"L\",\"cds_mask\":[true,true,false,false]}}') = 200 ]"

echo
echo "PUBLIC SITE"
chk "site returns 200"              "[ \$(curl -s -o /dev/null -w '%{http_code}' --max-time 30 $SITE) = 200 ]"
chk "JS bundle resolves (not SPA fallback)" "curl -s $SITE | grep -oE '/assets/index-[A-Za-z0-9_-]+\.js' | head -1 | xargs -I{} curl -sf --max-time 20 -o /dev/null $SITE{}"
chk "missing asset 404s properly"   "[ \$(curl -s -o /dev/null -w '%{http_code}' --max-time 20 $SITE/assets/nope.js) = 404 ]"
chk "site build matches live tunnel" "curl -s \$(curl -s $SITE | grep -oE '/assets/index-[A-Za-z0-9_-]+\.js' | head -1 | sed 's|^|$SITE|') | grep -q '$URL'"
chk "CORS allows the site origin"   "curl -s -I --max-time 20 -H 'Origin: $SITE' $URL/api/health | grep -qi 'access-control-allow-origin'"

echo
printf 'RESULT: %d passed, %d failed\n' "$pass" "$fail"
[ "$fail" -eq 0 ] || exit 1
