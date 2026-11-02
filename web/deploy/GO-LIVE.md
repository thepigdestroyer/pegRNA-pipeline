# Go live on the DSI box — production runbook

Target: `146-169-8-141.dsi.ic.ac.uk` as `mj1824`, storage `/mnt/data2/mj1824`.

Every step below is idempotent — safe to re-run.

---

## 0. Prerequisites you have to do yourself

| Thing | Why it can't be automated |
|---|---|
| **Be on the Imperial network** | `146.169.8.141:22` is filtered from outside. Imperial retired its central SSH gateway and the old VPNs; access is via **Unified Access**, which needs your Microsoft login + MFA. |
| **NVIDIA API key** | Personal credential tied to your NVIDIA account (see §2). |

Log in with your **Imperial short-form username and Microsoft account password** —
per Hiru, `.141` is on Imperial auth now. Type the password rather than pasting it:

```bash
ssh -vvv mj1824@146.169.8.141
```

If auth still fails, send Hiru the full `-vvv` output.

---

## 1. Deploy

**The repo is private**, so the box needs credentials before it can clone.
Easiest is a GitHub Personal Access Token with `repo` scope
(<https://github.com/settings/tokens>), used once and then cached:

```bash
ssh mj1824@146.169.8.141

export PEGRNA_ROOT=/mnt/data2/mj1824/pegrna
mkdir -p "$PEGRNA_ROOT"

# Clone once, with the token inline. git stores it in the remote URL, so the
# later `git pull` in deploy.sh just works.
read -rsp 'GitHub PAT: ' GH_PAT; echo
git clone "https://${GH_PAT}@github.com/thepigdestroyer/pegRNA-pipeline.git" \
  "$PEGRNA_ROOT/prime-editing"
unset GH_PAT

bash "$PEGRNA_ROOT/prime-editing/web/deploy/deploy.sh"
```

`deploy.sh` then builds `pegrna-web`, `pegrna-real` and `pridict2`, installs the
`systemd --user` unit, starts it on `127.0.0.1:8092`, and health-checks it.

> The frontend ships **pre-built** in `web/static/`, so the box needs no Node.

> If you would rather not leave a token in the remote URL, use a read-only
> **deploy key** instead: `ssh-keygen -t ed25519 -f ~/.ssh/pegrna_deploy -N ''`,
> add the `.pub` under the repo's *Settings → Deploy keys*, and clone the
> `git@github.com:` URL.

---

## 2. NVIDIA key — turns docking from a stand-in into a real model

1. Sign in at **<https://build.nvidia.com>** (free NVIDIA account).
2. Open **<https://build.nvidia.com/mit/boltz2>**.
3. Click **Get API Key**. The value starts with `nvapi-`.
4. Install it on the box, owner-only:

```bash
mkdir -p ~/.config/pegrna
printf '%s' 'nvapi-REPLACE_ME' > ~/.config/pegrna/nvidia.key
chmod 600 ~/.config/pegrna/nvidia.key
systemctl --user restart pegrna-web
```

Boltz-2 calls are heavy — minutes each, and they consume credits. The free pool
is fine for validation; sustained use needs credits on the NVIDIA account.

Never commit it, and never put it in the systemd unit. The app reads
`NVIDIA_API_KEY`, else `~/.config/pegrna/nvidia.key`.

**Check it took:**

```bash
curl -s localhost:8092/api/health | grep -o '"nvidia_key":[a-z]*'   # want: true
```

---

## 3. Genome FASTA — turns off-target into a real model

```bash
PEGRNA_ROOT=/mnt/data2/mj1824/pegrna \
  bash $PEGRNA_ROOT/prime-editing/scripts/install_genome.sh
systemctl --user restart pegrna-web
```

~900 MB download, ~3.1 GB unpacked, to `/mnt/data2` (not the home quota).
`PEGRNA_GENOME` is already set in `pegrna.env` to the path this writes.

---

## 4. Verify there are no stand-ins left

This is the acceptance test for "production". Run a real design and confirm
`stand_in_gates` is **empty**:

```bash
curl -s -X POST localhost:8092/api/design \
  -H 'Content-Type: application/json' \
  -d "{\"mode\":\"api\",\"config\":$(cat $PEGRNA_ROOT/prime-editing/examples/ctcf_cterm/config_autospacer_flag.json)}" \
  | python3 -c 'import json,sys; r=json.load(sys.stdin); print("ok:", r["ok"]); print("stand_in_gates:", r["gates"].get("stand_in_gates"))'
```

Expected on a complete production deployment:

```
ok: True
stand_in_gates: []
```

Anything listed there is a check that returned a **placeholder, not a verdict**.
The report says so in a banner and the UI shows a callout — do not let a design
go to synthesis off a stand-in.

---

## 5. Publish — without opening a port

Hiru's objection was explicit: *"it's not a good idea to expose ports on machines
which aren't dedicated towards hosting."* Tailscale Funnel satisfies that. The
app binds **loopback only**; the tunnel is an **outbound** connection to
Tailscale's relays, so **no inbound port is opened** and no root is needed
(userspace networking).

```bash
bash $PEGRNA_ROOT/prime-editing/web/deploy/funnel.sh up
bash $PEGRNA_ROOT/prime-editing/web/deploy/funnel.sh status   # prints the public URL
```

Needs, one time in the Tailscale admin console: HTTPS + MagicDNS enabled, and
the `funnel` node attribute granted to this node.

To take it down: `funnel.sh down`.

### If it should not be public at all

Skip §5 and use an SSH tunnel per user — nothing is exposed:

```bash
ssh -L 8092:127.0.0.1:8092 mj1824@146.169.8.141
# then browse http://localhost:8092
```

---

## 6. Lock it down before sharing the URL

A Funnel URL is public to anyone who has it. At minimum set an access token:

```bash
# in ~/.config/pegrna/pegrna.env
PEGRNA_ACCESS_TOKEN=<a long random string>
```

```bash
systemctl --user restart pegrna-web
```

Clients then need `X-Access-Token` (or `?token=`). Also worth reviewing in
`pegrna.env`: `PEGRNA_MAX_CONCURRENCY` (default 2 — this also protects your
NVIDIA credit pool), `PEGRNA_TIMEOUT`, `PEGRNA_KEEP_RUNS`.

---

## Operating it

```bash
systemctl --user status pegrna-web
systemctl --user restart pegrna-web
journalctl --user -u pegrna-web -f
curl -s localhost:8092/api/health | python3 -m json.tool
```

`loginctl enable-linger mj1824` is set by `deploy.sh` so the service survives
logout.

### Shipping a change

```bash
cd $PEGRNA_ROOT/prime-editing && git pull --ff-only
systemctl --user restart pegrna-web
```

If the change touched `web/frontend/`, rebuild it **before** committing (needs
Node, so do it on your laptop, not here):

```bash
cd web/frontend && npm ci && npm run build
git add ../static && git commit
```
