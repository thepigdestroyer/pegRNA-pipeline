# Production deployment

What is actually deployed, how to operate it, and the one step that still needs
a human.

```
Browser
  └─ https://pegrna-designer.vercel.app        Vercel: static React app only
       └─ direct XHR (CORS)                    NOT proxied — see "Why direct"
            └─ Cloudflare Tunnel               outbound only; no inbound port
                 └─ DSI 127.0.0.1:8093         FastAPI + the whole model stack
```

## Why direct, not proxied

A Vercel function is capped at **300 s** on this plan. A real `api` design with
live Boltz-2 docking takes **~15 minutes**, so a proxied design would always
time out. The browser therefore calls the backend directly, which makes it
cross-origin — hence `PEGRNA_CORS_ORIGINS` on the backend.

## Where things live (DSI box, as `mj1824`)

`/mnt/data2/mj1824` **does not exist** — Hiru's note pointed there, but all
three `/mnt` volumes are root-owned. Everything is under `$HOME`:

| Path | What |
|---|---|
| `~/pegrna/prime-editing` | the repo (pushed by rsync, no git credentials on the box) |
| `~/pegrna/envs/` | `pegrna-web`, `pegrna-real`, `pridict2` (~6.5 GB) |
| `~/pegrna/PRIDICT2/` | PRIDICT2 + trained weights |
| `~/pegrna/genomes/hg38.2bit` | off-target genome (835 MB) |
| `~/pegrna/runs/` | per-run outputs, pruned to `PEGRNA_KEEP_RUNS` |
| `~/.config/pegrna/pegrna.env` | all configuration: mode 600 (the directory 700; the audit checks) |
| `~/.config/pegrna/nvidia.key` | Boltz-2 key, mode 600 |

**Port 8093**, because `8092` is `jmoore3`'s separate instance. Do not take 8092.

## Services

Both are `systemd --user` units with lingering enabled, so they survive logout
and reboot.

```bash
systemctl --user status pegrna-web pegrna-tunnel
systemctl --user restart pegrna-web
tail -f ~/pegrna/cloudflared.log
curl -s localhost:8093/api/health | python3 -m json.tool
```

## Access control

There is none, by design: the designer is open to everyone, in both modes, with no token.
`api` mode spends NVIDIA credits and occupies the box for ~15 minutes a run, so the cost is
limited by caps rather than a gate:

- `PEGRNA_MAX_CONCURRENCY` (default 2): a further design gets a 429 until a slot frees.
- `PEGRNA_TIMEOUT` (1800 s) ends a run that goes on too long.
- `PEGRNA_MAX_CONFIG_BYTES` and `PEGRNA_MAX_GENOME_LEN` refuse oversized input.
- The thread and priority limits under "Resource limits" below.

If the site is ever abused, lower `PEGRNA_MAX_CONCURRENCY` and restart `pegrna-web`, or stop
it. `scripts/audit-production.sh` checks that nothing is gating the endpoints and that no
retired token is left in `pegrna.env`.

## The tunnel: named (permanent) or quick (fallback)

`tunnel-run.sh` is what the `pegrna-tunnel` unit runs, and it publishes the live URL to
`~/pegrna/tunnel-url.txt` so the audit and deploy scripts never have to guess.

**Named tunnel (permanent hostname).** If `~/.config/pegrna/tunnel.token` exists and
`PEGRNA_PUBLIC_URL` is set, the script does nothing but run that tunnel with
`cloudflared tunnel run --token-file`. One hostname, no rate-limited provisioning, no failover,
and `cloudflared` reconnects on its own after a network blip. It is outbound-only like the rest:
nothing listens on the box.

**Quick tunnel (the fallback, and what runs until the named one is set up).** Providers are tried
in order:

1. **cloudflared quick tunnel** — preferred, but Cloudflare rate-limits
   account-less provisioning (HTTP 429). Once tripped it stays tripped for a
   while, which would otherwise leave the site dark.
2. **localhost.run** (`ssh -R`) — no account, no extra software; used only when
   cloudflared cannot provision.

Both give a random hostname, so the failover keeps the site *up* but does not
make the URL *stable*. Deleting `tunnel.token` returns to this mode.

## Setting up the named tunnel

The quick tunnel is not production-grade, and Cloudflare says so itself: the hostname is random
and changes on every restart, there is no uptime guarantee, and provisioning is rate-limited
(restart it a few times in a row and you get HTTP 429 and no tunnel until it clears). Because the
frontend bakes the backend URL in at build time, a new hostname breaks the site until it is
rebuilt (`scripts/refresh-backend.sh` does that in one command, but it is a workaround).

**Status: done.** The site's backend is `https://pegrna.ukdri.org` (tunnel `pegrna-dsi`, set up on
2026-11-04), and `scripts/audit-production.sh` checks for it. The rest of this section is how it was
done and how to redo it, for another hostname or after a Cloudflare account change.

**Route A: one command (recommended).** Create a Cloudflare API token with
`Account > Cloudflare Tunnel > Edit` and `Zone > DNS > Edit` for the zone that holds the hostname
(dash.cloudflare.com > My Profile > API Tokens > Create Token > Create Custom Token), then store it
without pasting it anywhere:

```bash
read -rs T; security add-generic-password -a "$USER" -s cloudflare-tunnel-token -w "$T" -U; unset T
```

Then, from the repo (the DSI copy must be current, so rsync first):

```bash
scripts/setup-named-tunnel.py --dry-run     # shows the plan; reads no secret, changes nothing
scripts/setup-named-tunnel.py --deploy      # creates it, installs it on DSI, redeploys the site
```

Every Cloudflare step is idempotent (a second run creates and rewrites nothing that is already
right), but each run restarts the tunnel on DSI. The script creates or reuses the tunnel
`pegrna-dsi`, sets its ingress to `http://127.0.0.1:8093` (any other hostname gets a 404), creates
the proxied CNAME, installs the connector token on DSI over ssh stdin (mode 600; never on a command
line or in a log), sets `PEGRNA_PUBLIC_URL`, restarts the service and waits for
`/api/health`. It refuses to overwrite a DNS record of another type and never deletes anything.
Hostname defaults to `pegrna.ukdri.org`; pass `--hostname` to change it.

What the first real run showed:

- The account is taken from the zone's owner, so a token scoped to just a tunnel and a zone (which
  cannot list accounts) works. Pass `--account-id` only if the zone's account cannot be read.
- The hostname takes a few minutes to answer the first time (the connector registers and
  Cloudflare's edge learns it); the script waits five. If it still times out, everything is
  configured: check `ssh dsi tail ~/pegrna/cloudflared.log`, then run the deploy command it prints.
  Do not re-run the script just to wait, since each run restarts the tunnel and the site is
  unreachable until it is back.
- Until the Vercel build is redeployed against the new hostname, the site shows its backend-down
  banner. `--deploy` chains that redeploy; without it, run
  `PEGRNA_API_BASE=https://pegrna.ukdri.org bash web/vercel/deploy.sh --prod`.

**Route B: no API token.** In the Cloudflare dashboard (Zero Trust > Networks > Tunnels) create a
tunnel, add a public hostname `pegrna.<your-domain>` pointing at `http://localhost:8093`, and copy
the connector token it shows. On the DSI box:

```bash
umask 077; mkdir -p ~/.config/pegrna
printf '%s\n' '<connector-token>' > ~/.config/pegrna/tunnel.token
echo 'PEGRNA_PUBLIC_URL=https://pegrna.<your-domain>' >> ~/.config/pegrna/pegrna.env
systemctl --user restart pegrna-tunnel
```

then `PEGRNA_API_BASE=https://pegrna.<your-domain> bash web/vercel/deploy.sh --prod`.

After either route the URL never changes and `refresh-backend.sh` is never needed.

## Resource limits (this host is shared)

The box is a 64-core server other people use. A design is not allowed to take it:

| Setting | Default | What it limits |
|---|---|---|
| `PEGRNA_THREADS` | 4 | TensorFlow / BLAS / ViennaRNA thread pools in the pipeline |
| `PEGRNA_NICE` | 10 | scheduling priority of the whole pipeline (0 disables); other people's jobs win contention |
| `PRIDICT2_THREADS` | 2 | threads per PRIDICT2 worker process (`PRIDICT2_CORES` of them) |
| `PEGRNA_CASOFFINDER_THREADS` | 16 | cas-offinder's CPU backend (`POCL_MAX_PTHREAD_COUNT`) |
| `PEGRNA_MAX_CONCURRENCY` | 2 | simultaneous designs |

Why: uncapped, one design ran a PRIDICT2 worker at ~36 cores and a single off-target scan used
69 CPU-minutes; the load average hit 63. Capped, the same design runs at nice 10 with
cas-offinder at ~5 cores and a load average of ~3.5. `--cores` on PRIDICT2 limits its
*processes*, not the threads inside each, which is why the thread variables exist.

**Off-target results are cached** per guide under `PEGRNA_CACHE_DIR` (keyed on the genome file,
PAM and mismatch settings; replacing the genome invalidates it). The first design for a new set
of guides pays for a genome scan; repeats do not. After deploying, run the example once to warm
it, so the first visitor to try the demo does not pay:

```bash
bash scripts/warm-cache.sh
```

## Routine operations

**Ship a code change**

```bash
rsync -az --delete --exclude .venv --exclude node_modules --exclude web/runs \
  ./ dsi:/home/mj1824/pegrna/prime-editing/
ssh dsi 'systemctl --user restart pegrna-web'
```

**Ship a frontend change** (needs Node locally; the box has none)

```bash
cd web/frontend && npm run build          # commit web/static/ too
PEGRNA_API_BASE=https://<backend> bash web/vercel/deploy.sh --prod
```

**The site loads but cannot design** — the tunnel hostname changed:

```bash
scripts/refresh-backend.sh --prod
```

## Verifying it is really production

The meaningful check is that no gate is a placeholder:

```bash
ssh dsi 'set -a; . ~/.config/pegrna/pegrna.env; set +a
  export NVIDIA_API_KEY=$(cat ~/.config/pegrna/nvidia.key)
  export PATH=$HOME/pegrna/envs/pegrna-real/bin:$PATH
  ~/pegrna/envs/pegrna-real/bin/python -c "
import os
from pegrna_pipeline.hpc.scorer_factory import make_scorers
s = make_scorers(mode=\"api\", uniprot_acc=\"P49711\",
                 genome=os.environ[\"PEGRNA_GENOME\"],
                 nvidia_api_key=os.environ[\"NVIDIA_API_KEY\"])
print(\"stub_gates:\", s.stub_gates or \"() — all real\")"'
```

Expected: `stub_gates: () — all real`.

Anything listed there returned a **placeholder, not a verdict**. The report says
so in a banner and the UI shows a callout. Never order DNA off a stand-in.

## Known behaviour that is not a bug

- **`api` designs can still end in "no design".** Each gate is a real verdict, and a
  design that is genuinely poor (a guide with dozens of near-matches in the genome, an
  extension that folds over the spacer, one that sequesters the Cas9 scaffold) is
  rejected. The result is reproducible: the same inputs give the same verdict every
  time, because Boltz-2 returns identical scores for identical requests.
- **Docking is judged against a neutral RNA of the same length, not against 0.6.**
  Interface ipTM falls with RNA length regardless of design quality, so the usual
  AlphaFold-Multimer floor of 0.6 is out of reach for any tjPE-sized pegRNA. The gate
  averages 5 diffusion samples and asks whether the design looks like a neutral
  extension of its own length or like an RNA that has lost Cas9. See
  `gate_thresholds.DOCKING_REFERENCE` and `scripts/analysis/`.
- **Docking only judges pegRNAs up to 232 nt; longer designs get an advisory.** Past that
  the model's score for a normal RNA and for a non-binder converge until the gap is smaller
  than its own noise (measured to 336 nt), so the design is reported against references of
  its own length and a "Check these before you order" line says it was not judged. How much
  cargo fits depends on the second nick: up to about 60 bp stays inside with the nearest nick,
  and the farther nicks the 40-90 bp window allows cut that to about 25 bp (measured on the
  CTCF locus; a FLAG tag, 24 bp, fits almost everywhere). The range is derived from the data,
  not typed (`scripts/analysis/docking_reference.py`).
- **A design takes ~15 minutes.** cas-offinder alone is ~100 s (batched once per
  design, not per candidate); the rest is Boltz-2.
- **PRIDICT2 can pick a different spacer here than on a laptop.** The box pins
  torch 2.0.1, the version PRIDICT2 upstream supports; a newer torch gives
  slightly different outputs. The box is the reference.
