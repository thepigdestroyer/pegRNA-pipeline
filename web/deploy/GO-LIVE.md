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
| `~/.config/pegrna/pegrna.env` | all configuration |
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

`api` mode costs NVIDIA credits and occupies the box for ~15 minutes, so it is
token-gated. `fake` mode is open to everyone, so the site still works as a
public demo.

```
PEGRNA_API_MODE_TOKEN=<token>     # required for api mode only
PEGRNA_ACCESS_TOKEN=              # unset: would gate the entire site
```

The UI asks for the token only when `/api/health` reports
`api_mode_token_required: true`, and remembers it in `localStorage`.

## Tunnel providers and failover

`tunnel-run.sh` tries providers in order and publishes whichever comes up to
`~/pegrna/tunnel-url.txt`:

1. **cloudflared quick tunnel** — preferred, but Cloudflare rate-limits
   account-less provisioning (HTTP 429). Once tripped it stays tripped for a
   while, which would otherwise leave the site dark.
2. **localhost.run** (`ssh -R`) — no account, no extra software; used only when
   cloudflared cannot provision.

Both give a random hostname, so the failover keeps the site *up* but does not
make the URL *stable*.

## ⚠ The one outstanding step: a named tunnel

The deployment currently uses a **quick tunnel**, which is not production-grade
and Cloudflare says so itself:

- the hostname is **random and changes on every restart**
- there is **no uptime guarantee**
- provisioning is **rate-limited** — restart it a few times in a row and you get
  HTTP 429 and no tunnel at all until it clears

Because the frontend bakes the backend URL in at build time, a new hostname
breaks the site until it is rebuilt. `scripts/refresh-backend.sh` does that in
one command, but it is a workaround.

**The fix needs a Cloudflare login (about five minutes) and then never recurs:**

```bash
# on the DSI box
~/pegrna/bin/cloudflared tunnel login          # opens a browser URL
~/pegrna/bin/cloudflared tunnel create pegrna
~/pegrna/bin/cloudflared tunnel route dns pegrna pegrna.<your-domain>

mkdir -p ~/.cloudflared
cat > ~/.cloudflared/config.yml <<'EOF'
tunnel: pegrna
credentials-file: /home/mj1824/.cloudflared/<TUNNEL-ID>.json
ingress:
  - hostname: pegrna.<your-domain>
    service: http://127.0.0.1:8093
  - service: http_status:404
EOF
```

Then point the unit at `cloudflared tunnel run pegrna`, and deploy the site once
against that permanent hostname:

```bash
PEGRNA_API_BASE=https://pegrna.<your-domain> bash web/vercel/deploy.sh --prod
```

After that the URL never changes and `refresh-backend.sh` is never needed.

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
- **A design takes ~15 minutes.** cas-offinder alone is ~100 s (batched once per
  design, not per candidate); the rest is Boltz-2.
- **PRIDICT2 can pick a different spacer here than on a laptop.** The box pins
  torch 2.0.1, the version PRIDICT2 upstream supports; a newer torch gives
  slightly different outputs. The box is the reference.
