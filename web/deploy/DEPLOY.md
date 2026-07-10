# Deploying the pegRNA web app on the DSI server

Target: `146-169-8-141.dsi.ic.ac.uk` as `mj1824`, storage under `/mnt/data2/mj1824`,
public HTTPS via **Tailscale Funnel** (no inbound ports opened — per DSI's security
guidance). No local GPU: docking is the NVIDIA NIM API, everything else is CPU/API.

## 0. Prerequisites (one-time, human)

- **SSH access to the box.** `mj1824`'s Imperial/CX3 password is currently rejected on
  `.141` (that host is not on Active Directory yet, per the DSI email thread). Resolve one of:
  - DSI admin (Hiru/Danie) sets a local password or adds this deploy key to
    `~/.ssh/authorized_keys`, **or**
  - you log in interactively once and append the deploy key yourself.
- **Tailscale**: the node must be allowed to use Funnel in the tailnet ACLs
  (`nodeAttrs` → `funnel`), with HTTPS + MagicDNS enabled for the tailnet.
- **NVIDIA NIM key** for Boltz-2 docking (see `../KEYS.md`). Optional — without it the
  app runs every other gate for real and stubs docking.

## 1. One-command deploy

```bash
ssh mj1824@146-169-8-141.dsi.ic.ac.uk
export PEGRNA_ROOT=/mnt/data2/mj1824/pegrna
curl -Ls https://raw.githubusercontent.com/neurogenomics/prime-editing/main/web/deploy/deploy.sh -o /tmp/deploy.sh
bash /tmp/deploy.sh
```

`deploy.sh` will: clone to `$PEGRNA_ROOT/prime-editing`, build the `pegrna-web` and
`pegrna-real` conda envs (bootstrapping micromamba if needed), install the
`systemd --user` unit + `~/.config/pegrna/pegrna.env`, start the service on
`127.0.0.1:8092`, and health-check it.

## 2. Add the NVIDIA key (enables real docking)

```bash
mkdir -p ~/.config/pegrna
printf %s 'nvapi-XXXXXXXX...' > ~/.config/pegrna/nvidia.key
chmod 600 ~/.config/pegrna/nvidia.key
systemctl --user restart pegrna-web
```

## 3. Go public (Tailscale Funnel — no open ports)

```bash
PEGRNA_ROOT=/mnt/data2/mj1824/pegrna bash "$PEGRNA_ROOT/prime-editing/web/deploy/funnel.sh" up
```

Prints a public `https://pegrna-dsi.<tailnet>.ts.net` URL. That is the operable
website — anyone can open it, fill the form, and run a design.

## 4. Verify

```bash
curl -s http://127.0.0.1:8092/api/health            # local
# then open the Funnel URL; run "Load working example (real · api)" -> Design.
# Watch the Docking node show a real ipTM (~0.6–0.75). "docking_stub" = key missing.
```

## Operations

| Action | Command |
|---|---|
| status | `systemctl --user status pegrna-web` |
| logs | `journalctl --user -u pegrna-web -f` |
| restart | `systemctl --user restart pegrna-web` |
| update code | `git -C $PEGRNA_ROOT/prime-editing pull && systemctl --user restart pegrna-web` |
| funnel status/stop | `web/deploy/funnel.sh status` / `... down` |
| require a token | set `PEGRNA_ACCESS_TOKEN=...` in `~/.config/pegrna/pegrna.env`, restart |

## Config knobs (`~/.config/pegrna/pegrna.env`)

`PEGRNA_REPO`, `PEGRNA_WEB_PYTHON`, `PEGRNA_API_PYTHON`, `PRIDICT2_DIR`,
`PRIDICT2_PYTHON`, `PRIDICT2_CORES`, `PEGRNA_RUNS_DIR`, `PEGRNA_KEEP_RUNS`,
`PEGRNA_PORT`, `NVIDIA_KEY_FILE`, `PEGRNA_MAX_CONCURRENCY`, `PEGRNA_TIMEOUT`,
`PEGRNA_MAX_CONFIG_BYTES`, `PEGRNA_MAX_GENOME_LEN`, `PEGRNA_ACCESS_TOKEN`.
See `pegrna.env` for defaults.

## Why this shape

- **CPU/API only** — `api` mode needs no GPU: SpliceAI (TF-CPU), MaxEntScan (vendored),
  ViennaRNA (CPU), AlphaFold-DB (API), Boltz-2 (NVIDIA NIM API). `setup_envs.sh` installs
  exactly this stack in `pegrna-real` (no torch/boltz there).
- **Auto-spacer (PRIDICT2)** — `setup_envs.sh` also clones PRIDICT2 and builds a separate
  `pridict2` env (python 3.10, CPU `torch==2.0.1`) used for `api`-mode pegRNA design. Its
  location is carried by `pegrna.env` (`PRIDICT2_DIR` / `PRIDICT2_PYTHON` / `PRIDICT2_CORES`),
  which `web/app.py` injects into the `api` subprocess. Torch lives only in this env, so
  `pegrna-real` stays torch-free.
- **No open ports** — uvicorn binds `127.0.0.1`; Funnel tunnels outbound. Nothing to
  firewall on the shared box.
- **Bounded** — concurrency cap + payload/genome caps + run retention keep a public
  endpoint from exhausting the box or the NVIDIA credit pool.
