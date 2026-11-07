# pegRNA design — web frontend

A FastAPI web app that wraps the `pegrna design-recoded` (tj-PE) and `pegrna design-pa`
(Prime Assembly) CLIs (this repo). The browser posts a design config; the server
writes it to a per-run directory, runs the pipeline as a subprocess, streams
per-gate progress back, and renders the text report + parsed CSV components.

## What it designs

You supply the **target** (genome window, insert `site_s`, `rtt_aa`, strand) and
the **cargo** (insert). The pipeline designs the guides — you do **not** type the
spacer:

- **PRIDICT2** (local, CPU) designs the pegRNA spacer plus PBS/RTT at the insert
  site. The spacer is auto-designed by default; a **spacer override** is exposed
  under *Advanced overrides* in the form (and as an optional `spacer` config key)
  for the rare case a fixed guide is required.
- The pipeline owns *where* the ngRNA nicks — it enumerates and gates candidate
  nick positions, then tells PRIDICT2 to design pegRNA#2 there (whose
  spacer+scaffold becomes the ngRNA).
- Each candidate then runs the **gate cascade** (same names shown in the form's
  flowchart): **off-target** (cas-offinder), **INV-0** (recode / microhomology
  logic), **cryptic-splice** (MaxEntScan), **splice** (SpliceAI), **structural
  importance** (AlphaFold-DB), **2D structure** (ViennaRNA), and **docking**
  (Boltz-2 via NVIDIA). The first candidate passing all gates wins; `site_s`
  never moves.
- **Heterozygous / allele-specific**: supply phased `het_snps` (JSON array) to
  enable allele-specific design. Allele specificity is capped by natural
  heterozygosity and is never manufactured. `het_snps` already flows to the
  design driver; the het-specific *pair* pass (`design_recoded_with_het_edit`) is
  wired in a later stage — the form's "design het-specific pair" toggle records
  the request (`het_edit`) ahead of that wiring.

## Two ways to write the insert (`strategy`)

`/api/design` and `/api/design/stream` take a `strategy`. **`"tjpe"`** (the default) carries the insert inside the pegRNA
and runs `design-recoded`, which is what the section above describes. **`"pa"`** carries it in a separate donor DNA and runs
`design-pa` (Prime Assembly; [ADR 0030](../docs/decisions/0030-prime-assembly-design-pa.md)). Prime Assembly chooses its own
guides, so it takes none of the strand, arm-side, own-guide, allele-specific or anchor settings; it takes the flap length,
the number of off-target hits allowed per guide, the PBS lengths, the longest repeat, the distance between the nicks and the
number of guide pairs to try. A request is checked before a job starts, with the engine's own validator, and refused (400)
with its reason when it cannot work: a cargo shorter than two flaps (70 nt at the default flap), a setting out of range,
more than `PEGRNA_PA_MAX_PAIRS` pairs or more than 5 PBS lengths.

- `GET /api/health` lists the strategies this server has (`"strategies":["tjpe","pa"]`); a page built for Prime Assembly
  disables the choice against a server that does not list it.
- `GET /api/example/pa` serves a runnable Prime Assembly config (the 499 nt CTCF cassette).
- `GET /api/runs/{run_id}/{name}` serves a run's files. Prime Assembly writes `recoded_pa_report.txt`,
  `recoded_pa_components.csv`, `recoded_pa_changes.csv` and `recoded_pa_donor.fasta` (pegRNA A, pegRNA B and the donor).

## Modes

- **`fake`** — real design logic with stand-in scorers **and a fake PRIDICT
  runner** (no subprocess). CPU-only, no external services, always works. Runs in
  the light `pegrna-web` conda env.
- **`api`** — CPU-real gates, local PRIDICT2, plus cloud model calls, no local GPU:
  - auto-spacer: PRIDICT2 (local CPU subprocess, ~10–30 s/design; see the
    `PRIDICT2_*` env vars below)
  - splice: SpliceAI (TensorFlow-CPU, N-padded to 10 kb for short windows)
  - cryptic splice: MaxEntScan (`maxentpy`, vendored under `maxentpy_src/`)
  - structure: ViennaRNA (`vienna_structure.py`)
  - importance: live AlphaFold-DB API (`https://alphafold.ebi.ac.uk`) + biotite
  - docking: Boltz-2 via the NVIDIA-hosted NIM API
    (`https://health.api.nvidia.com/v1/biology/mit/boltz2/predict`)

  Runs in the heavier `pegrna-real` conda env. A full `api` design (PRIDICT +
  live docking) can take **5–20 min**; the app's concurrency cap (HTTP 429) and
  request timeout bound this.

## External services / credentials

- **PRIDICT2** (auto-spacer, `api` mode) runs **locally on CPU** as a subprocess —
  no network call, no GPU. Its location and interpreter are passed through from
  the env vars `PRIDICT2_DIR` / `PRIDICT2_PYTHON` / `PRIDICT2_CORES` (see the table
  below). `fake` mode uses a deterministic in-process fake runner instead.
- **AlphaFold-DB** and the **Boltz-2 NIM** are called at runtime; model weights
  are not vendored (the NIM runs remotely; AlphaFold structures are fetched and
  cached per accession). The MaxEntScan scoring matrices *are* vendored under
  `maxentpy_src/`, so the cryptic-splice gate needs no download.
- **NVIDIA API key** (for the Boltz-2 docking gate) is read at request time from
  `$NVIDIA_API_KEY` or the file `$NVIDIA_KEY_FILE`
  (default `~/.config/pegrna/nvidia.key`). It is never committed. Without it,
  `api` mode runs every other gate for real and stubs docking.
  **→ See [KEYS.md](KEYS.md) for how to obtain and install this key.**

## Running locally

```bash
# from the repo root, install the CLI + its api/model deps (see ../pyproject.toml)
pip install -e ..
# web layer:
pip install -r requirements.txt
# point at the repo (defaults to the deployment path) and serve:
PEGRNA_REPO="$(cd .. && pwd)" uvicorn app:app --host 127.0.0.1 --port 8092
```

Configuration via environment:

| var             | default                                   | meaning                          |
|-----------------|-------------------------------------------|----------------------------------|
| `PEGRNA_REPO`   | `/mnt/data3/jmoore3/prime-editing`        | path to this repo (finds CLI + examples) |
| `NVIDIA_KEY_FILE` | `~/.config/pegrna/nvidia.key`           | Boltz-2 NIM key file             |
| `NVIDIA_API_KEY`| —                                         | key via env (overrides the file) |
| `PRIDICT2_DIR`  | —                                         | PRIDICT2 checkout (auto-spacer, `api` mode) |
| `PRIDICT2_PYTHON` | —                                       | interpreter for the `pridict2` env |
| `PRIDICT2_CORES` | —                                        | CPU cores for PRIDICT2 batch design |
| `PEGRNA_PA_MAX_PAIRS` | `60`                                | most guide pairs one Prime Assembly request may try |

## Deployment

Live behind a systemd `--user` unit running uvicorn on `127.0.0.1:8092`,
exposed via a Tailscale Funnel. See the repo maintainer's notes; the code here
is deployment-agnostic.

## Layout

```
web/
  app.py            FastAPI app — endpoints, subprocess dispatch, streaming
  frontend/         React + TypeScript + Tailwind source (Vite)  -> see its README
  static/           BUILT frontend (committed, so the deploy box needs no Node)
  examples/         ready-to-load design configs
  maxentpy_src/     vendored MaxEntScan (pure-Python build) + scoring matrices
  requirements.txt  web-layer deps (FastAPI/uvicorn)
```

The UI is built from `frontend/` with `npm run build`, which writes
`static/index.html` + `static/assets/`. Those outputs are committed on purpose:
deployment stays a pure-Python affair. If you change anything under `frontend/`,
rebuild and commit the result or the server keeps serving the previous bundle.
