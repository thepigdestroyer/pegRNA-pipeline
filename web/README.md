# pegRNA design — web frontend

A FastAPI web app that wraps the `pegrna design-recoded` CLI (this repo). The
browser posts a design config; the server writes it to a per-run directory, runs
the pipeline as a subprocess, streams per-gate progress back, and renders the
text report + parsed CSV components.

## Modes

- **`fake`** — real design logic with stand-in scorers. CPU-only, no external
  services, always works. Runs in the light `pegrna-web` conda env.
- **`api`** — CPU-real gates plus cloud model calls, no local GPU:
  - splice: SpliceAI (TensorFlow-CPU, N-padded to 10 kb for short windows)
  - cryptic splice: MaxEntScan (`maxentpy`, vendored under `maxentpy_src/`)
  - structure: ViennaRNA (`vienna_structure.py`)
  - importance: live AlphaFold-DB API (`https://alphafold.ebi.ac.uk`) + biotite
  - docking: Boltz-2 via the NVIDIA-hosted NIM API
    (`https://health.api.nvidia.com/v1/biology/mit/boltz2/predict`)

  Runs in the heavier `pegrna-real` conda env.

## External services / credentials

- **AlphaFold-DB** and the **Boltz-2 NIM** are called at runtime; model weights
  are not vendored (the NIM runs remotely; AlphaFold structures are fetched and
  cached per accession). The MaxEntScan scoring matrices *are* vendored under
  `maxentpy_src/`, so the cryptic-splice gate needs no download.
- **NVIDIA API key** (for the Boltz-2 docking gate) is read at request time from
  `$NVIDIA_API_KEY` or the file `$NVIDIA_KEY_FILE`
  (default `~/.config/pegrna/nvidia.key`). It is never committed. Without it,
  `api` mode runs every other gate for real and stubs docking.

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

## Deployment

Live behind a systemd `--user` unit running uvicorn on `127.0.0.1:8092`,
exposed via a Tailscale Funnel. See the repo maintainer's notes; the code here
is deployment-agnostic.

## Layout

```
web/
  app.py            FastAPI app — endpoints, subprocess dispatch, streaming
  static/index.html vanilla-JS SPA (form, pipeline flowchart, report/CSV view)
  examples/         ready-to-load design configs
  maxentpy_src/     vendored MaxEntScan (pure-Python build) + scoring matrices
  requirements.txt  web-layer deps (FastAPI/uvicorn)
```
