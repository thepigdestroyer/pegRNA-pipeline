# Credentials guide — getting the keys for the pegRNA web app

You only need **one secret** to run the web app in `api` mode: an **NVIDIA API
key** for the Boltz-2 docking gate. Everything else the pipeline calls is
keyless:

- **AlphaFold-DB** (importance gate) — public API, no key.
- **MaxEntScan** (cryptic-splice gate) — vendored under `maxentpy_src/`, no download.
- **SpliceAI / ViennaRNA / biotite** — pip/conda packages (see `../pyproject.toml`).

`fake` mode needs **no keys at all** — use it to smoke-test the app first.

---

## The one key you need: NVIDIA API key (Boltz-2 NIM)

The docking gate offloads the Boltz-2 pegRNA+Cas9 complex prediction to NVIDIA's
hosted NIM, so no local GPU is required. The code calls:

```
https://health.api.nvidia.com/v1/biology/mit/boltz2/predict
```

### How to get one

1. Go to the NVIDIA API catalog: **https://build.nvidia.com** and sign in
   (create a free NVIDIA account if you don't have one).
2. Open the **Boltz-2** model page: **https://build.nvidia.com/mit/boltz2**.
3. Click **“Get API Key”** (a.k.a. *Generate API Key*). Copy the value — it
   starts with **`nvapi-`**.
4. New accounts get a pool of free credits. Boltz-2 complex predictions are
   heavy (each call can take minutes and consume credits), so for sustained use
   you may need to add credits / an NVIDIA enterprise plan.

> The key is a personal secret. Never commit it, paste it into the systemd
> unit, or share it. This repo is set up so it never has to live in tracked
> files.

### Where to put it

The app reads the key at request time, in this order:

1. the `NVIDIA_API_KEY` environment variable, else
2. the file named by `NVIDIA_KEY_FILE` (default `~/.config/pegrna/nvidia.key`).

**Recommended — key file (survives restarts, stays out of process listings):**

```bash
mkdir -p ~/.config/pegrna
printf '%s' 'nvapi-XXXXXXXX...' > ~/.config/pegrna/nvidia.key   # no trailing newline needed; app strips it
chmod 600 ~/.config/pegrna/nvidia.key                            # owner-only
```

**Or — environment variable** (e.g. for a shell run):

```bash
export NVIDIA_API_KEY='nvapi-XXXXXXXX...'
```

If neither is set, `api` mode still runs **every other gate for real** and just
**stubs docking** (the UI flags this as `docking_stub`). So the key is only
required for the docking step.

### Verify it works

```bash
# with the app running on :8092, load the working example (uses api mode + docking)
curl -s http://127.0.0.1:8092/api/example/real | head -c 200
# then run a design and watch the "docking" node in the pipeline flowchart light
# up with a real ipTM score (~0.6–0.75, ~40–75s). If you see docking_stub, the
# key wasn't found.
```

You can point at a different NIM endpoint with `BOLTZ_NIM_URL` if NVIDIA moves
the model.

---

## Optional / advanced env vars

| var                | needed?  | purpose |
|--------------------|----------|---------|
| `NVIDIA_API_KEY`   | api docking | key via env (overrides the file) |
| `NVIDIA_KEY_FILE`  | api docking | path to the key file (default `~/.config/pegrna/nvidia.key`) |
| `BOLTZ_NIM_URL`    | no       | override the Boltz-2 NIM endpoint |
| `PEGRNA_CAS9_SEQ`  | no       | supply a Cas9 sequence instead of resolving nCas9→SpCas9 (UniProt Q99ZW2 +D10A) from UniProt |
| `PEGRNA_REPO`      | no       | path to this repo (default is the deployment path; set to `..` when running from `web/`) |
| `PEGRNA_TIMEOUT`   | no       | per-run subprocess timeout (seconds) |
| `PEGRNA_GENOME`    | no       | genome FASTA for the off-target gate (otherwise clean-pass) |

The `RINALMO_SS_WEIGHTS`, `RNAFORMER_CKPT`, `RNAFORMER_CONFIG`, and `BOLTZ_CACHE`
vars belong to the **local-GPU (`real`) pipeline mode**, which this web app does
**not** offer — ignore them for web deployment.
