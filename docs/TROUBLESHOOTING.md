# Troubleshooting (HPC install + run)

Real failure modes hit while bringing the pipeline up on Imperial CX3, and how to
fix each. Ordered roughly by where they occur in the install → smoke → run flow.

If a gate is failing, first run the smoke to see *which* one and why:

```bash
# login-node gates (no GPU): splice, importance, structure(RNAformer/CPU)
python scripts/hpc/smoke_models.py --data $PEGRNA_DATA --only splice importance structure
# all 5 gates (GPU node):
qsub -v PEGRNA_DATA=$PEGRNA_DATA scripts/hpc/gpu_smoke.pbs   # then: tail -f pegrna_gpu_smoke.log
```

---

## Environment / staging

### `$EPHEMERAL` is empty over a non-login SSH channel
A non-interactive `ssh host 'cmd'` does **not** get the login shell, so `$EPHEMERAL`
(and other module-set vars) are unset — scripts then write to `/pegrna-data` or report
"file not found" for staged weights. **Fix:** wrap remote commands in a login shell
(`ssh host 'bash -lc "..."'`) or pass an absolute path
(`/rds/general/user/<user>/ephemeral/pegrna-data`).

### GRCh38 re-downloads every run (3 GB)
Older `stage_models.sh` guarded on `GRCh38.primary_assembly.genome.fa` but created
`GRCh38.fa`. **Fixed** — it now guards on `GRCh38.fa` (the name the off-target gate and
smoke expect via `PEGRNA_GENOME`). If you staged before this fix, just rename the file
to `GRCh38.fa`.

### conda/mamba not on the login-shell PATH
The env is job-script-provisioned; interactively you may need
`source $HOME/miniforge3/etc/profile.d/conda.sh && conda activate pegrna-recoded`.
`install.sh` sources conda itself.

---

## Structure gate (RiNALMo / RNAformer)

### `No module named 'flash_attn'`
RiNALMo hard-imports FlashAttention, which only builds on a GPU node with the CUDA
**dev** toolchain (`nvcc`). On a login node or a CPU-only run it is unavailable.
**This is not fatal:** `AutoStructureScorer` transparently falls back to **RNAformer on
CPU**. To actually use RiNALMo, run on a GPU node with `module load cuda/12.2` and
`pip install flash-attn --no-build-isolation` (the GPU-smoke PBS does this).

### RNAformer: `TypeError: expected str … not dict` / `'NoneType' object has no attribute 'items'`
The RNAformer config YAML was missing or was the empty `default_config.yaml` template.
RiboFormer is built from the config's dims and `load_state_dict(strict=True)`, so the
config that **matches the checkpoint** is required. **Fixed** — `stage_models.sh` now
fetches both `RNAformer_32M_state_dict_X.pth` **and** `RNAformer_32M_config_X.yml`, and
`_resolve_config` auto-pairs them by name. If you staged before this fix, fetch the
config: `wget -P $PEGRNA_DATA/rnaformer_weights <freiburg-url>/RNAformer_32M_config_biophysical.yml`.
Override explicitly with `RNAFORMER_CONFIG=/path/to/config.yml` if needed.

---

## Docking gate (Boltz-2)

### `Missing MSA's in input and --use_msa_server flag not set`
Boltz-2 needs an MSA per chain; the MSA server is unreachable on offline compute nodes.
**Fixed** — the adapter runs **single-sequence mode** (`msa: empty` per chain) by
default (`Boltz2DockScorer(offline_single_sequence=True)`, the default). If you want
server-generated MSAs on the login node instead, pass `use_msa_server=True` (needs
internet).

### `Boltz produced no confidence JSON`
Boltz couldn't find its weights on an offline node. **Fixed** — the adapter passes
`--cache $BOLTZ_CACHE`; make sure `install.sh`/`stage_models.sh` primed the cache on the
login node (`$PEGRNA_DATA/boltz_cache` should contain `boltz2_conf.ckpt`). Re-prime:
`boltz predict $PEGRNA_DATA/_prime.yaml --out_dir $PEGRNA_DATA/boltz_cache --cache $PEGRNA_DATA/boltz_cache --no_kernels`
(the `_prime.yaml` uses `msa: empty`).

### `ModuleNotFoundError: No module named 'cuequivariance_torch'` (Boltz)
Boltz-2 (>=2.2) accelerates its triangular-multiplicative updates with the
cuequivariance CUDA kernels, which it hard-imports only when kernels are enabled — and
that package needs a matching CUDA toolchain that is not in the base install.
**Fixed** — the adapter passes `--no_kernels` by default, so Boltz uses the
pure-PyTorch path (slower but portable, offline-safe). If you have the CUDA kernels
installed and want the speed-up, construct `Boltz2DockScorer(use_kernels=True)`.

---

## Off-target gate (Cas-OFFinder)

### off-target gate raises "requires the cas-offinder binary"
`cas-offinder` was not in the env. **Fixed** — it is now in `environment.yml`
(bioconda). For an env built before that: `mamba install -n pegrna-recoded -c bioconda cas-offinder`.
It runs on CPU (device `C`) — no OpenCL GPU needed. Set `PEGRNA_GENOME=$PEGRNA_DATA/GRCh38.fa`
(a full-genome CPU scan of one guide is slow; that is expected).

---

## Design run

### `NoValidPegRNA: no candidate survived all gates`
Fail-loud by design — the exception carries a per-candidate diagnostic (`diagnostics:
...`) naming the gate that killed each candidate. Common causes:
- **All candidates fail PBS2** — historically a bug (PBS2 length was fixed to PBS1's
  lengths); **fixed** by per-nick length scanning. If it still reports PBS2 misses, the
  nick's local GC genuinely admits no in-band PBS2 in 10-18 nt; widen
  `THRESHOLDS.pbs2_len_*` / band, or check that enough nicks survive off-target/arm
  screening (few surviving nicks ⇒ few PBS2 attempts).
- **Arm partition failed** — the re-templated arm crosses more than one coding boundary
  or a CDS run has no synonymous solution.
- **A real gate errored** (real mode) — each gate is guarded; the diagnostic shows the
  scorer error. Run `smoke_models.py` to confirm the gate works in isolation.

### config rejected: `unexpected keyword argument '_comment'`
The `design-recoded` CLI routes `_`-prefixed keys and provenance keys (`gene`,
`transcript`) into the report manifest; older builds passed the whole JSON as kwargs.
Update, or remove the metadata keys from the config.

---

## Git / auth (CX3 checkout)

### a token appeared in the CX3 git remote URL
Never embed a `gho_`/`ghp_` token in a remote URL. Use a **read-only SSH deploy key**:
generate a key on CX3, add it to the repo (`gh repo deploy-key add key.pub -R owner/repo`),
set `git remote set-url origin git@github.com:owner/repo.git` and pin it with
`git config core.sshCommand "ssh -i ~/.ssh/<key> -o IdentitiesOnly=yes"`. Revoke any
exposed token via GitHub → Settings → Applications → Authorized OAuth Apps.

---

## CI (before pushing)

CI runs **both** `ruff check src tests` and `black --check --line-length 100 src tests`
across py3.10/3.11/3.12, fail-fast. Run both locally before pushing:

```bash
ruff check src tests && black --check --line-length 100 src tests && pytest -q
```
