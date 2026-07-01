# Recoded TJ-PE pipeline — Imperial CX3 runbook

End-to-end run of the recoded TJ-pegRNA design pipeline with **every gate model
active** (SpliceAI, AlphaFold-importance, RiNALMo, RNAformer, Boltz-2) on the
Imperial CX3 HPC.

> **Status — read first.** The design logic, all gates, the orchestrator, the
> strand/TPRT verifier, T7/T7-rev elements, the CLI, and the report layer are
> implemented and unit-tested in CI (fakes for the heavy models). The **real
> model inference bodies are written against each tool's documented API but were
> NOT executable on the dev box** (no Ampere GPU / no weights). `smoke_models.py`
> is the gate that proves each one runs on CX3 before any production design.
> **The pipeline is production-ready ONLY after `smoke_models.py` passes 5/5 on
> CX3** and the calibration thresholds (below) are set. Until then it is
> validated-in-fake-mode, staged-for-real.

## 0. Connect

CX3 is behind ZScaler ZPA + SAML SSO (interactive — must be done by a human).
Use the `hpc-connect` skill or your usual `ssh` flow.

## 1. Install (one command, login node, has internet)

```bash
git clone https://github.com/thepigdestroyer/pegRNA-pipeline && cd pegRNA-pipeline
export PEGRNA_DATA=$EPHEMERAL/pegrna-data
bash scripts/hpc/install.sh
```

This does everything in steps 1-3 below (env create → `pip install -e '.[…all extras…]'`
→ stage weights/genome → login-node smoke). Idempotent. If you prefer the manual steps
(or need to debug one), they are:

```bash
# 1a. env
module load anaconda3/personal          # or: source ~/miniforge3/etc/profile.d/conda.sh
conda env create -f scripts/hpc/environment.yml   # includes cas-offinder
conda activate pegrna-recoded
pip install -e '.[recode,splice,importance,fullstructure,complex]'
# 1b. stage weights + reference genome (T2T-KOLF/GRCh38, SpliceAI, RiNALMo, RNAformer
#     + its matching config .yml, Boltz cache, AFDB dir)
export PEGRNA_DATA=$EPHEMERAL/pegrna-data
bash scripts/hpc/stage_models.sh
```

Resolve every `WARN` the staging script prints; pins go to `$PEGRNA_DATA/MANIFEST.txt`
for run provenance. If a model fails to load, see
[`docs/TROUBLESHOOTING.md`](../../docs/TROUBLESHOOTING.md).

## 3. Validate every real model (the gate)

```bash
# login-node CPU-capable gates (splice, importance, structure via RNAformer):
python scripts/hpc/smoke_models.py --data $PEGRNA_DATA --only splice importance structure
# all 5 gates on a GPU node (incl. RiNALMo + Boltz-2):
qsub -v PEGRNA_DATA=$PEGRNA_DATA scripts/hpc/gpu_smoke.pbs   # target: 5/5
```

`gpu_smoke.pbs` must report **5/5 passed**. This is where any mistake in the
documented-API inference bodies surfaces (wrong tensor shape, JSON key, checkpoint
path). `run_recoded_pipeline.pbs` re-runs the smoke and aborts the job if it fails.

## 4. Build a design config (or use the shipped CTCF example)

A complete, runnable example ships at
[`examples/ctcf_cterm/config.json`](../../examples/ctcf_cterm/config.json) (CTCF
C-terminal FLAG-P2A-ZeoR knock-in, real chr16+ locus). The config is a JSON of the
design inputs; `_`-prefixed keys and provenance keys (`gene`, `transcript`,
`uniprot_acc`) are routed into the report manifest:

```json
{
  "spacer": "….20nt fixed by site S….",
  "rtt_aa": "MKP…",
  "cargo": "<assemble via insert_elements.assemble_insert, may include T7/T7_rev>",
  "genome": "<T2T-KOLF locus window, sense 5'->3'>",
  "site_s": 150,
  "cds_mask": [true, true, "…per-base CDS vs non-coding…"],
  "pegrna_strand": "Rv",
  "pbs_lengths": [13, 15, 17],
  "allele_hunt": false,
  "het_snps": [],
  "uniprot_acc": "<target protein, for AlphaFold gate>",
  "cas9": "<nCas9 aa sequence, for Boltz-2 docking>",
  "manifest": {"genome_build": "T2T-KOLF2.1Jv1.1", "git_commit": "…"}
}
```

`cargo` with a T7 promoter, e.g. `FLAG-GGGGS-T7-P2A-ZeoR` (forward T7 transcribes
`P2A-ZeoR`) or `…-T7_rev-…` (transcribes `GGGGS-FLAG` upstream) — build with
`pegrna_pipeline.hpc.insert_elements.assemble_insert`.

## 5. Submit

```bash
qsub -v CONFIG=examples/ctcf_cterm/config.json,PEGRNA_DATA=$PEGRNA_DATA \
     scripts/hpc/run_recoded_pipeline.pbs
```

Output: `results/recoded_<ts>/recoded_tjpe_report.txt` (design + every gate
result + allele-specific options + provenance manifest).

## 6. Dry-run anywhere (no GPU)

```bash
pegrna design-recoded --config design_config.json --mode fake --output-dir results/dry
```

Runs the full design logic with stand-in scorers — for plumbing/CI, NOT a
biological result.

## Calibration parameters (set before trusting outputs — wet-lab)

| Parameter | Where | Default | Note |
|---|---|---|---|
| microhomology max-run `L` | `run_recoded_tjpe(max_run=)` | 8 | conservative; validate |
| RiNALMo pairing threshold | `RiNALMoScorer(threshold=)` | 0.5 | val-tune 0.01–0.30 |
| Boltz-2 interface-drop | `passes_docking_gate(max_drop=)` | 0.15 | no published cutoff |
| splice Δ | `passes_splice_gate(max_delta=)` | 0.2 | SpliceAI convention |
| nick distance | `ngrna_distance_*` | 40–70 | <70 bp untested (Zheng) |

## Definition of production-ready

1. `smoke_models.py` → 5/5 on CX3.
2. Stage-0 real T2T-KOLF diploid + GRCh38 chain wired (replacing fixtures).
3. Calibration thresholds set from a wet-lab calibration library.
4. A full `--mode real` run completes and the report's INV-0 / strand checks pass.

Until all four hold, treat outputs as engineering previews, not bench-ready designs.
