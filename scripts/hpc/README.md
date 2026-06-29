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

## 1. Environment (login node, has internet)

```bash
module load anaconda3/personal
conda env create -f scripts/hpc/environment.yml
conda activate pegrna-recoded
pip install -e '.[recode,splice,importance,fullstructure,complex]'
```

## 2. Stage weights + reference data (login node, once)

```bash
export PEGRNA_DATA=$EPHEMERAL/pegrna-data
bash scripts/hpc/stage_models.sh        # T2T-KOLF genome, SpliceAI, RiNALMo,
                                         # RNAformer, Boltz cache, AFDB dir
```

Resolve every `TODO`/`WARN` the script prints (exact Zenodo asset names, RiNALMo
SS-head checkpoint, RNAformer ckpt, GRCh38↔KOLF chain). They are pinned to
`$PEGRNA_DATA/MANIFEST.txt` for run provenance.

## 3. Validate every real model (the gate)

```bash
python scripts/hpc/smoke_models.py --data $PEGRNA_DATA
```

Must report **5/5 passed**. This is where any mistake in the documented-API
inference bodies surfaces (wrong tensor shape, wrong JSON key, wrong checkpoint
path) — fix against the installed versions until green. `run_recoded_pipeline.pbs`
re-runs this and aborts the job if it fails.

## 4. Build a design config

JSON consumed by `pegrna design-recoded`:

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
qsub -v CONFIG=/path/to/design_config.json,PEGRNA_DATA=$PEGRNA_DATA \
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
