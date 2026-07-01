# 0008. Thread real-scorer resource paths into the design driver

Date: 2026-07-01
Status: Accepted

## Context

A real-mode CX3 run (`pegrna design-recoded --config examples/ctcf_cterm/config.json
--mode real`) passed the 5/5 gate self-check and then aborted the actual design with:

```
RuntimeError: real off-target scoring requires genome= (cas-offinder + the
KOLF2.1J FASTA on CX3); pass genome=/path/to/genome.fa to make_scorers
```

Root cause is a wiring gap, not a gate fault. `run_recoded_pipeline.pbs` exports
`PEGRNA_GENOME`, `BOLTZ_CACHE`, `RINALMO_SS_WEIGHTS`, and `RNAFORMER_CKPT` into the design
step's environment, and `scripts/hpc/smoke_models.py` reads those same vars — which is why
the self-check passes. But `design_recoded` built `scorer_kwargs = {"uniprot_acc", "cas9"}`
and ignored every one of them, so `make_scorers(mode="real", ...)` reached
`_real_bundle(genome=None)` and the off-target closure raised. `_real_bundle` already
accepts `genome`, `boltz_cache`, `rinalmo_weights`, `rnaformer_model_dir`, `pam`, and
`max_mismatches`; only the design entry point failed to supply them. Fixing `genome`
alone is insufficient: the next gate (Boltz-2) would then fail on an offline compute node
for lack of `boltz_cache`, wasting another GPU cycle.

Naming hazard: `design_recoded` already has a `genome` parameter, but that is the *local
~240 bp design-context sequence* consumed by `run_recoded_tjpe` (nick search, arm recode,
INV-0). It is not the genome FASTA path the off-target scorer needs.

## Decision

Give `design_recoded` explicit optional parameters for the real-scorer resources —
`genome_fasta`, `boltz_cache`, `rinalmo_weights`, `rnaformer_model_dir`, `pam`,
`max_mismatches` — named to avoid the `genome` collision. In real mode only, resolve each
from the explicit parameter **or** the environment variable the smoke harness already uses
(explicit wins), and thread them into `scorer_kwargs`:

- `genome_fasta` or `$PEGRNA_GENOME` → `genome=`
- `boltz_cache` or `$BOLTZ_CACHE` → `boltz_cache=`
- `rinalmo_weights` or `$RINALMO_SS_WEIGHTS` → `rinalmo_weights=`
- `rnaformer_model_dir` or `dirname($RNAFORMER_CKPT)` → `rnaformer_model_dir=`

`make_scorers` / `_real_bundle` are unchanged. The CLI already routes any config key that
matches a `design_recoded` parameter, so a `genome_fasta` key in a config JSON now flows
through automatically; the committed PBS needs no edit because it already exports the env.
Fail-loud is preserved: if neither the parameter nor the env var supplies the genome,
`scorer_kwargs["genome"]` is `None` and the off-target scorer still raises its helpful
error. Fake mode is untouched (`scorer_kwargs = {}`).

## Alternatives considered

| Alternative | Pros | Cons | Why rejected |
|---|---|---|---|
| Config-key-only (`genome_fasta` in JSON, no env fallback) | Explicit, self-documenting | Forces a CX3-absolute path into a committed example; PBS already exports env | Duplicates state the PBS owns; brittle across machines |
| Make `_real_bundle` read the env vars itself | One place | Hides env reads in the "pure explicit-values" factory; harder to unit-test; couples factory to CX3 env names | Keeps the factory pure; resolution belongs at the entry point |
| Fix `genome=` only | Smallest diff | Boltz then fails on offline node for lack of cache | Would waste another 6-min + GPU-hour cycle; incomplete |
| Rename/overload the existing `genome` param | Fewer params | Collides with the local design-context sequence; would break `run_recoded_tjpe` | Two different genomes; must stay distinct |

## Consequences

- Positive: the real design run consumes exactly the resources the smoke self-check
  validated, via the same env contract; zero PBS edits required.
- Positive: explicit-over-env precedence keeps the wiring unit-testable (spy on
  `make_scorers`, monkeypatch env) without importing heavy models.
- Positive: all four real-scorer resources are wired at once, so the next run does not
  die at the Boltz cache after clearing off-target.
- Neutral: `design_recoded` now reads `os.environ` in real mode; documented and confined
  to the real branch.
- Negative: six new optional parameters on `design_recoded` (all default `None`,
  backward compatible).

## References

- `src/pegrna_pipeline/hpc/recoded_driver.py` (`design_recoded`), `scorer_factory.py`
  (`_real_bundle` contract, unchanged), `scripts/hpc/run_recoded_pipeline.pbs`,
  `scripts/hpc/smoke_models.py`.
- Diagnosis: CX3 job 3158373 (`results/recoded_20260701T203205Z` empty, exit 1;
  `pegrna_recoded.log`).
- `docs/decisions/0005-gate-and-tool-stack.md`, `docs/DESIGN_HISTORY.md` §5.
