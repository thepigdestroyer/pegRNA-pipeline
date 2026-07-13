# 0023. Run the reverse-targeting + het-specific designs in parallel (Stage G)

Date: 2026-07-13
Status: Accepted

## Context

`design_recoded_with_het_edit` produces two designs for one insert + site: the
reverse-targeting primary (antisense recode-into-CDS) and, when a discriminating phased
het is present, the allele-specific sense-window het pass. The unified plan Stage G (G3)
calls for running them in parallel. They were sequential: primary first, then het.

The two designs are independent -- each is a full `design_recoded` call with its own
genome/cds_mask/strand frame, its own scorer bundle, its own PRIDICT input (unique uuid),
and its own output subdir (`primary/` vs `het_edit/`). `find_allele_specific(het_snps)` (the
cheap, design-free check that decides whether the het pass runs at all) does not depend on
the primary result.

## Decision

- **Dispatch the two designs concurrently on a `ThreadPoolExecutor(max_workers=2)`** when
  both are needed. `find_allele_specific` is evaluated first (no design), so the primary-only
  paths (no het window, or no discriminating het) stay single-threaded and the het pass is
  never dispatched.
- **The primary is mandatory** -- its failure propagates (as before). **The het pass stays
  advisory** -- a `NoValidPegRNA` from it is caught and reported as a note, never killing the
  successful primary.
- Threads (not processes): the designs are largely I/O-bound (PRIDICT subprocess, Boltz-2
  NIM, AlphaFold-DB) with separate output dirs and separate scorer objects, so threads are
  safe and avoid process-spawn/pickling overhead. Each SpliceAI/`_PaddedSpliceAIScorer` is a
  distinct model object per design, so concurrent inference does not share TF model state.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Keep sequential | Leaves ~half the wall-clock on the table when both passes run (plan G3) |
| `ProcessPoolExecutor` | Heavier (pickling design kwargs + results, a second PRIDICT process tree); the work is I/O-bound so threads suffice |
| Parallelize inside a single design | Out of scope; the two designs are the natural independent unit |

## Consequences

- Positive: when both passes run, they overlap, roughly halving wall-clock for a het-edit
  request. Results are identical to the sequential order (independent designs, separate
  output dirs) -- pinned by a determinism test.
- Neutral: a primary failure now lets the (already-dispatched) het pass finish before the
  exception propagates -- wasteful only in the rare primary-failure case.
- Caveat (api mode): two designs load their own SpliceAI (TF) models and call Boltz-2 NIM
  concurrently; the web app's outer concurrency cap still bounds total load. Watch the first
  live api het run for TF-thread issues (none expected -- separate model objects).

## References

- `src/pegrna_pipeline/hpc/recoded_driver.py` (`design_recoded_with_het_edit`).
- `tests/test_het_parallel.py` (pair produced, determinism, primary-only paths).
  Builds on Stage 3c (het-edit wiring) and the unified plan Stage G.
