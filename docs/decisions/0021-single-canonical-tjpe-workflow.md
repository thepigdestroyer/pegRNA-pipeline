# 0021. One canonical tj-PE design workflow (Stage H); legacy orchestrator deprecated in place

Date: 2026-07-13
Status: Accepted

## Context

The unified plan (`2026-07-02-unified-tjpe-target-strand.md`) Stage H called for ONE
tj-PE workflow, deleting `recoded_tjpe_workflow.py` and folding its features into
`tjpe_workflow.py`. Two facts have since inverted that plan:

1. **The merge direction reversed.** `recoded_tjpe_workflow.run_recoded_tjpe` became the
   feature-rich, live design workflow: auto-spacer via PRIDICT2, PRIDICT-sourced ngRNA,
   the full gate cascade (off-target / cryptic-splice / SpliceAI / importance / structure /
   Boltz-2 docking), the CDS/UTR recode arm, and the het-edit pair. It is what the CLI
   (`design-recoded`) and the web app on the DSI box actually run. `tjpe_workflow.run_tjpe_workflow`
   is the older variant-ranking orchestrator on the HPC `driver`/PBS path.
2. **The HPC path is retired as a deploy target** (owner decision 2026-07-13: update only
   the DSI box, the website, and GitHub — never the HPC). `run_tjpe_workflow` is reachable
   only through `driver.py` (the `cli run` PBS-submission path).

`tjpe_workflow.py` is also a shared-primitives hub: `STANDARD_SGRNA_SCAFFOLD`, `TJPEResult`,
`rank_variants_by_pe2`, `select_tjpe_ngrna`, `derive_rc_pbs2`, `assemble_tj_pegrna_v2`,
`NgRNADesign` are imported by 8 modules (including the canonical path and `pe7_workflow`).

## Decision

- **`run_recoded_tjpe` is the single canonical tj-PE design workflow.** Documented as such
  in both module docstrings.
- **`run_tjpe_workflow` is deprecated in place** (module note + a `DeprecationWarning`),
  not deleted. It stays functional for the HPC `driver` path and, crucially, `tjpe_workflow.py`
  remains the home of the shared primitives the canonical path depends on.
- **A convergence / no-drift guard** (`tests/test_workflow_convergence.py`) pins that
  design_recoded drives `run_recoded_tjpe`, that the legacy entrypoint warns, and that BOTH
  paths construct through `assemble_tj_pegrna_v2` (the authoritative
  `Spacer.Scaffold.rcPBS2.RTT.PBS1.PolyU` layout, ADR 0014/0020) — so they cannot re-diverge.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Delete `recoded_tjpe_workflow.py`, fold into `tjpe_workflow.py` (plan's literal Stage H) | Backwards now: recoded is the feature-rich live path; folding it into the legacy module would move the live workflow onto the deprecated HPC orchestrator |
| Delete the legacy `run_tjpe_workflow` + `driver.py` + the `run/status/...` CLI (the whole HPC-submission subsystem) | Drastic and irreversible; removes an entire (HPC) feature the owner paused but did not ask to delete; large test-migration blast radius; against the "don't touch HPC" spirit |
| Rewire `driver.py` to call `run_recoded_tjpe` | Changes the HPC-submission path's behaviour and output contract; cannot be validated without running on the HPC (explicitly out of scope) |
| Leave both undocumented | Preserves the drift risk the unified plan was written to eliminate |

## Consequences

- Positive: there is one documented canonical design workflow; the legacy path is clearly
  marked and warns on use; the shared primitives stay put so nothing breaks.
- Positive: the no-drift guard fails CI if either path stops using the authoritative
  assembler, preventing a recurrence of the rcPBS2/PBS2-class divergence bug.
- Neutral: the HPC-submission subsystem (`driver.py`, the `run` CLI family, `pe7_workflow`)
  is untouched and still importable; retiring it, if ever wanted, is a separate, explicit
  change.
- Follow-up (only if the HPC subsystem is to be removed): delete `driver.py` +
  `run_tjpe_workflow` + the HPC job-management CLI and migrate the shared primitives into a
  neutrally-named module. Deliberately deferred.

## References

- `src/pegrna_pipeline/hpc/tjpe_workflow.py` (module note + `run_tjpe_workflow` deprecation);
  `src/pegrna_pipeline/hpc/recoded_tjpe_workflow.py` (`run_recoded_tjpe`, canonical).
- `tests/test_workflow_convergence.py`. Builds on ADR 0014 (TS/rcPBS2) and ADR 0020
  (authoritative assembler). Supersedes the direction of Stage H in the 2026-07-02 plan.
