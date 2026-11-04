# 0011. MaxEntScan cryptic-splice gate; SpliceAI demoted to advisory cross-check

Date: 2026-07-02
Status: Accepted. Amended by 0028: the SpliceAI cross-check binds only on REAL genomic context,
and the `api` path's N-padded number is advisory (it had been applied as a hard gate).

## Context

With the off-target gate fixed (ADR 0010), the real CTCF run reached the splice gate and
failed every candidate: `SpliceAIScorer.delta requires >= 10000 nt of context; got
ref=243 alt=746`. This is not a bug — SpliceAI is a 10 kb context-aware CNN and correctly
fail-loud-refuses on the ~243 bp local design window. Two facts make SpliceAI the wrong
primary instrument here:

1. It needs ≥5 kb of real genomic flanks each side; the recoded workflow carries only the
   local design window, and the config has no GRCh38 coordinate anchor to fetch more.
2. The dominant splice risk for this design is a **cryptic splice site inside the ~500 bp
   exogenous cassette** (FLAG-P2A-ZeoR-T7). SpliceAI/Pangolin/MMSplice are trained on SNVs
   and small indels; a 500 bp non-genomic insertion is out of distribution, and there is no
   genomic context for the inserted bases.

The right instrument for "did the edit create a strong new splice site?" is a motif-strength
scanner (MaxEntScan; Yeo & Burge 2004) that scores a candidate donor/acceptor from a short
window — no genomic context required.

## Decision

Add a **MaxEntScan cryptic-splice gate as the primary hard splice gate**, and **demote
SpliceAI to an advisory concordance cross-check** (field-standard: concordance across
predictors, RNA validation superseding).

- New `cryptic_splice.py`: `MaxEntScorer` Protocol; `MaxEntScanScorer` (real adapter over
  `maxentpy`, CX3-only, lazy import, staged in the `[splice]` extra); `MaxEntScanFake` (CI);
  `passes_cryptic_splice_gate(ref, edited, *, scorer, strong_threshold=8.0)`.
- The gate aligns ref/edited on common prefix/suffix to locate the variant region (inserted
  + recoded bases), scans the edited construct for donor (9-mer, `score5`) and acceptor
  (23-mer, `score3`) motifs whose window overlaps the variant region, and fails if a motif
  is at/above the strong MaxEnt threshold (8.0 bits) AND stronger than any reference motif
  (created or strengthened by the edit). Native motifs entirely within the unchanged flanks
  are ignored.
- SpliceAI (`passes_splice_gate`) now runs only as a cross-check: `ValueError` (inadequate
  context) → advisory skip (no longer aborts the candidate); if a caller supplies ≥10 kb
  context and SpliceAI flags disruption, that concordant call still fails the candidate.
- `maxent_scorer` is threaded through `GateScorers` / `make_scorers` / `run_recoded_tjpe`
  (default `None` → gate skipped when not injected, preserving direct-call unit tests).
  Real bundle → `MaxEntScanScorer()`; fake bundle → `MaxEntScanFake()` (all motifs below
  threshold → passes).

Backend choice (user decision): real MaxEntScan only, staged on CX3 (`maxentpy` in the
`[splice]` extra + `environment.yml`); CI exercises the gate logic via `MaxEntScanFake`.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Keep SpliceAI as the primary gate, wire 10 kb genomic context | Needs a GRCh38 coordinate anchor the config lacks (a separate feature) AND is still out-of-distribution for the 500 bp cassette; kept as a deferred cross-check instead |
| Self-contained PSSM splice-site scorer as the backend | User chose exact MaxEntScan (real `maxentpy` on CX3) over a locally-testable approximation |
| Embed the Yeo-Burge maximum-entropy matrices in-repo | Sizable static model data + careful re-implementation; `maxentpy` on CX3 is sufficient |
| Drop splice gating for last-exon C-terminal knock-ins | The cassette can still introduce a cryptic site regardless of native exon position; a motif scan is cheap and correct |

## Consequences

- Positive: the pipeline is unblocked end-to-end — the primary splice check is context-free
  and targets the real risk (a cryptic site in the inserted/recoded sequence).
- Positive: matches the field-standard concordance model (MaxEntScan + SpliceAI), and a
  short-context SpliceAI no longer aborts the run (advisory skip).
- Neutral / follow-up: the SpliceAI cross-check stays dormant on the CTCF example until a
  ≥10 kb genomic-context feature (coordinate anchor + FASTA slice) is added — a separate
  piece of work.
- Neutral: the real MaxEntScan backend is CX3-only (`maxentpy`); the gate logic is
  CI-tested via `MaxEntScanFake`, and validated on CX3 with `maxentpy` staged.
- Note: the `strong_threshold=8.0` (Yeo & Burge "strong" site) is conservative and tunable;
  the gate is deliberately conservative (fails on any created strong site).

## References

- `src/pegrna_pipeline/hpc/cryptic_splice.py`, `splice_screen.py` (SpliceAI, now advisory),
  `recoded_tjpe_workflow.py` (gate wiring), `scorer_factory.py` (`maxent_scorer`).
- MaxEntScan: Yeo & Burge 2004 *J Comput Biol* 11:377. SpliceAI: Jaganathan 2019 *Cell*.
  Concordance / ClinGen SVI: Walker 2023 *Am J Hum Genet*.
- `docs/decisions/0010-offtarget-gate-counts-near-tier-mismatches.md`; `docs/FEATURES.md`.
- Diagnosis: CX3 job 3163779 (`cand1..18: splice scorer error: ... requires >= 10000 nt`).
