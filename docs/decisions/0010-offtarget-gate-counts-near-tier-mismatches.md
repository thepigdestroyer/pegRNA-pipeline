# 0010. Off-target gate counts near-tier mismatches, not distant hits

Date: 2026-07-01
Status: Accepted

## Context

The first real-mode CTCF run to clear the genome-wiring fix (ADR 0008) still failed:
`NoValidPegRNA: no candidate survived all gates after 0 tried; diagnostics:` (empty). Root
cause, gathered on CX3 with the real `cas-offinder` + GRCh38:

- `enumerate_nick_sites` returns 6 nicks for the CTCF context; the pipeline works in fake
  mode, so enumeration is fine.
- `CasOffinderScorer` counted **every** genomic hit within `max_mismatches=3` as an
  off-target and `screen_nick_candidates` requires `max_offtargets=0`. Against 3.1 Gb
  GRCh38 every 20-mer + NGG has dozens of 3-mismatch sites, so every nick was dropped:

  | nick | mm0 | mm1 | mm2 | mm3 |
  |---|---|---|---|---|
  | pos=101 | 1 | 0 | 0 | 17 |
  | pos=67 | 1 | 0 | 5 | 36 |
  | pos=61 | 1 | 0 | 1 | 22 |
  | pos=70 | 1 | 2 | 23 | 423 |

  Every guide is unique at 0 mismatches (a good on-target); the load is almost entirely
  distant 3-mismatch hits, which are not functional off-targets. pos=101 is a textbook-clean
  guide (0 off-targets within 2 mismatches) yet was rejected over 17 distant 3mm sites.

- An empty nick pool means the candidate loop never runs, so the final `NoValidPegRNA`
  carried a blank diagnostic — the failure was unexplained.

The `max_offtargets=0` gate was calibrated against the fake scorer (which returns 0
off-targets); it is biologically wrong against a real genome searched to 3 mismatches.

## Decision

Count only **disqualifying** off-targets, decoupling the gate stringency from the search
depth:

- `CasOffinderScorer` gains `count_mismatch_threshold` (default **2**). Off-targets =
  extra EXACT (0-mismatch) duplicates + hits with `0 < mismatches <= count_mismatch_threshold`.
  Hits beyond the threshold are still **searched** up to `max_mismatches` (a full profile is
  available) but are **not counted** — so distant 3-mismatch noise no longer rejects a guide.
- `screen_nick_candidates` is unchanged (`max_offtargets=0`, `min_ontarget=0.5`); with the
  corrected count, "0 off-targets" now means "0 off-targets within 2 mismatches", a
  conservative, standard gate. pos=101 passes; the repetitive/low-complexity guides
  (pos=70/67/61, which have ≤2mm off-targets) are still rejected.
- `run_recoded_tjpe` appends an explanatory diagnostic when the nick pool is empty
  (distinguishing "0 nick sites in window" from "0 nicks after off-target screen") so a
  future 0-candidate run is never silent.

## Alternatives considered

| Alternative | Pros | Cons | Why rejected |
|---|---|---|---|
| Raise `max_offtargets` above 0 | One-line | A magic count that is mismatch-agnostic (2 exact duplicates vs 400 distant 3mm hits are not comparable) | Counting the right tier is more correct than allowing an arbitrary total |
| Search only to `max_mismatches=2` | Simplest | Conflates search depth with gate stringency; discards the 3mm profile that a report may want | Decoupling is cleaner and keeps characterisation |
| CFD specificity score (Doench 2016) | Field-standard, position/type-weighted | Needs the CFD matrix + per-hit mismatch positions; larger change; activity is already `score_guides.R`'s job | Deferred; the near-tier count is a conservative, sufficient proxy for the specificity gate |

## Consequences

- Positive: the real CTCF design gets a clean guide (pos=101) past the off-target gate, so
  candidates are actually evaluated instead of a silent 0-tried failure.
- Positive: the gate now reflects real off-target risk (near-tier hits), configurable via
  `count_mismatch_threshold`; fake mode and the smoke harness are unaffected (they report 0 /
  a single guide).
- Positive: an empty nick pool now explains itself.
- Neutral: this scorer covers SPECIFICITY only; on-target ACTIVITY (DeepSpCas9/DeepHF) is
  still `score_guides.R`'s job. A downstream gate (e.g. Boltz-2 docking ipTM floor) may still
  block a given locus — a separate concern, and now visible in the diagnostics.

## References

- `src/pegrna_pipeline/hpc/offtarget_scorer.py` (`CasOffinderScorer`),
  `nick_search.py` (`screen_nick_candidates`), `recoded_tjpe_workflow.py` (nick diagnostic).
- Evidence: CX3 `cas-offinder` on GRCh38 for the CTCF nicks (job 3158586 follow-up).
- `docs/decisions/0008-thread-real-scorer-resources-into-design-driver.md`; `docs/FEATURES.md`.
