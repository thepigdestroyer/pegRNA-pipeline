# 0016. Two-pegRNA PRIDICT adapter: no off-target/RTgDNA fields; RC(TS) orientation

Date: 2026-07-03
Status: Accepted

## Context

Stage C assembles the recoded tj-pegRNA from **two** PRIDICT2 designs and binds an ngRNA. Two
Stage-C assumptions (plan §15 items 1 and 4) were unverified: whether PRIDICT2 exposes a
per-spacer off-target-risk field, and the name of the column carrying the PBS-past-the-nick
("RTgDNA") homology. Both were verified against the real PRIDICT2 install on CX3
(`/rds/general/user/mj1824/home/PRIDICT2`).

**CX3 findings (grep of the PRIDICT2 source):**
- **No off-target-risk field** — no `offtarget`/`off_target`/`off-target` references anywhere in
  the PRIDICT2 `.py` sources.
- **No separate "RTgDNA" column** — no `rtgdna`/`homology`/`past-nick` references; `pridict_parser`
  parses only `RTrevcomp`/`RTlength`. The RT template *is* `RTrevcomp`; its genomic-homology tail
  (the part past the insert) is the "RTgDNA".

## Decision

- **`TwoPegDesign`** (frozen) carries the assembled components (all RC(TS)): `pbs1`, `first_end`
  (from pegRNA#1), `rcpbs2`, `last_end` (from pegRNA#2), the bound `ngrna`, lengths, editing
  position, target strand, and `offtarget_risk: Optional[float] = None`.
- **Orientation (`combine_two_peg`, C3):** pegRNA#1 (designed on TS) → `pbs1 = peg1.PBSrevcomp`,
  `first_end = peg1.RTrevcomp`, **untouched** (PRIDICT already emits RC(TS)). pegRNA#2 (designed on
  the opposite strand) → `rcpbs2 = RC(peg2.PBSrevcomp)`, `last_end = RC(peg2.RTrevcomp)` (both RC(TS)),
  and `ngrna = peg2.spacer + scaffold`. The RC(TS)-end-to-end property is asserted by a test with a
  known target strand.
- **ngRNA lifecycle (C4):** the ngRNA is a field of `TwoPegDesign` (bound); discarding a candidate
  discards its ngRNA (one immutable object). A new pegRNA#2 yields a new design with a new ngRNA.
- **`offtarget_risk` defaults to `None`** — PRIDICT2 has no such field, so spacer specificity is
  gated by Cas-OFFinder (Stage D keeps the Cas-OFFinder path; the "replace with PRIDICT off-target
  risk" D1 branch is not available).
- **`TwoPegPridictFake`** drives the pipeline in CI (deterministic, shaped like `parse_pegrna_csv`);
  the **real adapter (C2)** wraps `submit_pridict_batches` + `pridict_parser` and is validated on CX3
  in a real run (deferred).

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Use a PRIDICT off-target-risk field for spacer QC (plan D1) | The field does not exist in PRIDICT2 (CX3-verified); Cas-OFFinder is the real specificity gate |
| Parse a dedicated "RTgDNA" column for C5 exclusion | No such column exists; the genomic tail lives inside `RTrevcomp`, so C5 must operate on `RTrevcomp` |
| Re-orient pegRNA#1 as well | PRIDICT already emits pegRNA#1 as RC(TS); re-complementing would double-invert it |

## Consequences

- Positive: the two-pegRNA design is a small, testable transform. The RC(TS) property is asserted
  at the **unit level** (the transform + fake); the *biological* rule that PRIDICT-on-opposite-strand
  emits `PBSrevcomp = RC(OS) = TS` is verified against **real PRIDICT output in C2** (not here). The
  ngRNA is structurally bound (cannot drift from its tj-pegRNA). `TwoPegDesign` is a design
  representation only and is **not** wired into real construct assembly until C5 + C2 land.
- Positive: the off-target and RTgDNA unknowns are resolved from the real tool, not guessed.
- Follow-up (this Stage C, next sub-PR): **C5** RT-tail (RTgDNA) exclusion — where the genomic tail
  of `RTrevcomp` overlaps the recoded between-nicks arm, defer to the arm (it covers that span to
  the nucleotide). **C2** real adapter execution — validated on CX3.
- Neutral: `TwoPegPridict.design_pair` takes pre-split windows for now; the real adapter's editseq
  construction from an insert + locus is C2.

## References

- CX3: `/rds/general/user/mj1824/home/PRIDICT2` (grep: no offtarget/RTgDNA refs).
- `src/pegrna_pipeline/hpc/recoded_pridict.py` (`TwoPegDesign`, `TwoPegPridictFake`,
  `combine_two_peg`); `pridict_parser.py`, `pridict_submitter.py`, `data_types.PegRNADesign`.
- `tests/test_recoded_pridict.py`. Builds on ADR 0014 (target-strand convention). Plan Stage C
  in `docs/superpowers/plans/2026-07-02-unified-tjpe-target-strand.md`.
