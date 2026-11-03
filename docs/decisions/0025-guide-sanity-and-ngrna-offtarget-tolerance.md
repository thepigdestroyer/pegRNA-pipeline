# 0025. One definition of a usable guide; an explicit, disclosed off-target tolerance for the second nick

Date: 2026-11-03
Status: Accepted. Refines 0010 (its default is unchanged) and settles two findings from 0024.

## Context

ADR 0024 measured the gates and left two things alone because changing them needed a decision,
not a fix. Both were filed as issues (#24, #25) and are decided here.

**1. Spacer sanity was applied on one path only.** `select_tjpe_ngrna` (the PRIDICT path)
rejected a spacer containing `TTTT`, a homopolymer run of 5 or more, or GC outside 30-75%. The
genome-enumerated path applied none of it. So the same real guide, `GCTGCCACAGATGCCCCCAA`
(nick 101 of the CTCF example, the demo's only clean candidate, which contains `CCCCC`), was
rejected on one path and accepted on the other. And a guide with a `TTTT` run, which RNA
polymerase III reads as a stop signal, could be chosen on the enumerated path (nick 232's
`TTCTCCCTTCTTTCTTTTTT` would have been, had it been off-target clean).

**2. The ngRNA off-target tolerance was a dead parameter.** Both screens rejected on
`has_off or off > max_offtargets or on < min_ontarget`, and the scorers set `has_off` to
`off > 0`. Any guide with a single off-target was therefore rejected whatever `max_offtargets`
was, so the knob that was supposed to relax the zero tolerance of ADR 0010 could not. The
report also said "0 off-targets in tier" as a fixed string, true only because nothing else was
possible.

A third thing surfaced while checking how a reader would see any of this: the web UI declared a
`warnings` field and rendered none of it. Every advisory the pipeline emits, including 0024's
"docking gate is advisory-only for this design", was visible only inside the collapsed full
text report.

## Decision

**Guide sanity is one module, `guide_sanity`, called by every site that picks a guide.**

* HARD rules (a failing spacer is never chosen): `TTTT`, and GC outside 30-75%.
* ADVISORY (reported on the chosen design, never a rejection): a homopolymer run of 5 or more.

Only `TTTT` has a mechanism: Pol III terminates there. A run of five C, A or G has none, and no
primary evidence was found that it prevents a guide being made or working; general design
guides advise against runs of C and G but cite nothing specific. That is why it is advisory,
and why the two paths now agree. The rule is applied on the PRIDICT path (`select_tjpe_ngrna`),
on the genome path before the whole-genome scan (`filter_nick_sanity`, so an unusable guide
does not cost a scan), and the chosen pegRNA spacer is assessed too (a failure there is a
warning, since it is the user's or PRIDICT's choice). It only constrains guides expressed from
a U6 promoter; a synthetic pegRNA is unaffected by `TTTT`, but the pipeline stays conservative.

**The off-target rule is one function, `offtarget_rejection`, used by both paths.**

* The count is compared with the limit. The flag is not consulted on its own.
* A guide flagged with a count of zero is contradictory and is rejected.
* **Exact duplicates of the guide elsewhere in the genome are never tolerated.** The combined
  count lumps near-matches (1-2 mismatches: a weak, rare cutting site) with perfect copies
  (which cut as well as the intended site). A tolerance for the first must not admit the second,
  so `CasOffinderScorer` now returns an `OffTargetScore`: a real 3-tuple (it compares equal to
  and unpacks like the old triple, so no consumer or test changed) that also carries
  `.duplicates`.

**The tolerance is `ngrna_max_offtargets`, default 0 (ADR 0010 unchanged).** It is a parameter
of `design_recoded`, so a JSON config or the CLI sets it, and the web UI's Fine-tuning panel has
the field (0-10, with 0 marked as recommended). Whatever the chosen guide actually has is
reported: the off-target gate row states the count and the tolerance, and a warning appears
whenever it is above 0.

**Warnings are shown.** A "Check these before you order" callout lists them above the
collapsed report.

## Alternatives considered

| Alternative | Why not |
|---|---|
| Make homopolymers a hard rule on both paths | Consistent, but it rejects guides with no evidence of harm, including the demo's only clean one |
| Drop the sanity check on both paths | Lets `TTTT` through |
| Raise the default tolerance | 0010 chose zero deliberately and nothing here shows it harming a design. The knob exists; the default is a safety choice |
| Return a 4-tuple to carry duplicates | Breaks every consumer and `test_exact_duplicate_counts_as_offtarget`, which pins the 3-tuple |
| Treat a duplicate as "no unique on-target" (`on = 0`) | Same breakage, and the rejection message would be misleading |
| A position-aware CFD score instead of a count | Still the better measure, still deferred: it needs the mismatch positions of every hit |

## Consequences

* The demo is unaffected in outcome (nick 101 is still chosen) but is now chosen the same way
  on both paths.
* On the PRIDICT path a guide with a homopolymer is no longer dropped, so the DeepCas9-ordered
  list can differ from before; the chosen guide's advisories are reported.
* A single-strand nick at a near-match is generally repaired without mutation, which is the
  reasoning behind allowing a tolerance for a nicking guide at all. That is reasoning, not
  something measured here, and it is why the default stays 0 and a non-zero value is always
  disclosed.
* **Still open:** a count weighs a seed-region mismatch the same as a distal one, although they
  are very different risks. CFD scoring would fix that.
* Evidence that the tests can fail: restoring the old rule makes six of the tolerance tests
  fail; a guide with a perfect duplicate is rejected at every tolerance up to 10.

## References

* `src/pegrna_pipeline/hpc/guide_sanity.py`, `nick_search.py` (`offtarget_rejection`,
  `filter_nick_sanity`), `offtarget_scorer.py` (`OffTargetScore`), `pridict_ngrna.py`,
  `recoded_tjpe_workflow.py` (`ngrna_max_offtargets`), `recoded_driver.py`.
* Tests: `test_guide_sanity.py`, `test_ngrna_offtarget_tolerance.py`.
* ADR 0010 (the zero-tolerance default and its rejected alternatives), ADR 0024 (findings 2 and 3).
