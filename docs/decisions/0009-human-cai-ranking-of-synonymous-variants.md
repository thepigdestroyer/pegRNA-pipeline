# 0009. Rank synonymous codon variants by human CAI adaptation

Date: 2026-07-01
Status: Accepted

## Context

The recoded TJ-PE workflow enumerates synonymous DNA variants of the RTT / recoded
replacement arm and hands the survivors to PRIDICT2 and full-pegRNA assembly
(`enumerate_pbs_rtt_variants` -> `enumerate_variants`). Today those variants are ranked by
`(c_count DESC, |gc - 0.55| ASC)` — a synthesis-stability proxy, not a translation-fitness
one. The requirement is that the variants selected for the recoded arm be the ones most
likely to be well expressed in human cells, i.e. closest to human codon usage. The
codebase already ships the `h_sapiens_9606` table (`HUMAN_CODON_FREQUENCIES`,
`codon_frequency` return the within-family fraction), so the reference data exists; only a
scoring + ranking hook is missing.

Two hard constraints from the existing tests: `test_codon_variant.py` pins the exact field
set of `CodonVariant` (`{dna, c_count, gc_content, rank}`), and `test_enumerate_variants.py`
pins the current `(c_count, gc)` ordering. The driver's non-recoded PRIDICT window feeders
(`driver.py` tjpe/pa/pe7) also call `enumerate_variants` and are already validated.

## Decision

Deliver the human-codon preference at the lever that actually shapes the synthesized arm:
`run_recoded_tjpe` builds the RTT from `recode_arm` (the minimum-edit synonymous recoder),
**not** from the enumerated `PbsRttVariant.rtt_dna` (a multi-perspective review found that
field is never read in the orchestrator — ranking it alone is inert). So the primary change
is a human-frequency tiebreak inside `recode_arm`:

- `_prefer_human_codons(recoded, *, original, genome, frame, max_run)` — for every codon
  position the recode **changed** (to break microhomology), swap in the most frequent human
  synonymous codon that still keeps `longest_common_run(., genome) <= max_run`. Positions
  equal to `original` are left byte-exact, so this **never adds edits** and a no-op recode
  stays a no-op (min-edits principle and Decision 15 non-coding byte-exactness preserved).
  Iterated to a fixpoint (each swap strictly raises that codon's frequency → terminates), so
  the result is locally frequency-optimal at the changed positions given the run constraint.
- Applied in the `recode_arm` dispatcher so **both** backends (greedy `recode_min_edits`
  and DnaChisel `recode_min_edits_dnachisel`) get the preference uniformly, without tuning
  fragile DnaChisel objective weights.

Add the supporting primitive `human_cai(dna, frame=0)` to `codon_expander.py`: the Codon
Adaptation Index against the human genome-wide table — geometric mean of per-codon relative
adaptiveness `w = f(codon) / max_{c in family} f(c)`. Single-codon families (Met `ATG`, Trp
`TGG`) and stop codons carry `w = 1` and are excluded (no discriminating signal). The shipped
table is 2-decimal-rounded, so a rare codon can read `0.00`; `w` is floored
(`_CAI_ZERO_FLOOR = 0.01`) before the log so the index degrades smoothly instead of raising
`math domain error`.

Also add an opt-in `ranking` parameter to `enumerate_variants` / `enumerate_with_metadata`
(`"gc_c"` default = the existing `(c_count DESC, |gc-0.55| ASC, dna)` order, byte-for-byte
unchanged for the driver's tjpe/pa/pe7 feeders; `"cai"` = CAI DESC with the gc_c order as
tiebreak) and default `enumerate_pbs_rtt_variants(ranking="cai")`. This is a **complementary**
reusable primitive — correct for any future consumer and for surfacing CAI in the report —
but it is NOT the mechanism that changes the synthesized arm; `_prefer_human_codons` is.
`CodonVariant` is unchanged; CAI is recomputed on demand (`human_cai(v.dna)`).

## Alternatives considered

| Alternative | Pros | Cons | Why rejected |
|---|---|---|---|
| Rank `enumerate_pbs_rtt_variants` only (the first attempt) | Smallest change | `run_recoded_tjpe` never reads `rtt_dna`; the arm is built from `recode_arm`, so the ranking is inert in the production path (caught by review) | Does not change the synthesized arm; the real lever is `recode_arm` |
| Full max-CAI recode (DnaChisel `CodonOptimize` objective) | Maximally human-adapted arm | Adds edits beyond the microhomology break; supersedes the min-edits decision; changes the no-op case; DnaChisel objective-weight tuning is fragile | User chose the tiebreak-within-min-edits option; keeps min-edits + Decision 15 intact |
| Frequency-harmonisation (match the per-family distribution) | Most literal "same frequency as genome" | For a short arm each aa appears ~once, so it collapses to "prefer the frequent codon" anyway; less standard | User chose CAI adaptation ("most likely expressed"); converges for short arms |
| Add a `cai` field to `CodonVariant` | Report can read it directly | Breaks the pinned field-set test + changes the frozen dataclass shape | Recompute on demand; keep the dataclass and its tests intact |
| Highly-expressed-gene CAI reference (Sharp & Li) | Classic max-expression index | Needs a bundled reference gene set; not "the human genome" the request named | Genome-wide table is the requested reference and is self-contained (offline-safe) |
| Change the default `enumerate_variants` ranking to CAI | One behaviour | Breaks `test_enumerate_variants`; silently shifts validated tjpe/pa/pe7 driver outputs | Opt-in is backward compatible |

## Consequences

- Positive: the **synthesized** recoded arm now lands its forced microhomology-breaking
  edits on the most frequent human synonymous codons, so the recode is translation-competent
  — delivered at `recode_arm`, the lever `run_recoded_tjpe` actually uses.
- Positive: min-edits and Decision 15 (non-coding byte-exact) are preserved — the tiebreak
  never adds edits and never touches an unchanged position; a no-op recode stays a no-op.
- Positive: fully backward compatible — default `enumerate_variants` ranking, `CodonVariant`,
  and the tjpe/pa/pe7 driver feeders are unchanged; all existing tests pass untouched.
- Positive: `human_cai` is a pure, self-contained function (no reference FASTA, offline-safe
  on CX3), reusable for reporting and any future consumer.
- Neutral: the reference is the genome-wide `h_sapiens_9606` table (per the request), not a
  highly-expressed-gene subset; documented so a future max-expression variant is a known
  option. Cross-codon selection is not frequency-aware (the greedy still picks *which* codon
  to change on run-reduction; the tiebreak makes each *changed* codon human-optimal).
- Negative: `recode_arm` now runs a bounded fixpoint pass calling `human_cai`/
  `longest_common_run` over the changed codons; negligible for the short recoded arms.

## References

- `src/pegrna_pipeline/hpc/recode_arm.py` (`_prefer_human_codons`, `recode_arm`),
  `codon_expander.py` (`human_cai`, `enumerate_variants`), `recoded_tjpe_workflow.py`
  (`_partition_arm` -> `recode_arm`; `enumerate_pbs_rtt_variants`).
- `docs/decisions/0005-gate-and-tool-stack.md`; `docs/DESIGN_HISTORY.md`.
- CAI: Sharp & Li 1987, NAR 15:1281 (relative adaptiveness, geometric mean).
