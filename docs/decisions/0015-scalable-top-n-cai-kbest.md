# 0015. Scalable top-N-by-CAI via lazy k-best (not full enumeration)

Date: 2026-07-02
Status: Accepted

## Context

Stage B must hand PRIDICT2 the **top 1000-5000 synonymous variants of an insert ranked by
CAI**, and rank the insert middle by CAI independently. The existing `enumerate_variants`
walks the full `itertools.product` of per-codon choices. That is intractable at the scale the
pipeline targets: a 40 bp insert end is ~13 codons (10^5-10^7 synonymous combos), and a
<=500 bp middle is ~140 codons (~3^140 combos). Full enumeration to then take the top-N cannot
run.

Key structural fact: ranking synonymous variants by CAI is **separable**. `human_cai =
exp((1/K) * sum_i log w_i)` where `w = f / f_max` within each aa family and `K` (the number of
non-excluded codons -- single-codon families Met/Trp and stops are dropped) is fixed by the aa
sequence, not by codon choice. `exp` and dividing by the constant `K` are monotonic, so
**CAI-descending order == descending sum of per-codon log-w** -- a sum of independent per-
position terms.

## Decision

Implement a **lazy best-first (Dijkstra) k-best search over the codon-index lattice**:

- Each position's codons are sorted by CAI weight desc; per-choice `log w` precomputed.
- State = tuple of per-position indices; start = all-best (max score). A min-heap keyed by
  `(-score, -c_count, |gc-0.55|, dna)` pops states in CAI-descending order (with a deterministic
  tiebreak). Successors advance one index; a `visited` set dedupes the multiple lattice paths to
  a state. Because every edge only lowers the score, best-first pops in non-increasing CAI --
  provably correct and each state emitted once.
- Quality filters (internal stop, CpG, restriction site, homopolymer) are **non-separable**, so
  they gate *yielding* (not *exploration*): a filtered state is skipped but its successors are
  still expanded (a passing variant may sit lattice-adjacent to a filtered one). Because the
  stream is CAI-descending, filtering preserves order.
- `max_examined` (default `max(200_000, 200*max_variants)`) bounds exploration so a locus where
  filters reject almost everything terminates (returning fewer than requested) instead of
  walking the whole lattice.
- `select_top_n_for_pridict(ranked, n)` clamps `n` to `[1000, 5000]` and slices.
- `split_insert(insert, end_len=40)` is a plain byte split with an exact-reconstruction
  guarantee (`first + middle + last == insert`); codon-frame handling for CAI-ranking the middle
  is applied where the middle is scored (Stage C), not in the split.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Full `itertools.product` enumeration then top-N (as `enumerate_variants` does) | Intractable: ~3^140 combos for a 500 bp middle; cannot run |
| Greedy per-position argmax (single max-CAI variant) only | Gives the #1 variant but not the top-N stream PRIDICT and the INSERT lever need |
| Approximate beam search | Unnecessary: exact k-best is cheap here (O(N*L) heap ops); no need to trade correctness |
| Fold filters into the search (constrained k-best) | More complex for no gain: lazy post-filtering on a CAI-descending stream is exact and simpler |

## Consequences

- Positive: top-N by CAI is exact and scales to realistic insert lengths without enumerating the
  product; the #1 variant is the global max-CAI; results are deterministic.
- Positive: reuses the existing CAI weights (`_CAI_REL_ADAPT`, `_CAI_EXCLUDED`, floor) and filter
  predicates, so ranking is consistent with `human_cai` and `enumerate_variants`' filters.
- Neutral: ordering among *exactly-equal-CAI* variants follows heap arrival, not the full
  secondary tiebreak, across CAI ties -- deterministic but not a clean lexicographic secondary
  order. The consumer (PRIDICT pre-filter) is indifferent to order among equal-CAI variants.
- Negative / documented: under very restrictive filters the result can contain fewer than
  `max_variants` (bounded by `max_examined`); `split_insert` is byte-aligned, so the middle's
  codon frame must be handled by the Stage-C caller that CAI-ranks it.

## References

- `src/pegrna_pipeline/hpc/codon_expander.py` (`_cai_choices`, `_cai_kbest`,
  `rank_insert_variants_by_cai`, `select_top_n_for_pridict`, `split_insert`).
- `tests/test_cai_prefilter.py` (k-best-vs-bruteforce completeness, scalability, determinism).
- Builds on ADR 0009 (human-CAI codon preference). Plan: Stage B in
  `docs/superpowers/plans/2026-07-02-unified-tjpe-target-strand.md`.
