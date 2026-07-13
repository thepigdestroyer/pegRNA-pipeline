# 0019. Importance gate maps the insertion footprint to ABSOLUTE protein residues

Date: 2026-07-13
Status: Accepted

## Context

The live CTCF api run (ADR 0018) reached the AlphaFold importance gate and blocked
every candidate on residues [47, 48, 49]. For a CTCF C-terminal insertion those are
the wrong residues: `site_s=147` in the small local design window maps to the
window-local residue index ~49, and `_insertion_site_footprint` returned that local
index when no `protein_length` was supplied. The gate then queried the FULL CTCF
AlphaFold structure (P49711, 727 aa) at residue 49 -- an N-terminal core -- instead
of the true C-terminal insertion site near residue ~695. A C-terminal knock-in was
silently judged as an N-terminal core: a false block (and, symmetrically, a possible
false pass).

`_insertion_site_footprint` already supported an absolute mapping via `protein_length`
(it anchors the window's LAST coding residue to the C-terminus, `offset =
(protein_length-1) - last_local_residue`), but the example config omitted
`protein_length`, so the mapping silently fell back to window-local indices. The fake
CI scorer is index-agnostic (constant metrics), so the bug only ever manifested on a
real AFDB run.

## Decision

- **The importance scorer is the authoritative source of protein length.**
  `AFDBImportanceScorer.protein_length()` returns the residue count of the loaded AFDB
  model -- the exact structure `metrics()` indexes.
- **`run_recoded_tjpe` auto-derives the length** via `_resolve_protein_length`: an
  explicit `protein_length` wins, else the scorer's `protein_length()` is used, else
  `None` (a frame-agnostic fake -> footprint stays window-local, harmless because such
  a scorer ignores the index). No config change is required for a real run to map
  local -> absolute correctly.
- The C-terminal anchor in `_insertion_site_footprint` is retained (the design is for
  C-terminal knock-ins; the window's last coding base is the C-terminus).

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Require the caller to pass `protein_length` | Silent mis-index if omitted (exactly what happened); pushes a correctness-critical value onto every config |
| Fail loud when `protein_length` is unknown on a real gate | Safer than the silent bug, but unnecessary: the scorer already knows the length authoritatively |
| Drop the C-terminal anchor, use a supplied CDS offset | Larger change; the pipeline targets C-terminal insertions, for which the anchor is exact |

## Consequences

- Positive: a real run maps the insertion footprint onto the correct absolute residues
  with no config burden, removing the silent local-vs-absolute mis-index.
- Positive: fake/CI behaviour is unchanged (no `protein_length()` -> window-local ->
  index-agnostic fake), so existing footprint tests are untouched.
- Caveat: for a multi-fragment AFDB entry only the F1 fragment is loaded, so
  `protein_length()` returns the fragment length for proteins longer than the fragment
  size (~2700 aa). Documented on the method; out of scope for the current single-fragment
  targets.
- Verification: re-run the CTCF api design and confirm the importance footprint reports
  C-terminal residues (~695), not N-terminal (~49).

## References

- `src/pegrna_pipeline/hpc/structure_importance.py` (`AFDBImportanceScorer.protein_length`);
  `src/pegrna_pipeline/hpc/recoded_tjpe_workflow.py` (`_resolve_protein_length`,
  `_insertion_site_footprint`).
- `tests/test_recoded_tjpe_workflow.py`. Builds on ADR 0013 (importance gate targets the
  insertion site) and ADR 0018 (the run that surfaced this).
