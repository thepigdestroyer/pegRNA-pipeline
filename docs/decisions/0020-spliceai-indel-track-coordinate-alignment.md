# 0020. SpliceAI indel delta aligns in track coordinates (advisory splice cross-check)

Date: 2026-07-13
Status: Accepted

## Context

ADR 0017 wired `_PaddedSpliceAIScorer` into api mode so the SpliceAI advisory
cross-check could run on the web host's short design windows (N-pad up to the
>=10 kb context). Live verification of the FLAG example (FIX for the CTCF cascade)
showed the padded scorer WAS invoked (TensorFlow loaded) but the gate still recorded
`splice: advisory skip`. Reproduced on the box:

```
ValueError: operands could not be broadcast together with shapes (264,) (240,)
```

Root cause: `SpliceAIScorer.delta`'s indel branch computed the common prefix/suffix on
the raw INPUT sequences and then indexed the SpliceAI score TRACKS with those indices.
The model trims its context from each end, so a track is `len(input) - 2*context` long
(240 / 264 for the ref / alt here) while the input (and, once `_PaddedSpliceAIScorer`
N-pads it, the input + 2*5000) is far longer. The prefix computed on the padded input
(~5120) sliced the length-240 / length-264 tracks into full, mismatched arrays and
broadcast-errored. tjPE edits are always insertions (indels), so the advisory SpliceAI
cross-check broadcast-errored on EVERY candidate and silently degraded to advisory skip.
The branch is marked CX3-only (`# pragma: no cover`) and had never been exercised with
real weights, so the bug was latent until the first real indel run.

## Decision

- **Align the indel prefix/suffix in TRACK coordinates.** `delta` derives the per-side
  trim from the observed track length (`trim = (len(input) - len(track)) // 2`), slices
  each input down to the region its track covers (`*_core`), and computes the common
  prefix/suffix and the ref/alt lengths on those cores. This is robust to the model
  emitting full-length or context-trimmed output and to the N-padding added for short
  windows.
- The equal-length (SNV/MNV) branch already worked in track coordinates and is unchanged.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Leave it; MaxEntScan is the primary splice gate | The advisory SpliceAI cross-check (the whole point of ADR 0017's padded scorer) would remain dead for every tjPE insertion |
| Assume the model always trims exactly `context` | Fragile; deriving the trim from the actual track length also covers full-length-output models |
| Pad to exactly `2*context` and special-case | Produces zero-length tracks and does not fix the coordinate mismatch |

## Consequences

- Positive: the SpliceAI advisory cross-check now actually runs for tjPE insertions
  (padded on the web host, real flanks on CX3) instead of always advisory-skipping.
- Positive: the indel coordinate logic is now unit-tested in CI WITHOUT TensorFlow by
  injecting synthetic tracks (reproduces the old broadcast error; asserts a clean delta
  and that a signal in the inserted region is reported as a gain).
- Safety unchanged: MaxEntScan remains the PRIMARY splice gate (ADR 0011); SpliceAI is
  advisory either way.
- Verification: re-run the FLAG api design and confirm `gates.splice` reports a real
  SpliceAI delta rather than "advisory skip".

## References

- `src/pegrna_pipeline/hpc/splice_screen.py` (`SpliceAIScorer.delta` indel branch).
- `tests/test_api_splice_scorer_wiring.py` (track-coordinate indel tests). Builds on
  ADR 0017 (padded api splice scorer) and ADR 0011 (MaxEntScan primary gate).
- DSI reproduce: `_PaddedSpliceAIScorer().delta` broadcast error (240 vs 264).
