# 0017. api-mode splice gate uses the N-padding SpliceAI adapter

Date: 2026-07-13
Status: Accepted

## Context

A wiring audit of the deployed path (web app -> `design-recoded` -> `recoded_driver`
-> `make_scorers(mode="api")` -> `run_recoded_tjpe`) found a planned feature that
was defined but never reached.

`scorer_factory._PaddedSpliceAIScorer` exists specifically for the GPU-less web
host: the real `SpliceAIScorer.delta` fail-louds unless handed `>= 2 * context`
(>= 10 kb, i.e. >= 5 kb each side), but a local design window on the web host is
only hundreds of nt. `_PaddedSpliceAIScorer` N-pads both REF and ALT up to that
context (N is SpliceAI's own zero-signal symbol; the edit region is preserved
between the flanks) so the SpliceAI advisory cross-check can run there.

`_api_bundle` wired `splice_scorer=SpliceAIScorer()` (the raw adapter), leaving
`_PaddedSpliceAIScorer` orphaned. Consequence on the deployment: every candidate
hit the short-context `ValueError`, which `run_recoded_tjpe` catches and records
as `"advisory skip (needs >=10kb context; MaxEnt is primary)"`. The SpliceAI
cross-check therefore never fired on the web host. Not a safety hole (the
MaxEntScan cryptic-splice gate is the PRIMARY splice gate, ADR 0011, and always
runs), but a planned advisory gate silently did nothing.

## Decision

- **`_api_bundle` uses `_PaddedSpliceAIScorer()`** for `splice_scorer`. The unused
  `from ...splice_screen import SpliceAIScorer` is dropped from `_api_bundle`
  (`_PaddedSpliceAIScorer` is module-level and imports the raw scorer itself).
- **`_real_bundle` keeps the raw `SpliceAIScorer()`.** The HPC path supplies the
  true GRCh38 flanks (>= 5 kb each side), so padding is neither needed nor correct
  there (padding real context would dilute the signal).
- **The shared `run_recoded_tjpe` exception handling is left untouched.** The
  advisory-skip branch there catches `ValueError` (short context) only; a genuine
  scorer crash (e.g. missing weights) still skips the candidate. Broadening it to
  swallow `ImportError` would violate real mode's fail-loud contract ("no silent
  fallback in real mode"), so robustness for a mis-provisioned web host lives in
  the scorer selection, not the shared loop. On the actual deployment the
  `pegrna-real` env has TensorFlow + the SpliceAI weights (spliceai needs
  `setuptools<81` for `pkg_resources`), so the padded cross-check runs for real.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Leave `_api_bundle` on the raw `SpliceAIScorer` | The advisory SpliceAI cross-check never runs on the web host; a planned gate is dead code |
| Broaden `run_recoded_tjpe` to treat `ImportError` as advisory-skip too | Shared by real mode; would break real mode's fail-loud "no silent fallback" contract |
| Use the padded adapter in `_real_bundle` as well | HPC supplies true >= 5 kb flanks; N-padding real context dilutes the delta and hides real splice signal |
| Add an `inner=` param to `_PaddedSpliceAIScorer` for testability | Unnecessary; tests inject the inner by patching `splice_screen.SpliceAIScorer` (the constructor's own import target), leaving production untouched |

## Consequences

- Positive: on the web host the SpliceAI cross-check now runs (N-padded) instead
  of being silently skipped, so the api-mode splice gate matches its design.
- Positive: real mode is provably unchanged (regression test asserts
  `type(make_scorers(mode="real").splice_scorer) is SpliceAIScorer`).
- Positive: the padding math is unit-tested in CI without TensorFlow by injecting a
  recording inner (short window -> padded to >= 10 kb, N-flanked, edit preserved;
  adequate window -> passed through verbatim).
- Follow-up (out of scope here): the legacy `run` / `tjpe_workflow.assemble_tj_pegrna`
  path is NOT on the deployed critical path (the web app shells to `design-recoded`
  only; the legacy path is HPC-only and still carries the pre-`rcPBS2` orientation).
  It is left uncorrected and flagged for a separate change; the recoded path already
  derives `rcPBS2` correctly via `tjpe_workflow.derive_rc_pbs2` (ADR 0014).
- Verification: confirm on a live DSI api run that `gates.splice` reports a real
  SpliceAI delta rather than `"advisory skip"`.

## References

- `src/pegrna_pipeline/hpc/scorer_factory.py` (`_PaddedSpliceAIScorer`, `_api_bundle`,
  `_real_bundle`).
- `src/pegrna_pipeline/hpc/splice_screen.py` (`SpliceAIScorer.delta` 10 kb-context
  fail-loud); `src/pegrna_pipeline/hpc/recoded_tjpe_workflow.py` (advisory-skip branch,
  ~L766-796).
- `tests/test_api_splice_scorer_wiring.py`; builds on ADR 0011 (MaxEntScan primary
  splice gate) and ADR 0008 (real scorer resources).
