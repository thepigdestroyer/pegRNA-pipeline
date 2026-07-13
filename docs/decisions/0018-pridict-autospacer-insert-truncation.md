# 0018. PRIDICT auto-spacer editseq truncates the cargo to 40 bp

Date: 2026-07-13
Status: Accepted

## Context

Live api-mode verification on the DSI box (the real CTCF C-terminal
`config_autospacer.json`, a 503 bp cargo) failed: the auto-spacer PRIDICT run
produced an empty output directory and the pipeline died with a confusing
`FileNotFoundError` on the missing `peg1_pegRNA_Pridict_full.csv`.

A controlled probe on the box established the cause. PRIDICT2 hard-rejects large
insertions:

```
The chosen insertion is 100bp long. Only insertions up to 40bp are currently
supported by default! ... Invalid insertion length (100).
```

PRIDICT skips the offending row but still exits 0, so no CSV is written and the
subprocess "succeeds". The 503 bp cargo was far over the 40 bp ceiling.

Crucially, in the recoded tj-PE path PRIDICT is used ONLY as a spacer + nicking-guide
designer: `_design_peg1_via_pridict` returns just `.spacer` and the nicking rows and
discards PRIDICT's RTT/PBS. Spacer and nick placement are set by the genomic PAM
landscape around `site_s`, not by the insert length, so PRIDICT never needed the full
cargo. Feeding it the full cargo was the bug.

## Decision

- **Truncate the cargo to its leading `_PRIDICT_MAX_INSERTION = 40` bp when building the
  PRIDICT auto-spacer editseq** (`recoded_driver._design_peg1_via_pridict`). The FULL
  cargo still flows unchanged into `run_recoded_tjpe` for the actual construct; only the
  PRIDICT design editseq is truncated. A cargo <= 40 bp is passed whole.
- **`PridictLocalRunner._collect` fails loud with an actionable message** when a pegRNA
  CSV is absent after a zero-exit run: it names the likely causes (insertion > 40 bp, or
  < 100 bp flanking each side) instead of surfacing a bare `FileNotFoundError`.
- **Deploy-config fix (same session, no code):** the deployed
  `~/.config/pegrna/pegrna.env` predated the PR and lacked `PRIDICT2_DIR/PYTHON/CORES`
  (`deploy.sh` only copies the env file when absent), so api-mode auto-spacer had no
  runner at all. The block was appended to the live env and the service restarted.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Reject cargos > 40 bp | Defeats the feature; tj-PE exists precisely for large insertions |
| Feed PRIDICT a fixed dummy insert | The cargo's real leading window is faithful and free; a dummy adds a needless divergence |
| Raise PRIDICT2's 40 bp limit | Not a user-facing knob; the ceiling is a property of the released model |
| Only fix the runner error message | Leaves the auto-spacer path non-functional for any cargo > 40 bp |

## Consequences

- Positive: api-mode auto-spacer now works for arbitrary cargo sizes; the spacer/nick
  guides are identical to what a full-length insertion would yield.
- Positive: a skipped/oversized PRIDICT edit now fails with a clear, actionable error.
- Positive: unit-tested without launching PRIDICT (recording runner asserts the editseq
  carries at most 40 bp; a <=40 bp cargo is passed whole; `_collect` raises on a missing
  CSV).
- Verification: re-run the real CTCF api design on the DSI box and confirm a spacer is
  auto-designed, the ngRNA comes from PRIDICT nicking guides, and the gate cascade runs.

## References

- `src/pegrna_pipeline/hpc/recoded_driver.py` (`_PRIDICT_MAX_INSERTION`,
  `_design_peg1_via_pridict`); `src/pegrna_pipeline/hpc/pridict_local.py` (`_collect`).
- `tests/test_autospacer_wiring.py` (truncation), `tests/test_pridict_local.py`
  (missing-CSV error). Builds on ADR 0016 (two-peg PRIDICT adapter) and the 2026-07-10
  wiring spec.
- DSI probe: `PRIDICT2/predictions/probe_1783928961/probe.log` (40 bp ceiling).
