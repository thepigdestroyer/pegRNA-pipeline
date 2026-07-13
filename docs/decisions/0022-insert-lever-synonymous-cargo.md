# 0022. INSERT lever: synonymous-cargo retry on insert-side gate failures (Stage F)

Date: 2026-07-13
Status: Accepted

## Context

The unified plan (Stage F, gate->lever map §5) specifies two iteration levers: NICK (vary
the ngRNA nick) and INSERT (swap the insert for a CAI-ranked variant). The NICK lever was
already wired (the pipeline advances through admissible PRIDICT nicks). The INSERT lever
was not: `run_recoded_tjpe` took a fixed `cargo` and never varied it. The CAI machinery for
it already existed and was tested but unwired (`rank_insert_variants_by_cai`,
`select_top_n_for_pridict`, `split_insert` -- Stage B, ADR 0009/0015).

A synonymous cargo swap changes the pegRNA/RTT SEQUENCE, so it can only help the gates a
sequence affects: cryptic-splice (MaxEntScan), the SpliceAI cross-check, and RNA structure
(cross-region pairing). It cannot help cargo-INDEPENDENT gates: importance (judged on the
genomic insertion-site residues, ADR 0013/0019), off-target (spacer/PBS/ngRNA), INV-0/TPRT,
arm partition, or nick placement.

## Decision

- **`run_recoded_tjpe` becomes a thin wrapper** over the renamed core
  `_run_recoded_tjpe_single` (the full gate cascade, unchanged). The wrapper retries the
  core with the next CAI-ranked synonymous cargo variant.
- **Retry only on an insert-side failure.** `NoValidPegRNA` now carries `.diagnostics`; the
  wrapper retries iff a diagnostic contains an insert-side marker (splice / cryptic /
  cross-region / RiNALMo / structure). A cargo-independent failure breaks immediately -- a
  synonymous swap cannot rescue it, and re-running the cascade would just re-fail the same
  gate N times.
- **Variants are whole-cargo and synonymous**, generated only when the cargo is a clean
  in-frame CDS (length multiple of 3, translatable, no internal stop; a trailing stop codon
  is preserved). Otherwise the lever no-ops (original cargo only). This is valid for the
  recoded path because PRIDICT designs only the spacer/nicks from the site, not the cargo,
  so a synonymous swap leaves `site_s` and the spacer unchanged.
- **Bounded** by `max_insert_variants` (default 6; 1 disables the lever). The original cargo
  is always attempt 0.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Vary only the insert MIDDLE between two PRIDICT-designed 40 bp ends (plan's literal INSERT lever) | That is the orphaned two-pegRNA (Stage C) design; the recoded path does not lock the cargo ends, and the middle is frame-shifted (partial codons) -- whole-cargo synonymous variation is simpler and correct here |
| Reindent the 320-line candidate loop inside a cargo loop | High regression risk on the deployed core; the wrapper preserves the core verbatim |
| Retry on every failure | Wastes attempts re-failing cargo-independent gates (importance/off-target); the diagnostics gate the retry |
| Enumerate variants for any cargo | A synonymous swap is undefined for a non-CDS / frame-broken / internal-stop cargo; guarded to clean CDS |

## Consequences

- Positive: the INSERT lever is live; a cargo that trips splice/cryptic/structure can be
  rescued by a synonymous variant without moving the edit site or the spacer. Stage B's CAI
  functions are now wired into the design flow.
- Positive: cargo-independent failures still fail fast (no wasted iterations), and the lever
  no-ops safely for non-CDS cargos, so behaviour is unchanged where a swap cannot help.
- Neutral: `run_recoded_tjpe` is now a wrapper; the no-drift guard (ADR 0021) was updated to
  introspect `_run_recoded_tjpe_single`.
- F2: an explicit test now asserts every gate scorer (splice, cryptic, importance, structure,
  docking, off-target) is invoked on a full pass -- guarding against a silently-skipped gate.

## References

- `src/pegrna_pipeline/hpc/recoded_tjpe_workflow.py` (`run_recoded_tjpe` wrapper,
  `_run_recoded_tjpe_single`, `_cargo_variants`, `_is_insert_side_failure`, `NoValidPegRNA`).
- `src/pegrna_pipeline/hpc/codon_expander.py` (`rank_insert_variants_by_cai`, Stage B).
- `tests/test_insert_lever.py`. Builds on ADR 0009/0015 (CAI) and the §5 gate->lever map.
