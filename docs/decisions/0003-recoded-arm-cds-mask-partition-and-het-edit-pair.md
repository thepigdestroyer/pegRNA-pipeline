# ADR 0003 - cds_mask-aware arm partition + sense-window het-edit pair

Status: accepted
Date: 2026-06-30
Context: recoded TJ-PE workflow (`run_recoded_tjpe`), driven by the CTCF C-terminal
knock-in design request (antisense reverse-targeting + sense-window het-edit).

## Problem

`run_recoded_tjpe` unconditionally recodes the whole re-templated arm
`genome[S:nick)` and requires it protein-silent (`recoded_tjpe_workflow.py`
arm block). It never partitions the arm by `cds_mask` and never passes
`locked_span`/`locked_expected` to `assert_inv0`. Consequences:

1. For an antisense window (arm = CDS) the behaviour is correct.
2. For a sense window where the arm 3' of S is **3'UTR** (the C-terminal
   het-edit case), `recode_arm` edits the non-coding sequence (measured: 5 edits
   on a repetitive 3'UTR-like arm), **violating locked Decision 15** (non-coding
   re-templated byte-exact). The byte-exact protection in `assert_inv0`
   (`scarless.py:88` `locked_span`/`locked_expected`) is never engaged because the
   orchestrator does not partition the arm.
3. No sense-window (Fw) het-edit (tjpegRNA, ngRNA) pair is emitted at all - the
   allele machinery (`enumerate_nick_movements`, `build_restored_segment`) is
   implemented + unit-tested but not wired into the pool-walk.

## Decision

1. **Partition the arm by `cds_mask`** inside the pool-walk. Split `genome[lo:hi)`
   into maximal same-class runs; recode each CDS run (frame-aware, synonymous,
   microhomology <= max_run); keep each non-coding run **byte-exact**
   (`build_restored_segment`-equivalent slice). Reassemble the sense arm in order,
   then orient per strand. The realistic C-terminal arm is monotone (>=1 CDS run
   and/or 1 non-coding run, single boundary at S's stop), but the partition
   handles the general case.
2. Compute a single `recode_span` (union of CDS runs) and single `locked_span`
   (union of non-coding runs) in genome coords; pass both + `locked_expected =
   genome[locked_span]` to `assert_inv0`. Protein-silence check applies to CDS
   runs only.
3. A non-coding (locked) run raises a **microhomology advisory** warning (it
   cannot be recoded away) - never a hard gate.
4. **Het-edit pair**: a thin driver runs the primary pass and, when
   `allele_hunt` + a discriminating phased het lies on the opposite-strand
   (sense/3'UTR) side, runs a second pass over the reverse-complement window
   (`pegrna_strand` flipped) and surfaces it as `het_edit_pair`. Each pass yields
   a full (tjpegRNA, ngRNA) pair.

## Rejected alternatives

- *Recode the 3'UTR anyway* - violates Decision 15 (RBP/polyA/splice/miRNA
  base-sensitivity); rejected.
- *Move the nick 5' of S to keep the arm coding for a sense window* - the
  orchestrator forbids 5'-of-S placement (would duplicate the arm); rejected.
- *Emit only the primary pair and let the caller re-run for the het-edit* -
  leaves the locked-span path dormant and the byte-exact bug latent; rejected.

## Invariants preserved

- All-CDS arm (antisense primary) path is byte-for-byte unchanged
  (`locked_span=None` -> existing recode path).
- INV-0, TPRT strand self-check, splice/importance/structure/docking gates,
  fail-loud-on-empty, S-never-moves: all retained.
- 618-test green baseline must remain green.
