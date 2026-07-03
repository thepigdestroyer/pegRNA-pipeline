# 0014. Target-strand convention: the pegRNA carries rcPBS2 = RC(TS), 5′-of-RTT

Date: 2026-07-02
Status: Accepted

## Context

The recoded TJ-PE path emitted the **second-strand primer element in the wrong orientation**.
`strand_orientation.orient_pbs2`, for a reverse-targeting (`Rv`) pegRNA, returned the captured
target-strand window **literally** ("complement-of-complement = identity"). That is PBS2, not
rcPBS2. In template-jumping prime editing (tjPE) PBS2 does **not** sit on the pegRNA — it is
synthesised onto the 3′ end of the new DNA flap and captures the flap the nicking guide (ngRNA)
exposes on the non-target strand. The pegRNA carries the **reverse complement**, `rcPBS2 =
RC(PBS2)`. Because `simulate_tprt` only validated the RTT tail (not the rcPBS2), the wrong
orientation was never caught at runtime.

A second, separate confusion existed about **where** rcPBS2 sits in the pegRNA oligo.
`report_writer.py` and the older `assemble_tj_pegrna` claimed `[spacer][scaffold][insert][PBS1]
[RC-PBS2]` — RC-PBS2 3′-terminal — "per Zheng 2023 Fig 1a". This is a misreading: the paper's
prose lists the 3′-extension components as "insertion, PBS1, RC-PBS2", which is a **component
list, not a 5′→3′ layout**. Read as an order it wrongly puts RC-PBS2 at the 3′ terminus.

Two coordinate frames must not be conflated:
- **On the target strand / the edit:** PBS2's capture site is 3′ of PBS1's site.
- **On the pegRNA oligo:** the layout is `5′-spacer·scaffold·rcPBS2·RTT·PBS1·PolyU-3′` — rcPBS2
  is **5′-of-RTT** (immediately after the scaffold); PBS1 is the 3′-most *templating* element
  (it must be, to prime off the genomic-nick 3′-OH), and a PolyU / tevopreQ1 motif caps the 3′
  end. Mechanistically rcPBS2 **cannot** be 3′ of PBS1: RT extends away from PBS1 toward the
  5′ end, so anything 3′ of PBS1 is never copied into the flap and PBS1 could not prime.

## Decision

Adopt the **target-strand (TS) convention**: every pegRNA 3′-extension component is `RC(TS)` at
its locus (`PBS1`, `replacement_gDNA`, `insert`, `rcPBS2`). Concretely for the rcPBS2 element:

- `strand_orientation.orient_rcpbs2(ts_window)` returns `reverse_complement(ts_window)` (no
  strand branch; the caller supplies the target-strand window). `orient_pbs2` becomes a
  deprecated alias that warns and delegates (behaviour changed: it now returns `RC(TS)`).
- `_select_rcpbs2` (renamed from `_select_pbs2`) selects the length to hit the Tm band and emits
  `rcPBS2 = RC(TS)` (Rv: `RC(sense window)`; Fw: the sense window literally). Tm(NN), GC, and the
  repetitiveness screen are reverse-complement-invariant, so the length choice is unchanged by
  the orientation fix — only the emitted bases change. It also rejects non-ACGT windows
  (assembly-gap `N` runs), which `Tm_NN` silently ignores.
- `simulate_tprt` gains `simulate_rcpbs2_capture(rcpbs2, ts_capture_window)`: the flap-borne
  `PBS2 = RC(rcpbs2)` must recover the target-strand capture window; a TS-literal rcPBS2 raises
  `TPRTMismatch` (locus-dependent — blind on RC-palindromic windows).
- The public `pbs2` name is renamed to `rcpbs2` across `RecodedTJPEResult`, `_regions_for`,
  `passes_pbs2_gate`→`passes_rcpbs2_gate`, and the text/CSV/DOCX report + verification legend.

**Construct order is confirmed 5′-of-RTT** — `Spacer·Scaffold·rcPBS2·RTT·PBS1·PolyU` — which is
what `assemble_tj_pegrna_v2` already produces; the recoded path's assembly order was correct.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| Keep `orient_pbs2` = TS-literal for Rv (identity) | Biologically wrong: the pegRNA carries rcPBS2 = RC(TS), not PBS2; PBS2 lives on the flap (Zheng 2023) |
| Place rcPBS2 3′-terminal (per `report_writer.py:43` / `assemble_tj_pegrna`) | Mechanistically impossible: PBS1 must be the 3′-most templating element to prime; RT never copies anything 3′ of PBS1. The "3′-terminal" claim is a misreading of the paper's prose component-list |
| Wire the `simulate_rcpbs2_capture` guard into the candidate loop | It re-derives the TS window with the same Rv/Fw branch `_select_rcpbs2` uses, so `RC(RC(x))==x` is a near-tautology (blind to a branch-swap and to palindromes). Kept as a validated utility; a real independent check belongs in the PRIDICT stage where `Target-Strand` is supplied explicitly |
| Rename the shared `assemble_tj_pegrna_v2(pbs2=…)` param too | Out of scope; the assembler is a generic 6-slot primitive unified in the Stage-H merge. Renaming it now would churn the tjpe/scarless path + tests for no correctness gain |

## Consequences

- Positive: the recoded pegRNA now carries the correct `rcPBS2 = RC(TS)`; a mis-orientation
  fails loud via `simulate_rcpbs2_capture`; assembly-gap `N` windows can no longer produce a
  non-priming rcPBS2.
- Positive: naming is coherent (`rcpbs2` throughout result/report/CSV/verification), and the
  verification legend/report describe the correct 5′-of-RTT layout.
- Negative / follow-up: the **tjpe/scarless path is still wrong** — `assemble_tj_pegrna`
  assembles RC-PBS2 3′-terminal and the shared `_CAVEATS` in `report_writer.py` still documents
  that 3′-terminal architecture. These are consistent with each other (both wrong) and are left
  untouched here to avoid a tjpe report/code mismatch; **they must be corrected together in the
  unified-workflow merge (Stage H)** when `run_tjpe_workflow` adopts the 5′-of-RTT assembly.
- Neutral: `orient_pbs2` remains as a deprecated shim for any external caller.

## References

- Zheng C, Liu B, Dong X, Gaston N, Sontheimer EJ, Xue W. "Template-jumping prime editing
  enables large insertion and exon rewriting in vivo." *Nat Commun* 2023;14:3369
  (PMID 37291100). RC-PBS2 is 5′ of the 3′-extension; the prose component-list is not a 5′→3′
  order. Cross-check: Anzalone 2019 (*Nature* 576:149) canonical pegRNA (PBS 3′-terminal).
- `src/pegrna_pipeline/hpc/strand_orientation.py` (`orient_rcpbs2`, deprecated `orient_pbs2`);
  `tprt_sim.py` (`simulate_rcpbs2_capture`); `recoded_tjpe_workflow.py` (`_select_rcpbs2`,
  `RecodedTJPEResult.rcpbs2`, `passes_rcpbs2_gate`, `_regions_for`); `verification.py`,
  `report_writer.py`.
- Known-wrong tjpe-path order to fix in Stage H: `report_writer.py` `_CAVEATS`,
  `tjpe_workflow.assemble_tj_pegrna`.
