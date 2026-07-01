# Scarless tjPE — Design Spec

**Status:** Accepted (user-directed; adversarial review recorded in §10). 2026-06-09.

**Goal:** Insert 200–800 bp cargo via template-jumping PE with **zero endogenous
deletion** ("scarless"). Standard tjPE deletes the span between the pegRNA nick
(A) and the nicking-sgRNA nick (B); this design re-templates that span.

## 1. Construct
3′ extension (5′→3′): `insert(cargo) · restored[A→B] · PBS1 · RC-PBS2 · tevopreQ1`.
Net product = reference allele + cargo at site S, nothing deleted.

## 2. INV-0 (hard gate)
`predict_edited_allele()` must equal `reference[:S] + cargo + reference[S:]`,
except **at most one** base change, which must be inside the **nicking-guide PAM**
and **translationally silent**. Any other endogenous change, or any length
mismatch, is rejected. Implemented in `hpc/scarless.py::assert_inv0`.

## 3. Strand handling
- Forward pegRNA: pegRNA nick A anchors S; cargo in sense.
- Reverse pegRNA: nicking nick **B must equal S** (`require_b_at_s`); cargo is
  templated as the **reverse complement**; restored span is reverse-complemented
  (`build_restored_segment`). The predicted **sense** allele is identical
  regardless of synthesis strand.

## 4. Single silent PAM break
Byte-exact restoration regenerates the nicking-guide protospacer → re-cut risk.
Mitigation: one silent base change inside that PAM (`silent_pam_break`). If none
exists, emit the byte-exact design + HIGH re-cut warning (no fallback). No other
endogenous edit (seed, pegRNA-PAM, microhomology-break) is permitted.

## 5. Distance
Target ~40 bp nick-to-nick (window 36–44) to minimise the restored span. **Below
the 70 bp minimum shown to work in Zheng 2023** — advisory warning, validate.

## 6. PRIDICT2.0 — two-end decomposition (≤40 bp / in-distribution)
- End 1: enumerate synonymous variants of the first ~40 bp of cargo; PRIDICT-rank.
- End 2: enumerate distance-valid PBS2; hard thermo gates (len/GC/Tm/structure);
  PRIDICT PE2-surrogate rank at the ngRNA nick.
- Combine by `min()` (bottleneck), stitch, then full-pegRNA structure + INV-0.
PRIDICT never sees the template-jump/dual-PBS architecture → scores are
ranking proxies only, never absolute efficiency.

## 7. Structure pre-filter
ViennaRNA `RNAfold` (`hpc/rnafold.py`): per-component pre-screen + full-pegRNA
check (MFE + PBS/spacer accessibility). Thresholds must be calibrated on real
data (no published pegRNA ΔG cut-off). Optional dependency (`.[scarless]`).

## 8. Fail-loud (R4)
When no candidate passes all hard gates, raise `ScarlessNoValidDesign` with a
reason — never return a silent empty/garbage design.

## 9. Policy
"A-only": always force scarless, no fallback to a deleting design. Microhomology,
re-cut, OOD, and sub-70 bp distance are **advisory** warnings on `TJPEResult`.

## 10. Risk register (adversarial review — surfaced, consciously accepted)
Independent multi-model consensus **rejected this design as specified**; the user
elected to proceed. Risks remain advisory per the A-only policy.

| # | Sev | Risk | Evidence |
|---|---|---|---|
| R1 | CRITICAL | Microhomology trap — restored[A→B] is 100% homologous to the displaced strand → imperfect/partial edits. | Zheng 2023 lines 274–280: observed partial edits "due to the RTT being homologous to genomic sequence"; their fix was to AVOID microhomology + shorten templated span (the opposite of scarless). |
| R2 | CRITICAL | "Advisory-only" makes named key-risks non-blocking. | Consensus (Security, AI/ML, Contrarian). |
| R3 | HIGH | ~40 bp < 70 bp validated floor, stacked on untested scarless mechanism; close opposite-strand nicks raise DSB/translocation risk. | Below published minimum. |
| R4 | HIGH | No-fallback + strict gates ⇒ may return nothing. | Consensus. Mitigated: fail loud with reason (§8). |
| R5 | MEDIUM | pegRNA's own protospacer at nick A is also regenerated (only the nicking-guide PAM is disrupted). | Mechanism. |
| R6 | MEDIUM | PRIDICT MIN-ranking over thousands of OOD codon variants may sort noise. | Consensus (AI/ML). |

**Required before bench use:** wet-lab validation of (a) scarless precision vs a
deleting tjPE control, and (b) the sub-70 bp distance. See vault `50-Open-Questions`.

## Implementation
Plan: `docs/superpowers/plans/2026-06-09-scarless-tjpe.md`. Code: `hpc/scarless.py`,
`hpc/rnafold.py`, `hpc/tjpe_workflow.py`. Tests: `tests/test_scarless*.py`,
`tests/test_rnafold.py`.
