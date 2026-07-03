# Unified PRIDICT2-driven, target-strand-centric tjPE workflow — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development or
> superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.
>
> **STATUS:** design being finalised with the user (merge-feature table under revision). Only **Stage A**
> (target-strand convention + rcPBS2, self-contained) is cleared to start. Stages B+ await sign-off on the
> two-pegRNA/CAI mechanics below.

**Goal:** ONE tjPE workflow (merge `recoded_tjpe_workflow` into `tjpe_workflow`; **no separate recoded
module**) that uses PRIDICT2, under the target-strand (TS) convention, with the recode-arm + full gate
stack + verification, and an iterative search that returns the best gate-passing pegRNA.

**Convention (canonical everywhere):** the whole pegRNA 3′ extension = **RC(target strand)**. `rcPBS2 =
RC(PBS2) = RC(TS)` (PBS2 lives on the synthesised flap, not the pegRNA; Zheng 2023). Insert as it reads on
the TS = `coding_insert` if TS is the coding strand, else `RC(coding_insert)`; pegRNA carries `RC(that)`.

---

## The unified design (user-authoritative, 2026-07-02)

**Insert / CAI (CAI is NOT a mere tiebreak):**
1. CAI **enumerates and ranks all codon variants** of the insert by similarity to human-genome codon usage.
2. PRIDICT2 cannot process every variant → the pipeline selects the **top 1000–5000 CAI-ranked**
   (human-codon-similar) candidates and passes those to PRIDICT2, for the **two 40 bp ends** and for
   **both** pegRNAs (the PBS1-containing and the PBS2-containing designs).
3. The **middle** of the insert (beyond the two 40 bp ends) is enumerated + ranked by **CAI only**
   (PRIDICT never sees it).

**Two-pegRNA PRIDICT design + ngRNA:**
4. PRIDICT designs **pegRNA#1** → `PBS1` + first-40 bp of insert; PRIDICT already returns these as RC(TS)
   → **use untouched**.
5. PRIDICT designs a **separate pegRNA#2** for the PBS2 end → `rcPBS2 = RC(pegRNA#2 PBS)` and
   `last-40 bp = RC(pegRNA#2 insert)`, both fed back into the tj-pegRNA (RC'd because we target the
   opposite strand to what PRIDICT assumes).
6. **The ngRNA = pegRNA#2's spacer + scaffold** (that pegRNA's spacer targets the 2nd nick = the ngRNA
   nick). The ngRNA is **not** separately selected or genome-enumerated.

**Spacer QC (replaces `_spacer_passes_sanity`):**
7. Filter out only spacers **PRIDICT2 flags as high off-target risk** (verify PRIDICT2 exposes this).

**Between-nicks arm:**
8. **Recode-arm (Mechanism 2, `_partition_arm`) is canonical for CDS** spans (synonymous recode, human-CAI,
   `run ≤ max_run` → breaks microhomology). For **UTR / non-coding it is byte-exact (locked, Decision 15)**
   — critical for the fwd-targeting het-allele-aware pegRNA (edits in the UTR, must preserve the het allele
   byte-for-byte). `_partition_arm` already does CDS-recode + non-coding-lock. The arm is EXACTLY
   `genome[min(S,nick):max(S,nick)]`.

**Assemble → gate → iterate:**
9. Assemble the full tj-pegRNA, then run **ALL** hard gates (structure RiNALMo/RNAformer, AlphaFold
   insertion-site importance, cryptic-splice MaxEntScan + SpliceAI advisory, Boltz docking, off-target) —
   **verify every gate is active (check thoroughly).**
10. **Iterative search:** on a gate failure, depending on which gate failed, either pick a **new insert
    middle** (a CAI variant not in the two PRIDICT-designed 40 bp ends) OR choose a **new pegRNA#1/pegRNA#2
    (ngRNA) pair** from PRIDICT, and restart — iterating until the best pegRNA that passes **all** gates.
11. Do all of the above for **both** the reverse-targeting and the het-specific (fwd, allele-aware)
    approaches, ideally **in parallel**.

---

## Stages (Stage A cleared; B+ await sign-off)

- **Stage A — TS convention + rcPBS2 (self-contained, no PRIDICT).** `orient_pbs2`→`orient_rcpbs2`=RC(TS);
  `_select_pbs2`→rcPBS2; rename `pbs2`→`rcpbs2` in result/report/CSV/verification; extend `simulate_tprt`
  to validate the rcPBS2 capture; ADR + docs. Ship as its own PR; sync origin+neuro+CX3. **[details in
  the previous revision of this file; unchanged.]**
- **Stage B — CAI enumerate/rank + top-1000–5000 selection** feeding PRIDICT for the two ends (both
  pegRNAs) + CAI-only ranking for the middle.
- **Stage C — PRIDICT2 two-pegRNA adapter** (fake + real): pegRNA#1 (PBS1+first-40, untouched) and
  pegRNA#2 (rcPBS2+last-40 RC'd; spacer+scaffold → ngRNA). Reuse `pridict_submitter`/`pridict_parser`.
- **Stage D — spacer QC via PRIDICT off-target-risk** (replace `_spacer_passes_sanity`).
- **Stage E — recode-arm canonical for CDS, byte-exact for UTR** (confirm `_partition_arm` matches; wire
  UTR-lock for the fwd het-aware pegRNA).
- **Stage F — assemble + full gate stack + iterative search** (gate→action mapping: middle-insert vs
  new PRIDICT pair).
- **Stage G — reverse + het-specific in parallel;** unified result + verification DOCX (rcPBS2 labels).
- **Stage H — merge cleanup:** delete the separate recoded module; migrate its tests; CX3 end-to-end;
  ADRs; PR → merge → sync origin+neuro+CX3.

## Gate → lever map (LOCKED)

Two levers: **NICK** (vary the ngRNA nick 40–70 bp = a new pegRNA#2/ngRNA PRIDICT pair; S is fixed) and
**INSERT** (swap the insert *middle* for a CAI variant; the two PRIDICT-designed 40 bp ends are fixed).

- **HARD STOP (fail loud, no iteration):** `_partition_arm` no synonymous CDS solution / arm crosses >1
  coding boundary; INV-0 violation; TPRT strand mismatch.
- **off-target** (spacer/PBS/ngRNA, or PRIDICT off-target risk) → **NICK**.
- **rcPBS2 selection fail** → **NICK**.
- **cryptic-splice (MaxEntScan) + SpliceAI advisory** → site in the **cassette/insert** → **INSERT**;
  site in the **recode arm** → **NICK** (no separate re-recode lever).
- **structure** (RiNALMo/RNAformer cross-region pairing) → localized to the recode-arm portion → **NICK**;
  else → **INSERT**.
- **AlphaFold insertion-site importance** → **INSERT** (linker/cargo) first, then **HARD STOP** if no
  insert variant passes. **Implies an ADR-0013 upgrade:** the gate must model the *cargo's* structural
  impact at the fixed S (fold with the tag/linker), not just bare-site residue burial.
- **docking** (Boltz-2 ipTM) → **NICK** (new PRIDICT pair) first, then **INSERT**.
- **Winner:** highest PRIDICT2 PE2 among candidates that pass **all** gates. CAI top-N pre-filter size
  (1000–5000, default ~2000) configurable.

## Open items to verify DURING build (not blocking the plan)
1. PRIDICT2 exposes a per-spacer off-target-risk field (vs needing Cas-OFFinder) — verify on CX3.
2. Construct order (rcPBS2 5′-of-RTT vs 3′-terminal per Zheng) — confirm at Stage A Task A3.
3. `ts_is_coding` per-region (CDS vs UTR) source — config flag vs MANE strand inference.
4. PRIDICT2 field carrying the PBS-past-the-nick homology ("RTgDNA") — confirm on CX3 (Stage C).
