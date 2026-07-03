# Unified PRIDICT2-driven, target-strand-centric tjPE workflow — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development or
> superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax. Read the whole plan before
> starting a stage; each stage ends with a checkpoint that must pass before the next stage begins.
>
> **STATUS (2026-07-02):** design finalised and locked with the user. **Stage A is cleared to start.**
> Stages B–H are fully specified here; a few PRIDICT2-schema/CX3 facts are flagged to verify *during*
> the relevant stage (they do not block planning).

**Goal:** Replace the two divergent tjPE code paths with **one** tjPE workflow that (a) is
target-strand-centric, (b) sources PBS1 / insert-ends / rcPBS2 / the ngRNA from **PRIDICT2**, (c) recodes
the between-nicks arm to break microhomology (CDS) or keeps it byte-exact (UTR), (d) runs the full hard-gate
stack on the assembled tj-pegRNA, and (e) iterates intelligently (two levers) until it returns the best
gate-passing pegRNA + its paired ngRNA — for both the reverse-targeting and the het-specific forward designs.

**Architecture:** The whole pegRNA 3′ extension is the reverse complement of the **target strand (TS)**.
CAI enumerates/ranks all insert codon variants and pre-filters the top 1000–5000 to PRIDICT2, which designs
two pegRNAs (the PBS1 end and the PBS2 end); the second pegRNA's spacer+scaffold *is* the ngRNA. The arm
between the two nicks is recoded (CDS) or locked (UTR). Assembled candidates run every gate; failures are
routed to one of two levers (move the nick, or change the insert middle) by a locked gate→lever map.

**Tech Stack:** Python 3.11, Biopython, `python_codon_tables`, PRIDICT2 (CX3; conda env `pridict2`;
`pridict2_pegRNA_design.py batch`), DnaChisel, maxentpy (MaxEntScan), SpliceAI, RiNALMo/RNAformer, Boltz-2,
cas-offinder, python-docx, pytest, ruff, black (line-length 100).

---

## 1. Background & motivation

Two tjPE implementations diverged in the repo:

- **`tjpe_workflow.py`** — the older, **PRIDICT2-driven** path: ranks codon variants by PRIDICT PE2,
  does a two-ended PRIDICT decomposition, derives `rc_pbs2` correctly (`derive_rc_pbs2`, 3′-terminal per
  Zheng 2023), selects an ngRNA (`select_tjpe_ngrna`), and has a scarless mode
  (`build_restored_segment` — restores the between-nicks span **verbatim**, leaving microhomology intact).
- **`recoded_tjpe_workflow.py`** — the newer, **genome-derived** path (all our CTCF work): recodes the
  between-nicks arm to **break microhomology** (`_partition_arm`), runs a full gate stack
  (off-target / cryptic-splice / importance / structure / docking), emits the verification DOCX, and
  uses human-CAI codon preference — but it **does not use PRIDICT2** (it derives PBS1/PBS2 from the
  genome) and it carried a **PBS2 orientation bug** (put PBS2 where rcPBS2 belongs).

The user's decision: **there must be ONE workflow.** Merge the recoded path's recode-arm + gate-stack +
verification into the PRIDICT2-driven path, governed by a single target-strand convention. This plan is
that merge.

The CTCF real-mode gate cascade that motivated the correctness work (all fixed and merged):
genome-wiring (PR #11) → off-target near-tier counting (PR #13) → splice via MaxEntScan/maxentpy
(PR #14/#15) → insertion-site importance (PR #17) → verification DOCX (PR #16).

---

## 2. The target-strand (TS) convention (canonical — apply everywhere)

- **TS = the strand the TJ-pegRNA targets** (NOT the ngRNA). It may be the coding or the non-coding strand.
- **Every pegRNA 3′-extension component = RC(TS)** at its locus:
  - `PBS1 = RC(TS)` at the pegRNA nick (S).
  - `replacement_gDNA = RC(TS over the between-nicks span)`.
  - `insert = RC(desired new TS)`.
  - `rcPBS2 = RC(TS at the capture locus)`.
- **rcPBS2, not PBS2.** PBS2 exists on the newly-synthesised DNA flap (it captures the flap the ngRNA
  generates on the non-target strand). The tj-pegRNA carries the reverse complement, `rcPBS2 = RC(PBS2)
  = RC(TS)` (Zheng 2023 naming).
- **Insert coding rule.** The insert *as it reads on the TS* = `coding_insert` if TS is the coding strand,
  else `RC(coding_insert)`. The pegRNA then carries `RC(that)` (like every component). Net pegRNA insert =
  `RC(coding_insert)` if TS=coding, else `coding_insert`.

---

## 3. The unified workflow (end-to-end)

1. **CAI enumerate + rank.** Enumerate all synonymous codon variants of the insert; rank by CAI
   (similarity to the human-genome codon table). CAI is a first-class ranker, **not** a tiebreak.
2. **CAI pre-filter → PRIDICT.** PRIDICT2 cannot process every variant; hand it the **top 1000–5000**
   CAI-ranked candidates (configurable) — for the **two 40 bp ends** and for **both** pegRNAs (PBS1-end
   and PBS2-end designs). The insert **middle** (beyond the two ends) is ranked by **CAI only** (PRIDICT
   never sees it).
3. **Two-pegRNA PRIDICT design.**
   - **pegRNA#1** → `PBS1` + first-40 bp of insert; PRIDICT already returns these as RC(TS) → **use
     untouched**.
   - **pegRNA#2** → `rcPBS2 = RC(pegRNA#2.PBS)` and `last-40 bp = RC(pegRNA#2.insert)`, fed into the
     tj-pegRNA (RC'd because we target the opposite strand to what PRIDICT assumes). **The ngRNA =
     pegRNA#2.spacer + pegRNA#2.scaffold, UNTOUCHED** (only its PBS + insert are RC'd).
4. **ngRNA lifecycle.** The ngRNA is **stored bound** to its tj-pegRNA. If the tj-pegRNA passes all
   gates → presented as a **pair** (tj-pegRNA + ngRNA). If it fails → **both discarded/regenerated**.
5. **Spacer QC** = filter out spacers PRIDICT2 flags as **high off-target risk** (replaces the
   TTTT/homopolymer/GC heuristic).
6. **Between-nicks arm** = exactly `genome[min(S,nick):max(S,nick)]`. **CDS → recode** (`_partition_arm`,
   synonymous, human-CAI, `run ≤ max_run` → breaks microhomology). **UTR/non-coding → byte-exact/locked**
   (Decision 15) — critical for the fwd het-allele-aware pegRNA that edits in the UTR.
7. **Assemble the full tj-pegRNA** (`Spacer·Scaffold·rcPBS2·RTT·PBS1·PolyU`, RTT = insert +
   replacement_gDNA) and **run EVERY hard gate** (verify each is active).
8. **Iterative search** (two levers) until the best gate-passing pegRNA is found (see §5).
9. **Two designs in parallel:** reverse-targeting and het-specific forward (allele-aware).
10. **Winner** = highest PRIDICT2 PE2 among candidates passing **all** gates. Emit the verification DOCX
    for the winning pair.

---

## 4. Decision log (every decision + rationale)

| # | Decision | Rationale |
|---|---|---|
| D1 | ONE tjPE workflow; **no separate `recoded_tjpe`** | Two divergent paths caused the PBS2 bug + PRIDICT gap; a single source of truth prevents drift |
| D2 | The unified workflow **still uses PRIDICT2** | PE2 efficiency ranking + validated PBS/RT designs are the reason to use PRIDICT at all |
| D3 | Target-strand-centric convention (everything = RC(TS)) | Removes sense/antisense/coding ambiguity; one rule for all components |
| D4 | `rcPBS2 = RC(TS)` (PBS2 lives on the flap, not the pegRNA) | Zheng 2023 mechanism; the recoded path's `orient_pbs2`=TS-literal was a bug |
| D5 | Insert coding rule (coding_insert vs RC by TS) | The protein must read correctly whichever strand is targeted |
| D6 | CAI enumerates/ranks all variants; top 1000–5000 → PRIDICT; middle CAI-only | PRIDICT can't process every variant (compute); the middle isn't PRIDICT's concern; CAI keeps codons human-like |
| D7 | Two-pegRNA PRIDICT design; **ngRNA = pegRNA#2 spacer+scaffold untouched** | pegRNA#2's spacer already targets the 2nd nick; only its PBS+insert are RC'd into the tj-pegRNA |
| D8 | ngRNA stored **bound** to its tj-pegRNA; presented/discarded **together** | They are one physical design; a discarded tj-pegRNA invalidates its ngRNA |
| D9 | Spacer QC = PRIDICT off-target risk | PRIDICT's model is more informative than the TTTT/GC heuristic |
| D10 | recode-arm canonical for **CDS**; byte-exact for **UTR** | Breaking microhomology needs synonymous recode (CDS only); UTR must stay byte-exact (het-allele) |
| D11 | Winner = highest PRIDICT PE2 among gate-passers; CAI top-N configurable | PE2 is the efficiency signal; gates are hard constraints; CAI is the pre-filter |
| D12 | Two levers: **NICK** (40–70 bp) and **INSERT** (middle CAI variant); **S is user-fixed** | The pipeline may move the nick and change the insert middle, but never the insert location |
| D13 | Locked gate→lever map (see §5) | Deterministic, debuggable iteration |
| D14 | AlphaFold gate upgraded to model **cargo impact** at S (not bare-site burial) | So the INSERT lever can actually help before hard-stop |
| D15 | Correctness invariants (`_partition_arm`, INV-0, TPRT) = **hard stop** | These shouldn't fail for a well-formed locus; surface, don't silently iterate |
| D16 | "Sensitive nucleotides" = union of existing gates over the nick/recode-arm region | No new annotation needed; the gates already define sensitivity |
| D17 | Run reverse-targeting + het-forward designs **in parallel** | Independent designs; parallelism saves wall-clock |

---

## 5. Gate → lever map (LOCKED)

**Two levers.** **NICK:** vary the ngRNA nick within 40–70 bp = pick a new pegRNA#2/ngRNA PRIDICT pair
(S is fixed). **INSERT:** swap the insert *middle* for a CAI variant (the two PRIDICT-designed 40 bp ends
are fixed).

| Gate | Failure meaning | Lever |
|---|---|---|
| `_partition_arm` no synonymous CDS solution / arm crosses >1 coding boundary | invariant | **HARD STOP** |
| INV-0 (`assert_inv0`) violation | invariant | **HARD STOP** |
| TPRT (`simulate_tprt`) strand mismatch | invariant | **HARD STOP** |
| Off-target (spacer/PBS/ngRNA; PRIDICT off-target risk) | placement | **NICK** (new PRIDICT pair) |
| rcPBS2 selection fail (no Tm/GC/len-valid) | placement | **NICK** |
| Cryptic-splice (MaxEntScan) + SpliceAI advisory — site in **cassette/insert** | molecule | **INSERT** |
| Cryptic-splice + SpliceAI advisory — site in **recode arm** | placement | **NICK** |
| Structure (RiNALMo/RNAformer cross-region pairing) — in **recode arm** | placement | **NICK** |
| Structure — elsewhere in the pegRNA | molecule | **INSERT** |
| AlphaFold insertion-site importance | cargo impact | **INSERT** first (linker/cargo) → **HARD STOP** if none pass |
| Docking (Boltz-2 ipTM floor/drop) | complex | **NICK** first → **INSERT** fallback |

**Winner:** highest PRIDICT2 PE2 among candidates passing **all** gates. **CAI top-N** pre-filter size
(1000–5000, default ~2000) is a config knob.

**Iteration control (implement as explicit, capped loops):** an outer loop over PRIDICT pegRNA#1/#2 pairs
(NICK lever) and an inner loop over CAI insert-middle variants (INSERT lever). Each lever has a max-attempt
cap; on exhaustion, emit a precise diagnostic naming the gate + lever tried. Never silently return nothing.

---

## 6. File structure

- `src/pegrna_pipeline/hpc/strand_orientation.py` — TS-centric; `orient_pbs2`→`orient_rcpbs2` = RC(TS).
- `src/pegrna_pipeline/hpc/tjpe_workflow.py` — becomes the single canonical workflow (absorbs the recoded
  features); `recoded_tjpe_workflow.py` is deleted at Stage H after migration.
- `src/pegrna_pipeline/hpc/recoded_pridict.py` (new) — two-pegRNA PRIDICT adapter (fake+real), CAI
  pre-filter, ngRNA lifecycle, RTgDNA handling.
- `src/pegrna_pipeline/hpc/codon_expander.py` — CAI enumerate/rank + top-N selection (extend existing).
- `src/pegrna_pipeline/hpc/structure_importance.py` — upgrade the importance gate to model cargo impact.
- `src/pegrna_pipeline/hpc/verification.py`, `report_writer.py` — rcPBS2 labels; targeting flag.
- `src/pegrna_pipeline/hpc/scorer_factory.py`, `recoded_driver.py` — adapter wiring; parallel reverse+het.
- `src/pegrna_pipeline/hpc/tprt_sim.py` — validate rcPBS2 capture.
- `examples/ctcf_cterm/config.json` — `coding_insert` + per-region `ts_is_coding`.
- `docs/decisions/0014..NNNN-*.md`, `docs/FEATURES.md`, `docs/DESIGN_HISTORY.md`.

---

## 7. Stage A — Target-strand convention + rcPBS2 (self-contained; NO PRIDICT)

**Checkpoint A (exit criteria):** the pegRNA 3′ extension is RC(TS) end-to-end incl. `rcPBS2`; the field is
named `rcpbs2` in the result/report/CSV/verification; `simulate_tprt` validates the rcPBS2 capture; full
suite green; ruff+black clean; ADR 0014 written; shipped as its own PR; origin+neuro+CX3 synced.

### Task A1 — `orient_rcpbs2` returns RC(TS)
**Files:** Modify `strand_orientation.py`; Test `tests/test_strand_orientation.py`.
- [ ] Failing test: `orient_rcpbs2("ACGTACGTACGT")` == reverse-complement of that window.
- [ ] Run → fail (undefined). Implement `def orient_rcpbs2(ts_window): return _rc(ts_window)`; keep
  `orient_pbs2` as a deprecated alias that emits `DeprecationWarning` and returns `orient_rcpbs2(...)`.
- [ ] Run → pass. Commit `refactor(strand): rcPBS2 = RC(target strand)`.

### Task A2 — `_select_pbs2` → `_select_rcpbs2` uses `orient_rcpbs2`
**Files:** `recoded_tjpe_workflow.py` (until merged) / `tjpe_workflow.py`; Test `tests/test_pbs2_selection.py`.
- [ ] Update the pinning test to expect `RC(captured)` (was the literal window for Rv); show one worked value.
- [ ] Run → fail. Rename + call `orient_rcpbs2`; Tm/GC/repeat now run on the RC. Run → pass. Commit.

### Task A3 — rename `pbs2`→`rcpbs2` across result/report/CSV/verification
**Files:** `recoded_tjpe_workflow.py` (`RecodedTJPEResult`, `_regions_for`), `report_writer.py`,
`verification.py` (labels + CATEGORY_STYLE), and their tests.
- [ ] Rename field + all references; update tests asserting `pbs2`→`rcpbs2`. Commit.
- [ ] **CHECKPOINT DECISION (confirm with user):** construct order — rcPBS2 5′-of-RTT (current recoded,
  "locked B1") vs 3′-terminal (`spacer·scaffold·insert·PBS1·rcPBS2`, Zheng 2023 + `tjpe_workflow`). Record
  the choice in ADR 0014.

### Task A4 — extend `simulate_tprt` to validate rcPBS2 capture
**Files:** `tprt_sim.py`; Test `tests/test_tprt_rcpbs2.py`.
- [ ] Failing test: a correct rcPBS2 reconstructs the 2nd-strand capture; a non-RC'd rcPBS2 raises
  `TPRTMismatch`. Implement the capture reconstruction (PBS2-on-flap = RC(rcPBS2) must equal the TS window
  it captures). Run → pass. Commit.

### Task A5 — ADR 0014 + docs
- [ ] Write `docs/decisions/0014-target-strand-convention-and-rcpbs2.md` (the §2 rules, the PBS2-vs-rcPBS2
  clarification, rejected alternatives). Update `FEATURES.md` + verification legend. Commit.
- [ ] **Ship Stage A as its own PR → merge → sync origin + neuro + CX3.** Full suite + lint gates green.

---

## 8. Stage B — CAI enumerate/rank + top-N pre-filter

**Checkpoint B:** CAI ranks all insert variants by human-codon similarity; the top-N (config, 1000–5000)
are selectable for PRIDICT; the insert middle is CAI-ranked independently; unit-tested; no PRIDICT yet.

- [ ] **B1** — `codon_expander`: `rank_insert_variants_by_cai(insert_aa|insert_dna) -> ranked list`
  (reuse `human_cai`). TDD: ranking is by descending CAI; deterministic tiebreak.
- [ ] **B2** — `select_top_n_for_pridict(ranked, n)` with `n` clamped to [1000, 5000] (default 2000). TDD.
- [ ] **B3** — split the insert into `first_end (≤40)`, `middle`, `last_end (≤40)`; TDD reconstruction
  (`first+middle+last == insert`). The middle is CAI-ranked; the ends go to PRIDICT.

---

## 9. Stage C — PRIDICT2 two-pegRNA adapter (fake + real) + ngRNA lifecycle

**Checkpoint C:** a fake adapter drives the full pipeline in CI; the real adapter submits PRIDICT on CX3;
pegRNA#1 ends used untouched; pegRNA#2 PBS+insert RC'd; ngRNA = pegRNA#2 spacer+scaffold, stored bound to
its tj-pegRNA; RTgDNA overlap with the between-nicks arm excluded.

- [ ] **C1** — `recoded_pridict.py`: `TwoPegDesign` dataclass (`pbs1`, `first_end`, `rcpbs2`, `last_end`,
  `ngrna` (=spacer+scaffold), `pbs1_len`, `rt_len`, `editing_position`, `target_strand`, `offtarget_risk`).
  `TwoPegPridictFake` returns deterministic values shaped like `parse_pegrna_csv`
  (`PBSrevcomp`/`RTrevcomp`/`PBSlength`/`RTlength`/`Target-Strand`/`Editing_Position`). TDD.
- [ ] **C2** — real adapter `TwoPegPridictReal` submits `pridict2_pegRNA_design.py batch` via
  `submit_pridict_batches`, parses with `pridict_parser`. **VERIFY ON CX3:** (a) the PRIDICT2 column that
  carries the PBS-past-the-nick homology ("RTgDNA") — currently only `RTrevcomp`/`RTlength` are parsed;
  (b) that PRIDICT exposes a per-spacer off-target-risk field (else fall back to Cas-OFFinder for spacer QC).
  Record findings in ADR.
- [ ] **C3** — orientation: pegRNA#1 `PBS1`+first-end **untouched**; pegRNA#2 `rcpbs2 = RC(PBS)`,
  `last_end = RC(insert)`; `ngrna = pegRNA#2.spacer + scaffold` (untouched). TDD: assembled extension =
  RC(TS) end-to-end for a known TS.
- [ ] **C4** — ngRNA lifecycle: bind ngRNA to its tj-pegRNA candidate; on discard, drop both. TDD.
- [ ] **C5** — RTgDNA exclusion: where PRIDICT's PBS-past-the-nick homology overlaps the between-nicks arm,
  exclude it (the arm already covers that span, to the nucleotide). TDD with a fake carrying an RTgDNA tail.

---

## 10. Stage D — Spacer QC via PRIDICT off-target risk

**Checkpoint D:** spacer QC uses PRIDICT's off-target risk (or Cas-OFFinder fallback if PRIDICT lacks it,
per C2); the TTTT/GC heuristic is removed from the selection path.
- [ ] **D1** — replace `_spacer_passes_sanity` usage with a PRIDICT-off-target-risk filter (threshold
  configurable). TDD both accept/reject. Keep the heuristic only as a last-resort fallback if the field is
  absent.

---

## 11. Stage E — recode-arm canonical for CDS, byte-exact for UTR

**Checkpoint E:** the between-nicks arm recodes CDS (microhomology-broken) and is byte-exact for UTR; the
fwd het-allele-aware pegRNA edits the UTR byte-for-byte; the arm span is asserted equal to the exact
between-nicks span.
- [ ] **E1** — confirm `_partition_arm` recodes CDS runs + locks non-coding runs; add an explicit test that
  the arm span == `genome[min(S,nick):max(S,nick)]` to the nucleotide (the user's hard requirement).
- [ ] **E2** — wire the UTR byte-exact path for the het-allele-aware fwd pegRNA (locked_span + het allele).
  TDD: a het SNP in the UTR arm is preserved byte-for-byte.

---

## 12. Stage F — Assemble + full gate stack + iterative search

**Checkpoint F:** the assembled tj-pegRNA runs **every** gate (verified active by an explicit test that
each gate is invoked); failures route to the correct lever per §5; the search is capped and emits precise
diagnostics; the winner is the highest-PE2 all-gate-passer.
- [ ] **F1** — assembly in TS convention (`Spacer·Scaffold·rcPBS2·RTT·PBS1·PolyU`, order per A3).
- [ ] **F2** — a test that asserts **each** gate function is actually called for a candidate (guard against
  a silently-skipped gate — the user's explicit requirement "MAKE SURE ALL OF THEM ARE ACTIVE").
- [ ] **F3** — implement the §5 gate→lever routing; TDD each mapping (a fake scorer forces each gate to
  fail and the test asserts the correct lever fires: NICK = new PRIDICT pair; INSERT = new CAI middle;
  HARD STOP for invariants + exhausted-AlphaFold).
- [ ] **F4** — locate splice/structure failures as **cassette vs recode-arm** to pick INSERT vs NICK
  (reuse the verification annotation model's span logic). TDD both.
- [ ] **F5** — AlphaFold gate upgrade (D14): fold/assess the cargo's impact at S (with tag/linker), INSERT
  lever first, hard-stop after insert variants exhausted. TDD.
- [ ] **F6** — winner selection: highest PRIDICT PE2 among all-gate-passers. TDD.

---

## 13. Stage G — Reverse + het-forward designs in parallel; unified output

**Checkpoint G:** for one insert+location the pipeline produces the reverse-targeting design AND the
het-specific forward design, each flagged, each with its bound ngRNA, each with a verification DOCX; runs
in parallel where possible.
- [ ] **G1** — orchestrate both designs (reverse-targeting; het-specific fwd, allele-aware UTR-locked).
- [ ] **G2** — flag outputs (`targeting = "reverse" | "forward_het"`) in the result, report, CSV, and
  verification DOCX. TDD both produced + flagged.
- [ ] **G3** — parallel execution (thread/process) with per-design isolation. TDD determinism.

---

## 14. Stage H — Merge cleanup + CX3 end-to-end

**Checkpoint H:** `recoded_tjpe_workflow.py` is deleted; all its tests migrated to the unified workflow;
the full suite + lint green; CX3 runs the real unified workflow end-to-end and produces a paired
tj-pegRNA+ngRNA with a verification DOCX; ADRs + docs updated; PR merged; origin+neuro+CX3 synced.
- [ ] **H1** — migrate/rename all `recoded_*` tests to the unified module; delete the separate module.
- [ ] **H2** — update `driver`/`recoded_driver` entry points to the single workflow.
- [ ] **H3** — CX3: stage the `pridict2` env; run `run_recoded_pipeline.pbs`; confirm RC(TS) end-to-end,
  exact between-nicks arm, ngRNA pairing, gate pass, DOCX render. Report the PRIDICT2 RTgDNA/off-target
  findings.
- [ ] **H4** — ADR(s), `FEATURES.md`, `DESIGN_HISTORY.md`; PR → merge → sync origin + neuro + CX3.

---

## 15. Open items to verify DURING build (non-blocking)
1. PRIDICT2 per-spacer off-target-risk field exists (else Cas-OFFinder fallback) — Stage C2/D.
2. Construct order (rcPBS2 5′-of-RTT vs 3′-terminal) — Stage A3, ADR 0014.
3. `ts_is_coding` per region (CDS vs UTR) — config flag vs MANE-transcript strand inference — Stage E.
4. PRIDICT2 "RTgDNA"/PBS-past-the-nick homology field name — Stage C2/C5.

## 16. Risk register
- **RC-flip errors** (the class this whole convention prevents): every extension component is RC(TS);
  `simulate_tprt` validates RTT + rcPBS2; a wrong flip fails loud. Mitigation: TDD every orientation.
- **PRIDICT compute** for 1000–5000 candidates × two ends × two designs: batch + cap; CX3 job time.
- **Silently-skipped gate:** Stage F2 asserts every gate is invoked.
- **Iteration non-termination:** capped loops per lever + precise exhaustion diagnostics.
- **UTR byte-exactness for het alleles:** locked_span + INV-0 byte-exact check (Decision 15).

## 17. References
- ADRs 0008–0013 (genome-wiring, CAI recode, off-target near-tier, cryptic-splice, verification DOCX,
  insertion-site importance) + 0014 (this convention).
- Zheng et al. 2023 (tjPE / rcPBS2). Yeo & Burge 2004 (MaxEntScan). Jaganathan 2019 (SpliceAI).
  Sharp & Li 1987 (CAI). PRIDICT2 (Mathis et al. 2023/2024).
- `docs/DESIGN_HISTORY.md`; prior plan `2026-07-02-pridict2-target-strand-recoded.md` (superseded by this).
