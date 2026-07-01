# TJ-PE ngRNA Safety Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make TJ-PE select a paper-consistent, unique-hit nicking sgRNA partner for the best safe TJ-pegRNA, and make reports expose every required TJ-pegRNA/ngRNA component.

**Architecture:** Keep PRIDICT2 PE2 as the pegRNA/variant ranking engine, but add a safety-gated ngRNA-selection layer before final TJ-pegRNA assembly. The workflow will try variants in PE2 order, filter compatible ngRNAs by TJ-PE distance/strand/length/sanity, require exact-hit uniqueness when off-target data is supplied, and rank only surviving ngRNAs. Reporting will include ngRNA selection metadata and all TJ-pegRNA components.

**Tech Stack:** Python dataclasses and pure functions in `src/pegrna_pipeline/hpc/`; R/Biostrings off-target scanner in `scripts/score_guides.R`; pytest unit tests.

---

### Task 1: Extend TJ-PE result metadata and selection primitives

**Files:**
- Modify: `src/pegrna_pipeline/hpc/tjpe_workflow.py`
- Test: `tests/test_run_tjpe_workflow.py`
- Test: `tests/test_tjpe_edge_no_ngrna.py`

- [ ] **Step 1: Write failing tests for TJ-PE ngRNA filters**

Add tests that prove:
- default TJ-PE distance window is 70-90 bp;
- 20 nt ngRNAs outrank 19/21 nt when scores are otherwise comparable;
- candidates with `mm0_count != 1` are rejected when exact-hit metadata is available;
- the workflow falls back to the next PE2-ranked variant if the top variant has no safe ngRNA.

- [ ] **Step 2: Run targeted tests and confirm failure**

Run: `python -m pytest tests/test_run_tjpe_workflow.py tests/test_tjpe_edge_no_ngrna.py -q --tb=short`

Expected before implementation: failures showing missing exact-hit filtering / old 40-100 default.

- [ ] **Step 3: Implement metadata and selection**

Add optional fields to `NgRNADesign`: `exact_hit_count`, `deepspcas9`, `deephf`, `cfd_specificity` defaulting to `None` where unavailable. Add a pure `rank_tjpe_ngrnas()` helper that filters by strand, sanity, 19-21 nt, distance 70-90 by default, and exact-hit uniqueness when `require_unique_exact_hit=True`. Ranking key: safe exact hit first, 20 nt preference, available CFD/DeepHF/DeepSpCas9/PRIDICT DeepCas9 in descending order.

- [ ] **Step 4: Make workflow try variants in PE2 order**

Change `run_tjpe_workflow()` so it iterates ranked variants and returns the first PE2-ranked pegRNA with a safe ngRNA. If none survive, return/raise an explicit diagnostic result with no assembled oligo and reasons preserved.

- [ ] **Step 5: Verify targeted tests pass**

Run: `python -m pytest tests/test_run_tjpe_workflow.py tests/test_tjpe_edge_no_ngrna.py -q --tb=short`

Expected: all pass.

### Task 2: Wire exact-hit scoring before final TJ-PE selection

**Files:**
- Modify: `src/pegrna_pipeline/hpc/driver.py`
- Modify: `src/pegrna_pipeline/hpc/offtarget_scorer.py`
- Modify: `scripts/score_guides.R`
- Test: `tests/test_offtarget_scorer.py`
- Test: `tests/test_score_guides_script.py`

- [ ] **Step 1: Write failing tests for pre-selection exact-hit scan**

Add tests showing the driver submits a candidate ngRNA exact-hit scoring pass before final TJ-PE selection, and that `score_guides.R --max-mismatches 0` is supported.

- [ ] **Step 2: Run targeted tests and confirm failure**

Run: `python -m pytest tests/test_offtarget_scorer.py tests/test_score_guides_script.py -q --tb=short`

Expected before implementation: failure because there is no pre-selection scoring path.

- [ ] **Step 3: Add exact-hit scoring API**

Add a reusable `score_guides()` / `submit_offtarget_job(..., max_mismatches=0)` path that can score candidate ngRNAs before final result assembly. It must write `mm0_count` and leave higher mismatch buckets blank when `max_mismatches=0`.

- [ ] **Step 4: Update TJ-PE driver flow**

After PRIDICT batch parsing, rank variants by PE2, collect candidate ngRNAs for a bounded top-N variant set, exact-scan those candidates, attach `exact_hit_count` to candidates, then call `run_tjpe_workflow(require_unique_ngrna_exact_hit=True)`.

- [ ] **Step 5: Verify targeted tests pass**

Run: `python -m pytest tests/test_offtarget_scorer.py tests/test_score_guides_script.py -q --tb=short`

Expected: all pass.

### Task 3: Fix TJ-PE report clarity and caveats

**Files:**
- Modify: `src/pegrna_pipeline/hpc/annotator.py`
- Modify: `src/pegrna_pipeline/hpc/report_writer.py`
- Test: `tests/test_report_writer.py`

- [ ] **Step 1: Write failing report tests**

Assert that TJ-PE reports include: ngRNA selection window, spacer length, exact-hit count, PBS1, RC-PBS2, full ngRNA oligo, and no caveat implying a 19 nt multi-hit ngRNA is acceptable.

- [ ] **Step 2: Run report tests and confirm failure**

Run: `python -m pytest tests/test_report_writer.py -q --tb=short`

Expected before implementation: failures for missing exact-hit/window fields and stale caveat.

- [ ] **Step 3: Update annotations and CSV columns**

Add TJ-PE CSV/report fields for `ngrna_spacer_length`, `ngrna_exact_hit_count`, `ngrna_selection_window`, `ngrna_rank_reason`, and `pbs2_source_nick_position`.

- [ ] **Step 4: Verify report tests pass**

Run: `python -m pytest tests/test_report_writer.py -q --tb=short`

Expected: all pass.

### Task 4: Full verification and CTCF rerun decision

**Files:**
- Modify only if tests expose required follow-up.

- [ ] **Step 1: Run full test suite**

Run: `python -m pytest tests/ -q --tb=short`

Expected: exit code 0.

- [ ] **Step 2: Run lint checks**

Run: `python -m black --check . && python -m ruff check .`

Expected: exit code 0.

- [ ] **Step 3: Summarize whether existing CTCF result is invalidated**

The previous CTCF report is invalidated for synthesis because its ngRNA has `mm0_count = 4`. After code fix, rerun TJ-PE to find a safe ngRNA or explicitly report no safe ngRNA in the 70-90 bp TJ-PE window.

---

## Self-review

- Spec coverage: paper-backed TJ-PE architecture, ngRNA uniqueness, 19-21 nt length, 20 nt preference, default 70-90 bp window, report annotations, and tests are covered.
- Placeholder scan: no `TBD`/`TODO` placeholders remain.
- Type consistency: names match existing `TJPEResult`, `NgRNADesign`, `OffTargetReport`, and report writer patterns.
