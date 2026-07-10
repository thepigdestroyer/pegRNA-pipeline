# Wiring PRIDICT2 into the live design flow (auto-spacer, local CPU) — Design Spec

Date: 2026-07-10
Status: Proposed (awaiting user review)
Builds on: `docs/superpowers/plans/2026-07-02-unified-tjpe-target-strand.md` (the unified
target-strand two-peg design — deep mechanics live there; this spec covers the *wiring* delta and
the DSI/CPU execution model).

## 1. Problem

PRIDICT2 is **not wired into the live design flow**. Evidence: `recoded_driver.py` has zero PRIDICT
references; the only `src/` file referencing `TwoPegPridict`/`combine_two_peg` is `recoded_pridict.py`
itself (an orphaned Stage-C module — built and unit-tested, never connected). The live path
(`design-recoded` CLI → `recoded_driver.design_recoded` → `run_recoded_tjpe`) takes a **hand-supplied
`spacer`** and only auto-enumerates the ngRNA nick / PBS / RTT / recode-arm around it
(`recoded_tjpe_workflow.py`: *"the spacer is FIXED by site S. The spacer is never varied."*).

The tool should design its own guides: given the target + insert location, **PRIDICT2 designs the
pegRNA** (spacer + scaffold + PBS + RTT), and the pipeline consumes it — especially important for the
het-edit path, where the spacer must not be guessed.

## 2. Goal

Given target inputs (genome window, `site_s`, insert/cargo, `rtt_aa`, strand), the pipeline:
1. Calls **PRIDICT2 (local, CPU, on the DSI box)** to design the pegRNA at the insert site.
2. **Keeps ownership of *where* the ngRNA (pegRNA#2) targets** — the pipeline enumerates candidate
   nick positions and gates them; it then *tells PRIDICT2 where to design pegRNA#2*, whose
   spacer+scaffold becomes the ngRNA. PRIDICT designs the guide; the **pipeline decides the target
   and gates it**.
3. Retains the current **recode-arm** (CDS synonymous / UTR byte-exact) and the full **gate cascade**
   (off-target, cryptic-splice, splice, importance, structure, docking).
4. Drops the required manual `spacer` input (kept as an advanced override).

Non-goals: replacing the gate cascade; changing the TS convention (already specced in the unified
plan); running PRIDICT on GPU/HPC.

## 3. Execution model (DSI, CPU) — VERIFIED

PRIDICT2 is installed and confirmed running on `.141` CPU:
- Location `/mnt/data2/mj1824/pegrna/PRIDICT2`; env `/mnt/data2/mj1824/pegrna/envs/pridict2`
  (micromamba from `pridict2_repo.yml`, python 3.10, `torch==2.0.1+cpu`).
- Deep model **~10s/sequence** on CPU (batch amortizes further).
- Invocation: `python pridict2_pegRNA_design.py batch --input-fname <BASENAME> --output-dir <RELDIR>
  --nicking --cores N`.
- **Quirk 1**: `--input-fname` is resolved relative to `PRIDICT2/input/` → pass the *basename*, write
  the CSV into `input/`.
- **Quirk 2**: `--output-dir` is relative to the repo dir and must **pre-exist** (PRIDICT does not
  `mkdir` it).
- Input CSV: `sequence_name,editseq` (editseq bracket notation: insertion `(+INSERT)`, sub `(REF/ALT)`,
  del `(-SEQ)`, with adequate flanks each side). CONFIRMED on the box 2026-07-10: `(+SEQ)` works,
  `(/SEQ)` errors.
- Outputs per row: `<name>_pegRNA_Pridict_full.csv` (parsed by `pridict_parser.parse_pegrna_csv`:
  `Spacer-Sequence, PBSrevcomp, RTrevcomp, PBSlength, RTlength, Target-Strand, Editing_Position,
  PRIDICT2_0_editing_Score_deep_HEK, …`) and `<name>_nicking_guides.csv` (`Nicking-Protospacer,
  Nicking-Position-to-edit, Nicking-PAMdisrupt, Target_Strand, DeepCas9score, …`).

## 4. Architecture

### 4.1 Local PRIDICT runner (`pridict_local.py`, new)
A `PridictLocalRunner` that mirrors the `TwoPegPridict` contract but runs a **local subprocess** (not
the SSH/PBS `submit_pridict_batches`, which is CX3-only):
- Config via env: `PRIDICT2_DIR`, `PRIDICT2_PYTHON` (the pridict2 env python), `PRIDICT2_CORES`.
- `run_batch(rows: list[{sequence_name, editseq}]) -> dict[name, (pegRNA rows, nicking rows)]`:
  writes `input/<uuid>.csv`, `mkdir -p <reldir>`, runs the batch subprocess, parses both output CSVs
  via `pridict_parser`, cleans up. **Batch mode**: one PRIDICT call designs pegRNA#1 *and* all
  candidate pegRNA#2 targets at once (amortizes model load, one subprocess per design request).
- Timeout + fail-loud with a captured-stderr diagnostic (no silent stub in `api`/`real`).
- A `PridictRunnerFake` (deterministic, no subprocess) for CI — shaped like the real parser output.

### 4.2 Editseq construction (`input_formatter`, extend)
- `format_bracket_insertion` already exists (referenced by the CX3 adapter). Add/verify:
  - pegRNA#1 editseq = insertion `(+INSERT)` at `site_s` in the genome window (adequate flanks;
    PRIDICT wants ≥ ~100 bp each side — validate window length, fail loud if too short).
  - pegRNA#2 editseq = the *nick-defining* edit at each pipeline-chosen nick position (a silent/marker
    edit that makes PRIDICT design a guide nicking there). Exact encoding to verify during
    implementation against PRIDICT semantics.

### 4.3 Flow (the wiring) — pipeline owns "where", PRIDICT owns "how"
In `run_recoded_tjpe` (kept, refactored to source guides from PRIDICT):
1. **pegRNA#1**: PRIDICT designs at `site_s` → `spacer`, `PBSrevcomp`, `RTrevcomp`, scaffold. This is
   the auto-spacer (replaces the fixed input).
2. **ngRNA target selection (pipeline)**: `enumerate_nick_sites` (existing, `ngrna_distance_min/max`
   window, opposite strand) → candidate nick positions → `screen_nick_candidates` off-target gate
   (existing). The **pipeline chooses** which positions are admissible.
3. **pegRNA#2 design (PRIDICT, told where)**: for the admissible nick positions, batch-PRIDICT designs
   a pegRNA#2 targeting each; its `spacer+scaffold` becomes the ngRNA.
4. **Assemble** via `combine_two_peg` (existing: `rcpbs2 = orient_rcpbs2(peg2.PBSrevcomp)`,
   `last_end = RC(peg2.RTrevcomp)`, `ngrna = peg2.spacer + scaffold`; ADR 0014/0016).
5. **Recode arm** (existing `_partition_arm`: CDS synonymous / UTR byte-exact) + **full gate cascade**
   (existing, unchanged order): recode/INV-0/rcPBS2/structure → splice → cryptic → importance →
   docking. First candidate passing all gates wins; `site_s` never moves.
6. Gate→lever iteration (per unified plan): a nick-side failure → try the next pipeline-chosen nick
   (new pegRNA#2 target); an insert-side failure → CAI insert variant. The pipeline drives both.

Reuse (no rewrite): `recoded_pridict.combine_two_peg`, `TwoPegDesign`, `strand_orientation.orient_rcpbs2`,
`pridict_parser.parse_pegrna_csv`, `enumerate_nick_sites`, `screen_nick_candidates`, `_partition_arm`,
the gate scorers, verification DOCX.

### 4.4 CLI / config / web
- `design_recoded`: `spacer` becomes **optional**; add `auto_spacer: bool = True`. When auto, the
  runner supplies the spacer; when a spacer is given, keep the current fixed-spacer path (override).
- CLI/config docs updated; `spacer` no longer required in `config.json`.
- Web (`web/app.py` + `static/index.html`): move the spacer field into an "Advanced override"
  section (default hidden); add het-edit inputs (`het_snps`, "design het-specific pair" toggle) wired
  to `design_recoded_with_het_edit`; update the help text + pipeline flowchart labels so the page
  **accurately documents** what is auto-designed vs user-supplied.

### 4.5 Het-edit
Wire the CLI/web to call `design_recoded_with_het_edit` (antisense primary + allele-specific
sense-window pass) when `het_snps` + a het-edit window are supplied; the sense-window pegRNA is also
PRIDICT-designed. Document that allele specificity is capped by natural heterozygosity (never
manufactured), per the existing driver.

### 4.6 Deploy
- Add the `pridict2` env build to `web/deploy/setup_envs.sh` (clone PRIDICT2, create env from
  `pridict2_repo.yml`, `pip install torch==2.0.1 --index-url .../cpu`).
- Add `PRIDICT2_DIR` / `PRIDICT2_PYTHON` / `PRIDICT2_CORES` to `web/deploy/pegrna.env`; `web/app.py`
  injects them into the `api`-mode subprocess env.
- Redeploy to `.141`, restart the service, re-run an `api` design end-to-end (auto-spacer + live
  Boltz-2), confirm a real design with no manual spacer.

## 5. Performance
PRIDICT ~10–30 s/design (CPU, batched) + Boltz-2 docking minutes ⇒ an `api` design can be 5–20 min.
The web app's concurrency cap (429) + timeout already bound this; document expected latency in the UI.

## 6. Testing (TDD)
- `pridict_local.py`: unit tests with `PridictRunnerFake` + a captured real-output fixture for the
  parser mapping; subprocess-invocation test (monkeypatched) covering the two path quirks.
- editseq builder: unit tests (insertion/flank/short-window fail-loud).
- `run_recoded_tjpe` PRIDICT path: integration test with fakes (auto-spacer produces a design;
  pipeline-chosen nick + gate cascade still gates; manual-override path unchanged).
- Real PRIDICT smoke on DSI (marked, not CI).
- Full suite green; ruff/black clean.

## 7. Rollout (stages, each its own PR + CI + DSI redeploy)
1. `pridict_local.py` runner + `PridictRunnerFake` + editseq builder (self-contained, TDD).
2. Wire into `run_recoded_tjpe`/`design_recoded` (auto-spacer; pipeline-owned ngRNA target).
3. CLI/config: spacer optional; deploy env (`setup_envs.sh`, `pegrna.env`).
4. Web: form (spacer→advanced), het-edit inputs, accurate feature docs.
5. Het-edit end-to-end wiring.
6. DSI redeploy + real `api` validation.

## 8. Open items to verify during implementation
- Exact pegRNA#2 editseq encoding that makes PRIDICT nick at a chosen position.
- PRIDICT minimum flank length (window-size validation).
- Whether to pre-filter CAI insert variants before PRIDICT (unified plan §CAI) — deferred to a later
  stage unless needed.
