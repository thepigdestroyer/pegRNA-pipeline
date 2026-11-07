# Prime Assembly (A+) Implementation Plan

> **Status: executed, 2026-11-07.** Every task was done. What was built, the real-model evidence, and the places where the
> build differs from this plan and the spec are in [ADR 0030](../../decisions/0030-prime-assembly-design-pa.md) and in the
> spec's "As built" section. This file is kept as the record of how the work was ordered.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `pegrna design-pa`: two pegRNAs plus one donor whose arms restore the DNA between the nicks, every gate a real model, reachable from the CLI, the API and the site.

**Architecture:** Pure geometry and an in-silico simulator (`pa_geometry.py`, `pa_simulator.py`), then a ranked pair-walk (`recoded_pa_workflow.py`) that reuses the tj-PE arm partition, INV-0 and gates, with per-pegRNA gates memoized by the pegRNA's sequence. Driver, reports and CLI mirror `design-recoded`. Spec: `docs/superpowers/specs/2026-11-04-prime-assembly-design.md`. Structure: `docs/STRUCTURE.md`.

**Tech Stack:** Python 3.10+, Biopython, DnaChisel (optional `[recode]`), ViennaRNA (optional), pytest; FastAPI; React + TypeScript + Vite.

## Global Constraints

- Python 3.10 grammar (CI runs 3.10, 3.11, 3.12): an f-string must not reuse its own quote character inside a replacement field.
- `ruff check src tests` and `black --check src tests` are the CI lint. Line length is black's default for the repo.
- CI installs only the `[dev]` extras: no DnaChisel, no ViennaRNA. Tests that need either use `pytest.importorskip`.
- Everything is 5'→3'; antiparallel relationships use `reverse_complement` (`docs/STRUCTURE.md` §0). A cut `b` is the index of the first base right of the cut, for both strands.
- Flap length default 35, scanned 30 to 50; PBS lengths 13, 15, 17; span 30 to 100 with a hard floor of 23; minimum cargo `2 * F`; `max_run` 8; tail `UUU`; off-target tolerance 0 and exact duplicates never tolerated (spec §5).
- A PA pegRNA is at most 20 + 76 + 50 + 17 + 3 = 166 nt: always inside the calibrated docking range (96 to 232 nt).
- Defaults for the owner's four decisions (spec §8): PolyU tail; an unrecodable PAM warns and ranks lower; TTTT in an RTT or PBS rejects; a genomic run under a flap above `max_run` is an advisory.
- Never print or commit a secret. Deploy to Vercel only through `web/vercel/deploy.sh`.
- Commit messages avoid "Closes #N" wording; issues are closed by hand after CI.

## File structure

| File | Responsibility |
|---|---|
| `src/pegrna_pipeline/hpc/recode_arm.py` (modify) | deterministic constraint order; multi-start with fixed orderings (#30) |
| `src/pegrna_pipeline/hpc/nick_search.py` (modify) | one cut convention for both strands (#31) |
| `src/pegrna_pipeline/hpc/pa_geometry.py` (new) | guides, cuts, PAM-in pairs, arms, PAM breaking, donor, flaps, pegRNA assembly |
| `src/pegrna_pipeline/hpc/pa_simulator.py` (new) | `simulate_pa`: rebuild the outcome from the ordered pegRNAs and the donor |
| `src/pegrna_pipeline/hpc/recoded_pa_workflow.py` (new) | `run_recoded_pa`: ranked pair-walk, gates, `RecodedPAResult` |
| `src/pegrna_pipeline/hpc/pa_report.py` (new) | text report, components CSV, donor FASTA |
| `src/pegrna_pipeline/hpc/recoded_pa_driver.py` (new) | `design_pa`: config to scorers to walk to files |
| `src/pegrna_pipeline/agent/cli.py` (modify) | `design-pa`; `run --mode pa` points to it |
| `examples/ctcf_cterm/config_pa.json` (new) | CTCF + GGGGS-FLAG-P2A-ZeoR |
| `web/app.py`, `web/frontend/src/*`, `web/static/*` (modify) | `strategy` field, PA result, the site |
| removed | `pa_workflow.py`, `pair_scorer.py`, the `pa` branch of `driver.py` and their tests |
| tests | `test_recoder_determinism.py`, `test_nick_cut_convention.py`, `test_pa_geometry.py`, `test_pa_simulator.py`, `test_recoded_pa_workflow.py`, `test_pa_report_driver_cli.py`, `test_web_pa.py`; shared helper `tests/pa_helpers.py` |

---

## Stage 0: the two defects

### Task 0.1: a deterministic recoder (#30)

**Files:** Modify `src/pegrna_pipeline/hpc/recode_arm.py:98-156`. Test: `tests/test_recoder_determinism.py`.

**Interfaces:**
- Produces: `_avoid_kmers(genome: str, max_run: int) -> List[str]` (sorted, unique, length `max_run + 1`); `recode_min_edits_dnachisel(...)` unchanged in signature, now deterministic and multi-start.

- [ ] **Step 1: Failing tests.** (a) `_avoid_kmers` returns a sorted list and drops short tails (pure Python, runs in CI). (b) With `dnachisel`: the 91-nt CTCF arm `genome[56:147)`, frame 2, `max_run` 8, gives one identical arm across subprocesses with `PYTHONHASHSEED` 0, 1, 2, 5, 9, 11 (the seeds that disagreed in #30), and the edit count is at most the best seen then (14).
```python
def test_avoid_kmers_are_sorted_and_unique():
    assert _avoid_kmers("ACGTACGTAC", 3) == sorted(set(_avoid_kmers("ACGTACGTAC", 3)))
    assert all(len(k) == 4 for k in _avoid_kmers("ACGTACGTAC", 3))
```
- [ ] **Step 2: Run, expect FAIL** (`_avoid_kmers` not defined): `pytest tests/test_recoder_determinism.py -q`.
- [ ] **Step 3: Implement.** `_avoid_kmers`; solve under several fixed orderings (sorted; reversed; three `random.Random(seed).shuffle` orderings, seeds 1 to 3); keep the fewest edits, ties broken by the recoded string. Wrap in the existing `try`.
- [ ] **Step 4: Run, expect PASS**, then the whole suite for regressions (`pytest -q`).
- [ ] **Step 5: Real-model check, before and after.** On DSI run the CTCF demo before the change (record pegRNA and gate results) and after; both must pass every gate. The arm may differ: that is the point.
- [ ] **Step 6: Commit** `fix(recoder): deterministic constraint order, multi-start with fixed orderings`.

### Task 0.2: one cut convention for both strands (#31)

**Files:** Modify `src/pegrna_pipeline/hpc/nick_search.py:165-169`; modify `src/pegrna_pipeline/hpc/pridict_ngrna.py` only if the real-data check below shows the same skew. Test: `tests/test_nick_cut_convention.py`.

- [ ] **Step 1: Failing tests.** Build the genome around a known protospacer on each strand and compute the cut independently (the probe in issue #31): the forward guide's `nick_position` is `p0 + 17`; the reverse guide's is `q0 + 6`, where `q0` is the sense index of the PAM complement `CCN`.
- [ ] **Step 2: Run, expect FAIL** for the reverse case (25 returned, 26 expected).
- [ ] **Step 3: Fix.** The reverse scan maps `x` to `n - x` (not `n - 1 - x`). Update any existing test that encoded the old reverse-strand value.
- [ ] **Step 4: A forward-strand pegRNA design base by base.** Run `run_recoded_tjpe` with `pegrna_strand="Fw"` on a synthetic locus and assert the product equals `genome[:S] + cargo + genome[S:]` except inside the recoded arm, and that the arm spans exactly the bases between the two cuts.
- [ ] **Step 5: PRIDICT-path convention on real data (DSI).** Run PRIDICT2 on the CTCF window; for each nicking guide of either strand locate its protospacer in the window, compute the true cut, and compare with `site_s + distance_to_edit`. If they disagree by strand, normalise in `build_pridict_ngrnas` with a test; if they agree, record that in the test docstring.
- [ ] **Step 6: Commit** `fix(nick): one cut convention for both strands`.

---

## Stage 1: geometry and the simulator (pure, no models)

### Task 1.1: guides, cuts, PAM-in pairs

**Files:** Create `src/pegrna_pipeline/hpc/pa_geometry.py`, `tests/pa_helpers.py`. Test: `tests/test_pa_geometry.py`.

**Interfaces (produces):**
```python
MIN_NICK_SEPARATION_BP = 23; PAM_FOOTPRINT_CLASH_BP = 12
@dataclass(frozen=True)
class Guide:
    strand: str                         # "top" | "bottom"
    spacer: str                         # 20 nt, 5'->3', as in the guide RNA
    pam: str                            # 3 nt, 5'->3' on the PAM strand
    cut: int                            # genome[:cut] | genome[cut:]
    pam_span: Tuple[int, int]           # sense [lo, hi)
    protospacer_span: Tuple[int, int]   # sense [lo, hi)
def find_guides(genome: str) -> List[Guide]          # both strands; sorted by (cut, strand, spacer)
@dataclass(frozen=True)
class GuidePair:
    a: Guide; b: Guide
    @property
    def span(self) -> int: ...
def pam_in_pairs(guides: Sequence[Guide], *, site_s: int, min_span: int, max_span: int) -> List[GuidePair]
```

- [ ] **Step 1: Failing tests.** For each strand, construct a genome around a known protospacer and PAM and assert `cut` equals the hand-computed first-base-right-of-cut; `find_guides` finds it with the right spacer, pam and spans. `pam_in_pairs`: only top-left and bottom-right pairs; `a.cut <= S <= b.cut`; the span window; footprints not overlapping; `min_span < 23` raises `ValueError`.
```python
def test_cut_convention_both_strands():
    proto, pam = "GATCGATCGGCTAGCTTACG", "TGG"
    pad = "ACACACACACACACACACAC"
    top = pad + proto + pam + pad
    g = next(g for g in find_guides(top) if g.spacer == proto)
    assert g.cut == len(pad) + 17 and g.strand == "top"
    bottom = pad + rc(proto + pam) + pad
    g = next(g for g in find_guides(bottom) if g.spacer == proto)
    assert g.cut == len(pad) + 6 and g.strand == "bottom"
```
- [ ] **Step 2: Run, expect FAIL** (module missing).
- [ ] **Step 3: Implement** with the formulas of `scripts/analysis/probe_pa_geometry.py` (`cut = i - 3` for a top guide whose PAM starts at `i`; `cut = j + 6` for a bottom guide whose `CCN` starts at `j`).
- [ ] **Step 4: Run, expect PASS. Step 5: Commit** `feat(pa): guide enumeration, one cut convention, PAM-in pairs`.

### Task 1.2: arms, PAM breaking, the donor

**Files:** Modify `pa_geometry.py`. Test: `tests/test_pa_geometry.py`.

**Interfaces:**
```python
@dataclass(frozen=True)
class Arm:
    lo: int; hi: int; sense: str; original: str
    recode_span: Optional[Tuple[int, int]]; locked_span: Optional[Tuple[int, int]]
    locked_expected: Optional[str]; n_edits: int; max_run_after: int; microhomology_flag: bool
def build_arm(genome, lo, hi, cds_mask, max_run) -> Optional[Arm]      # empty Arm when lo == hi
@dataclass(frozen=True)
class PamOutcome:
    side: str                # "A" | "B"
    status: str              # "broken-by-recode" | "broken" | "intact-noncoding" | "intact-no-synonymous-change"
def break_pams(genome, arm5, arm3, guides: Sequence[Tuple[str, Guide]], cds_mask, max_run
              ) -> Tuple[Arm, Arm, List[PamOutcome]]
def assemble_donor(arm5: str, cargo: str, arm3: str) -> str
def product_allele(genome: str, n1: int, n2: int, donor: str) -> str  # genome[:n1] + donor + genome[n2:]
```
`build_arm` wraps `recoded_tjpe_workflow._partition_arm` (coding runs recoded, non-coding byte-exact) and caches by arguments.

- [ ] **Step 1: Failing tests.** `build_arm` of an all-coding span keeps the protein and the longest identical run is at most `max_run`; of a non-coding span is byte-exact with `microhomology_flag`; an empty span is an empty arm. `break_pams`: a PAM in a coding run whose GG can change synonymously is broken (status `broken`) and the edit is synonymous; a PAM in a locked run stays `intact-noncoding` and the arm is byte-identical; a Trp/Met codon (no synonym) gives `intact-no-synonymous-change`. `product_allele` equals `genome[:S] + cargo + genome[S:]` apart from the recoded bases; `len(product) == len(genome) + len(cargo)`.
- [ ] **Step 2: Run, expect FAIL. Step 3: Implement** (PAM positions: top guide bases `pam_span[0] + 1` and `+ 2` must stay `G`; bottom guide bases `pam_span[0]` and `+ 1` must stay `C`; codon from `_frame_at`; alternatives ordered by human codon frequency then lexicographically; a change must keep the arm's longest identical run at most `max_run`).
- [ ] **Step 4: Run, expect PASS. Step 5: Commit** `feat(pa): restored arms, PAM breaking, donor`.

### Task 1.3: flap scan, extensions, the pegRNA

**Files:** Modify `pa_geometry.py`. Test: `tests/test_pa_geometry.py`.

**Interfaces:**
```python
def flap_lengths(default: int = 35, lo: int = 30, hi: int = 50) -> List[int]   # 35, 36, 34, 37, 33, ... within [lo, hi]
@dataclass(frozen=True)
class Extension:
    side: str; flap_len: int; pbs_len: int; pbs: str; rtt: str
def derive_extension(side: str, genome: str, cut: int, donor: str, flap_len: int, pbs_len: int) -> Optional[Extension]
def choose_extension(side, genome, cut, donor, pbs_len, *, default=35, lo=30, hi=50) -> Optional[Extension]
def genomic_run_under_flap(side, genome, n1, n2, donor, flap_len) -> int
def assemble_pa_pegrna(*, spacer: str, rtt: str, pbs: str, tail: str = "UUU", scaffold: str = STANDARD_SGRNA_SCAFFOLD) -> str
def pa_regions(spacer: str, rtt: str, pbs: str, tail: str = "UUU") -> Dict[str, Tuple[int, int]]
```
Derivations (side `"A"` is the top-strand guide at `n1`; `"B"` the bottom-strand guide at `n2`):
`A: pbs = RC(genome[n1-L:n1]), rtt = RC(donor[:F])`; `B: pbs = genome[n2:n2+L], rtt = donor[-F:]`. `derive_extension` returns `None` when the PBS window leaves the genome or `F > len(donor)`. `choose_extension` takes the first `F` with no `TTTT` in RTT or PBS and a first RTT base other than `C`; failing that the first with no `TTTT`; else `None`.

- [ ] **Step 1: Failing tests** for `flap_lengths` order, both derivations on a hand-built locus, `TTTT` rejection, the C-start preference, the `None` cases, `genomic_run_under_flap`, `assemble_pa_pegrna` length (`20 + 76 + F + L + 3`) and `pa_regions` offsets. **Step 2: FAIL. Step 3: Implement. Step 4: PASS. Step 5: Commit** `feat(pa): flap scan, extensions, pegRNA assembly`.

### Task 1.4: the simulator

**Files:** Create `src/pegrna_pipeline/hpc/pa_simulator.py`. Test: `tests/test_pa_simulator.py`.

**Interfaces:**
```python
class PAMismatch(RuntimeError): ...
def simulate_pa(*, genome: str, a: Guide, b: Guide, donor: str, peg_a: str, peg_b: str,
                flap_a: int, pbs_a: int, flap_b: int, pbs_b: int, tail_len: int = 3,
                expected_product: Optional[str] = None) -> str     # returns the product
```
It cuts each ordered pegRNA into spacer, scaffold, RTT, PBS and tail by the given lengths, then checks: the spacer is the guide's protospacer on its strand and the PAM is NGG; `n1 < n2`; `RC(PBS_A) == genome[n1-L:n1]` and `PBS_B == genome[n2:n2+L]`; the flap each RTT would write equals the donor end it must match (`RC(RTT_A) == donor[:F_A]`, `RC(RC(RTT_B)) == donor[-F_B:]`); then returns `genome[:n1] + donor + genome[n2:]` and compares it with `expected_product` when given.

- [ ] **Step 1: Failing tests.** The happy path on a planted locus; a property test over many seeded random genomes, cargos and geometries; mutation tests: for each of A and B, replace RTT or PBS by its reverse, its complement, or (where a reverse complement belongs) its identity: every one raises `PAMismatch`; a wrong cut raises; a wrong spacer raises.
- [ ] **Step 2: FAIL. Step 3: Implement. Step 4: PASS. Step 5: Commit** `feat(pa): simulate_pa`.

---

## Stage 2: workflow, reports, command line

### Task 2.1: the walk

**Files:** Create `src/pegrna_pipeline/hpc/recoded_pa_workflow.py`. Test: `tests/test_recoded_pa_workflow.py` (synthetic locus from `tests/pa_helpers.py`, fake scorers; the real CTCF window behind `importorskip("dnachisel")`).

**Interfaces (produces):**
```python
class NoValidPADesign(NoValidPegRNA): ...
@dataclass
class PAGuideDesign:
    side: str; spacer: str; cut: int; pam: str; pam_status: str
    flap_len: int; pbs_len: int; pbs: str; rtt: str; tail: str
    pegRNA_sequence: str; offtargets: int; on_target: float
@dataclass
class RecodedPAResult:
    guide_a: PAGuideDesign; guide_b: PAGuideDesign
    site_s: int; span: Tuple[int, int]; cargo: str
    arm5: str; arm3: str; donor: str; product: str
    donor_parts: List[Tuple[str, str]]          # [("arm5", ..), ("cargo", ..), ("arm3", ..)]
    recode_edits: int; deletion_size: int       # always 0
    candidates_tried: int; pairs_considered: int
    warnings: List[str]; gate_results: List[Tuple[str, str]]
    recode_span: Optional[Tuple[int, int]]; locked_span: Optional[Tuple[int, int]]
def run_recoded_pa(*, cargo: str, genome: str, site_s: int, cds_mask: Sequence[bool],
    splice_scorer, full_structure_scorer, importance_scorer, dock_scorer, offtarget_scorer,
    maxent_scorer=None, pbs_lengths=(13, 15, 17), flap_length=35, flap_min=30, flap_max=50,
    min_span=30, max_span=100, max_run=8, tail="UUU", cas9="nCas9", structure_threshold=None,
    max_interface_drop=0.15, protein_length=None, insertion_importance_context=1,
    max_offtargets=0, min_ontarget=0.5, on_gate=None) -> RecodedPAResult
```
Walk: refuse `len(cargo) < 2 * flap_length`; enumerate guides; drop those failing `assess_spacer` or the batched off-target screen (`offtarget_rejection`); pair; rank by `(advisories, run_total, restored_bp, -on_target_sum, span, a.cut, b.cut)`; per pair build arms, break PAMs, donor, product, `assert_inv0`, choose extensions; then product-level gates (cryptic, splice, importance, the last memoized), then per-pegRNA gates (folding with the scaffold advisory, docking), memoized by pegRNA sequence; `simulate_pa` on the winner; fail loud with per-pair diagnostics.

- [ ] **Step 1: Failing tests:** refuses a short cargo; finds the planted pair; the result's donor is `arm5 + cargo + arm3` and the product equals `product_allele`; deletion size 0; INV-0 holds; each gate that rejects leaves a diagnostic (fake scorers that reject cryptic, splice, structure, docking); **each pegRNA is docked at most once** across pairs (count calls on a recording dock scorer); output identical across two runs; a fail-loud case with no pair; the standing caveat and PAM/flap-identity warnings present when they apply.
- [ ] **Step 2: FAIL. Step 3: Implement. Step 4: PASS. Step 5: Commit** `feat(pa): ranked pair-walk with memoized per-pegRNA gates`.

### Task 2.2: reports

**Files:** Create `src/pegrna_pipeline/hpc/pa_report.py`. Test: `tests/test_pa_report_driver_cli.py`.

**Interfaces:** `render_pa_report(result, manifest) -> str`, `render_pa_csv(result, manifest) -> str` (columns `index,component,role,part_of,sequence,length`), `render_pa_fasta(result) -> str`, `pa_gates_payload(result, scorers_stub, mode) -> dict` (the JSON the API serves, including a `pa` object with both guides, the donor parts and the span).

- [ ] Tests: the stand-in banner when gates are stand-ins; both pegRNAs, the donor with its parts, the order sheet and the standing caveat in the text; CSV rows in order with lengths; FASTA with three records; JSON keys. Implement, then commit `feat(pa): reports, CSV, FASTA, gates JSON`.

### Task 2.3: driver, command line, example

**Files:** Create `src/pegrna_pipeline/hpc/recoded_pa_driver.py`, `examples/ctcf_cterm/config_pa.json`; modify `src/pegrna_pipeline/agent/cli.py`. Test: `tests/test_pa_report_driver_cli.py`.

**Interfaces:** `design_pa(*, cargo, genome, site_s, cds_mask, mode="fake", output_dir=None, manifest=None, uniprot_acc="", genome_fasta=None, ..., progress_cb=None) -> RecodedPADesignOutput`; CLI `pegrna design-pa --config PATH --mode fake|real|api --output-dir DIR [--progress]`; files `recoded_pa_report.txt`, `recoded_pa_components.csv`, `recoded_pa_donor.fasta`, `recoded_pa_gates.json`. On failure the CLI writes the stand-in record as `design-recoded` does.

- [ ] Tests: CLI in fake mode writes the four files for the CTCF example (behind `importorskip("dnachisel")`) and for the synthetic config; a failing design exits non-zero with the diagnostics; `--progress` emits `@@P` lines. Implement; commit `feat(pa): design_pa driver and design-pa command`.

### Task 2.4: remove the legacy prototype

**Files:** delete `pa_workflow.py`, `pair_scorer.py`, the `pa` branch of `driver.py`, `tests/test_pa_*`, `tests/test_pair_scorer.py`, `tests/fixtures/pridict/pa_*`, the PA parts of `report_writer.py`/`annotator.py` if nothing else uses them; `pegrna run --mode pa` prints a pointer to `design-pa`; README and `docs` updated.

- [ ] Grep every reference first; remove only what nothing else uses; run the full suite; commit `chore(pa): remove the legacy prototype`.

### Task 2.5: real models on DSI, ADR 0030

- [ ] Rsync, run `design-pa --mode api` for CTCF + GGGGS-FLAG-P2A-ZeoR, assert: a design; no stand-in gates; both pegRNAs at most 166 nt; docking is a verdict (not advisory); INV-0; `simulate_pa`. Record the numbers in ADR 0030; mark the spec Implemented; commit.

---

## Stage 3: API, site, acceptance

### Task 3.1: API

**Files:** Modify `web/app.py`. Test: `tests/test_web_pa.py`.

- [ ] `DesignRequest.strategy: Literal["tjpe", "pa"] = "tjpe"`; both design endpoints pick the subcommand and the file prefix from it; `_validate` checks the PA rules (cargo at least `2 * flap`); `/api/example/pa` serves the example. Tests with the pipeline subprocess mocked, as `test_web_app.py` does. Commit `feat(web): strategy=pa on the design endpoints`.

### Task 3.2: the site

**Files:** Modify `web/frontend/src/*`, rebuild `web/static` (default `/static/` base only).

- [ ] A choice of how the insert is delivered; the PA pipeline steps; a result view with both pegRNAs, the donor map, copyable sequences, the FASTA download, the warnings callout; `tsc --noEmit` clean; the committed-bundle test passes. Commit `feat(ui): Prime Assembly mode`.

### Task 3.3: acceptance and deploy

- [ ] Redeploy DSI, run `scripts/audit-production.sh`, deploy the site with `web/vercel/deploy.sh --prod`, extend `scripts/qa-browser.mjs` to run a PA design end to end, run it on the live site, update `docs/FEATURES.md`, the README modes table and `web/deploy/GO-LIVE.md`, close #29 with the evidence. Commit.

---

## Self-review against the spec

- §2 acceptance 1 (real models) is Task 2.5; acceptance 2 (property and mutation tests) is Task 1.4; acceptance 3 is Task 3.3.
- §4.1 geometry, arms, PAM breaking, donor, flaps, the cargo rule: Tasks 1.1 to 1.3 and the refusal in 2.1. §4.2 walk, gates and ranking: 2.1. §4.3: 1.4. §4.4: 2.2, 2.3, 3.1, 3.2. §4.5: 2.4.
- §7 stage 0 is Stage 0 here.
- **Spec correction:** §4.2 described two independent half-designs. A guide's PAM can lie inside the other guide's arm (the CTCF best pair has B's PAM in arm5), so breaking it is a pair-level edit. The saving the spec wanted comes instead from memoizing per-pegRNA gates by the pegRNA's sequence. Spec §4.2 and §6 are updated to say so, and the "next pairs" item in §4.4 is dropped.
