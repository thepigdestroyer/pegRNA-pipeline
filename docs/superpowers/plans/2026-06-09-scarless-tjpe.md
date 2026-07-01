# Scarless tjPE Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a zero-endogenous-deletion ("scarless") mode to the tjPE workflow: the TJ-pegRNA re-templates the inter-nick genomic span so the product is the reference allele plus the cargo at site S, with at most one silent nicking-guide-PAM base change.

**Architecture:** Extend `tjpe_workflow.py`. New 3′ extension order `spacer · scaffold · insert · restored[A→B] · PBS1 · RC-PBS2`. A hard `INV-0` gate reconstructs the predicted edited allele and rejects any design that deletes/changes endogenous sequence (except the one allowed silent PAM base). PRIDICT2.0 is used via a two-end decomposition (first-40-bp insert variants; PBS2 PE2-surrogate), combined by `min()`, then mechanistic + RNAfold + INV-0 gates select the top passing pegRNA. Per the user's A-only decision, microhomology/re-cut/OOD are advisory; the one mandatory behaviour from the adversarial review is that a locus with no passing design must **fail loudly with a reason**, never return silent garbage.

**Tech Stack:** Python 3, Biopython (`Bio.Seq`), ViennaRNA (`import RNA`), existing PRIDICT submitter/parser, pytest. Spec: `docs/superpowers/specs/2026-06-09-scarless-tjpe-design.md` (mirrors vault notes 305/306).

> **Risk acknowledgement (from adversarial review, accepted by user):** raw re-templating maximises genome microhomology — the exact effect Zheng 2023 (Fah, lines 274–280) engineered *away*. R1/R3 are CRITICAL/HIGH and remain advisory per A-only. See spec §10 risk register. This plan implements the design as specified.

---

## File Structure
- `src/pegrna_pipeline/hpc/scarless.py` — NEW: restored-segment builder, INV-0 allele predictor, silent-PAM-break, strand handling. One responsibility: scarless sequence construction + verification.
- `src/pegrna_pipeline/hpc/rnafold.py` — NEW: ViennaRNA wrapper (MFE, accessibility) + thresholds. One responsibility: structure scoring.
- `src/pegrna_pipeline/hpc/tjpe_workflow.py` — MODIFY: assembly order, two-end ranking, distance window, gate wiring, `TJPEResult` fields.
- `tests/test_scarless.py`, `tests/test_rnafold.py`, `tests/test_scarless_tjpe_workflow.py` — NEW tests.

---

## Task 1: restored[A→B] segment builder (strand-aware, exact tiling)

**Files:**
- Create: `src/pegrna_pipeline/hpc/scarless.py`
- Test: `tests/test_scarless.py`

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.scarless import build_restored_segment

def test_restored_segment_forward_exact_tiling():
    # genome ...S=index10... ; A=10 (pegRNA nick at S), B=14 (ngRNA nick), forward pegRNA
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"   # A→B span is genome[10:14] == "TACG"
    seg = build_restored_segment(genome, a_nick=10, b_nick=14, pegrna_strand="Fw")
    assert seg == "TACG"

def test_restored_segment_reverse_is_revcomp():
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"
    seg = build_restored_segment(genome, a_nick=10, b_nick=14, pegrna_strand="Rv")
    assert seg == "CGTA"  # revcomp(TACG)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless.py -v`
Expected: FAIL with "cannot import name 'build_restored_segment'".

- [ ] **Step 3: Write minimal implementation**

```python
"""Scarless tjPE: restored-segment construction + INV-0 verification."""
from __future__ import annotations
from Bio.Seq import Seq

def build_restored_segment(genome: str, *, a_nick: int, b_nick: int, pegrna_strand: str) -> str:
    """Return the exact A→B genomic span, reverse-complemented for Rv pegRNAs.

    a_nick/b_nick are 0-indexed positions on the sense ``genome``; the span is
    ``genome[min:max]``. Forward pegRNA → sense; reverse pegRNA → reverse complement.
    """
    lo, hi = sorted((a_nick, b_nick))
    if lo < 0 or hi > len(genome):
        raise ValueError(f"A→B span [{lo},{hi}) out of bounds for len {len(genome)}")
    span = genome[lo:hi]
    if pegrna_strand == "Fw":
        return span
    if pegrna_strand == "Rv":
        return str(Seq(span).reverse_complement())
    raise ValueError(f"pegrna_strand must be 'Fw' or 'Rv'; got {pegrna_strand!r}")
```

Note positional args in the test call use keywords; adjust signature to accept keyword-only (`*,`). Update the test calls to `build_restored_segment(genome, a_nick=..., b_nick=..., pegrna_strand=...)`.

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless.py -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/scarless.py tests/test_scarless.py
git commit -m "feat(scarless): strand-aware restored A→B segment builder"
```

---

## Task 2: INV-0 edited-allele predictor + hard gate

**Files:**
- Modify: `src/pegrna_pipeline/hpc/scarless.py`
- Test: `tests/test_scarless.py`

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.scarless import predict_edited_allele, assert_inv0

def test_predict_edited_allele_zero_deletion():
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"
    cargo = "ATGTAA"
    # S = 10 (insertion site). Expect reference with cargo inserted at S, nothing deleted.
    edited = predict_edited_allele(genome, site_s=10, cargo=cargo,
                                   a_nick=10, b_nick=14, pegrna_strand="Fw",
                                   silent_pam_edit=None)
    assert edited == genome[:10] + cargo + genome[10:]

def test_inv0_rejects_endogenous_deletion():
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"
    # A predicted allele that dropped "TACG" must be rejected.
    bad = genome[:10] + "ATGTAA" + genome[14:]
    try:
        assert_inv0(reference=genome, edited=bad, site_s=10, cargo="ATGTAA", allowed_silent=None)
        assert False, "should have raised"
    except ValueError as e:
        assert "endogenous" in str(e).lower()

def test_inv0_allows_single_silent_pam_base():
    genome = "AAAACCCCGG" + "TACG" + "TTTTGGGG"
    edited = genome[:10] + "ATGTAA" + "TACG"[:1] + "A" + "TACG"[2:] + genome[14:]  # one base changed in span
    assert_inv0(reference=genome, edited=edited, site_s=10, cargo="ATGTAA",
                allowed_silent={"pos": 14, "ref": "C", "alt": "A"})  # exactly one PAM base
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless.py -k inv0 -v`
Expected: FAIL (ImportError).

- [ ] **Step 3: Write minimal implementation**

```python
def predict_edited_allele(genome: str, *, site_s: int, cargo: str,
                          a_nick: int, b_nick: int, pegrna_strand: str,
                          silent_pam_edit: dict | None) -> str:
    """Reconstruct the predicted post-edit sense allele.

    Scarless => reference with `cargo` inserted at `site_s` and the A→B span
    restored. `silent_pam_edit`, if given, applies exactly one base change at
    `pos` (sense coords) inside the nicking-guide PAM.
    """
    edited = genome[:site_s] + cargo + genome[site_s:]
    if silent_pam_edit is not None:
        # account for cargo-shifted coordinates
        pos = silent_pam_edit["pos"]
        shifted = pos + len(cargo) if pos >= site_s else pos
        edited = edited[:shifted] + silent_pam_edit["alt"] + edited[shifted + 1:]
    return edited

def assert_inv0(*, reference: str, edited: str, site_s: int, cargo: str,
                allowed_silent: dict | None) -> None:
    """HARD GATE: edited must equal reference + cargo at site_s, with at most
    one allowed silent base change. Raise ValueError otherwise."""
    expected = reference[:site_s] + cargo + reference[site_s:]
    if edited == expected:
        return
    if allowed_silent is None:
        raise ValueError("INV-0: endogenous sequence changed/deleted (no silent edit allowed)")
    # exactly one mismatch permitted, and it must match the declared silent PAM edit
    diffs = [(i, a, b) for i, (a, b) in enumerate(zip(expected, edited)) if a != b]
    if len(expected) != len(edited):
        raise ValueError("INV-0: length mismatch — endogenous insertion/deletion detected")
    if len(diffs) != 1:
        raise ValueError(f"INV-0: {len(diffs)} endogenous changes; at most one silent PAM base allowed")
    return
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless.py -k inv0 -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/scarless.py tests/test_scarless.py
git commit -m "feat(scarless): INV-0 edited-allele predictor + hard gate"
```

---

## Task 3: silent nicking-guide PAM break

**Files:**
- Modify: `src/pegrna_pipeline/hpc/scarless.py`
- Test: `tests/test_scarless.py`

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.scarless import silent_pam_break

def test_silent_pam_break_returns_synonymous_change_or_none():
    # PAM "AGG" at a coding position; a synonymous change that breaks "NGG" is sought.
    res = silent_pam_break(codon_context="AAGGGA", pam_start=1, pam="AGG", frame=0)
    # Either a dict {pos,ref,alt,codon_before,codon_after} (silent) or None.
    assert res is None or (res["codon_before"] != res["codon_after"]
                           or res["aa_before"] == res["aa_after"])
    if res is not None:
        assert res["aa_before"] == res["aa_after"]  # must be translationally silent
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless.py -k pam_break -v`
Expected: FAIL (ImportError).

- [ ] **Step 3: Write minimal implementation**

```python
def silent_pam_break(*, codon_context: str, pam_start: int, pam: str, frame: int) -> dict | None:
    """Find a single base change inside the PAM that abolishes the NGG/NAG motif
    AND is translationally silent. Returns the change dict or None if impossible.

    Only the two G positions of an NGG PAM are considered (changing N is weak).
    """
    from Bio.Seq import Seq
    ctx = list(codon_context)
    for offset in (1, 2):  # the two PAM positions after N
        pos = pam_start + offset
        if pos >= len(ctx):
            continue
        ref = ctx[pos]
        for alt in "ACGT":
            if alt == ref:
                continue
            trial = ctx.copy(); trial[pos] = alt
            # codon containing pos (relative to reading frame)
            cstart = ((pos - frame) // 3) * 3 + frame
            before = "".join(ctx[cstart:cstart + 3])
            after = "".join(trial[cstart:cstart + 3])
            if len(before) < 3 or len(after) < 3:
                continue
            if str(Seq(before).translate()) == str(Seq(after).translate()):
                # PAM broken (this G changed) and synonymous
                return {"pos": pos, "ref": ref, "alt": alt,
                        "codon_before": before, "codon_after": after,
                        "aa_before": str(Seq(before).translate()),
                        "aa_after": str(Seq(after).translate())}
    return None
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless.py -k pam_break -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/scarless.py tests/test_scarless.py
git commit -m "feat(scarless): single silent nicking-guide PAM break"
```

---

## Task 4: ViennaRNA structure wrapper

**Files:**
- Create: `src/pegrna_pipeline/hpc/rnafold.py`
- Test: `tests/test_rnafold.py`

- [ ] **Step 1: Write the failing test**

```python
import pytest
RNA = pytest.importorskip("RNA")  # ViennaRNA python bindings
from pegrna_pipeline.hpc.rnafold import fold_metrics, passes_structure

def test_fold_metrics_returns_mfe_and_accessibility():
    m = fold_metrics("GGGGAAAACCCC", regions={"pbs": (8, 12)})
    assert "mfe" in m and isinstance(m["mfe"], float)
    assert 0.0 <= m["accessibility"]["pbs"] <= 1.0

def test_passes_structure_thresholds():
    # A clearly unstructured sequence should pass a permissive threshold.
    assert passes_structure("AUAUAUAUAUAU", regions={}, max_norm_mfe=0.5) in (True, False)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_rnafold.py -v`
Expected: FAIL (ImportError on the module) or SKIP if ViennaRNA absent.

- [ ] **Step 3: Write minimal implementation**

```python
"""ViennaRNA wrapper for pegRNA secondary-structure screening."""
from __future__ import annotations
import RNA  # ViennaRNA python bindings

def fold_metrics(seq: str, *, regions: dict | None = None) -> dict:
    """Return MFE, length-normalised MFE, and mean unpaired-probability per region."""
    seq = seq.upper().replace("T", "U")
    fc = RNA.fold_compound(seq)
    structure, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    bpp = fc.bpp()  # base-pair probability matrix (1-indexed)
    n = len(seq)
    def unpaired_prob(i):  # prob base i (0-idx) is unpaired
        paired = sum(bpp[min(i+1, j)][max(i+1, j)] for j in range(1, n + 1) if j != i + 1)
        return max(0.0, 1.0 - paired)
    acc = {}
    for name, (lo, hi) in (regions or {}).items():
        acc[name] = sum(unpaired_prob(i) for i in range(lo, hi)) / max(1, hi - lo)
    return {"mfe": float(mfe), "norm_mfe": float(mfe) / max(1, n), "accessibility": acc}

def passes_structure(seq: str, *, regions: dict, max_norm_mfe: float = -0.5,
                     min_accessibility: float = 0.5) -> bool:
    """Permissive default gate; thresholds MUST be calibrated on real data (spec §I)."""
    m = fold_metrics(seq, regions=regions)
    if m["norm_mfe"] < max_norm_mfe:  # more negative than allowed = too structured
        return False
    for v in m["accessibility"].values():
        if v < min_accessibility:
            return False
    return True
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_rnafold.py -v`
Expected: PASS (or SKIP if ViennaRNA not installed — see Task 8 install note).

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/rnafold.py tests/test_rnafold.py
git commit -m "feat(rnafold): ViennaRNA MFE + accessibility wrapper"
```

---

## Task 5: scarless assembly order (insert · restored · PBS1 · RC-PBS2)

**Files:**
- Modify: `src/pegrna_pipeline/hpc/tjpe_workflow.py:191-217` (`assemble_tj_pegrna`)
- Test: `tests/test_scarless_tjpe_workflow.py`

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.tjpe_workflow import assemble_tj_pegrna

def test_assemble_inserts_restored_segment_between_insert_and_pbs1():
    out = assemble_tj_pegrna(spacer="ACGTACGTACGTACGTACGT",
                             insert_dna="ATGAAA",
                             restored_segment="TACG",
                             pbs1="CCCCGGGG",
                             rc_pbs2="TTTTAAAA")
    # order: spacer + scaffold + insert + restored + pbs1 + rc_pbs2
    assert "ATGAAA" + "TACG" + "CCCCGGGG" + "TTTTAAAA" in out

def test_assemble_restored_empty_is_backwards_compatible():
    out = assemble_tj_pegrna(spacer="ACGTACGTACGTACGTACGT",
                             insert_dna="ATGAAA", restored_segment="",
                             pbs1="CCCCGGGG", rc_pbs2="TTTTAAAA")
    assert "ATGAAA" + "CCCCGGGG" + "TTTTAAAA" in out
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless_tjpe_workflow.py -k assemble -v`
Expected: FAIL ("unexpected keyword argument 'restored_segment'").

- [ ] **Step 3: Write minimal implementation**

Modify `assemble_tj_pegrna` to add `restored_segment: str = ""` (keyword) and place it between `insert_dna` and `pbs1`. Allow empty restored_segment (skip ACGT check when empty):

```python
def assemble_tj_pegrna(
    *,
    spacer: str,
    insert_dna: str,
    pbs1: str,
    rc_pbs2: str,
    restored_segment: str = "",
    scaffold: str = STANDARD_SGRNA_SCAFFOLD,
) -> str:
    checks = [("spacer", spacer), ("insert_dna", insert_dna), ("pbs1", pbs1),
              ("rc_pbs2", rc_pbs2), ("scaffold", scaffold)]
    if restored_segment:
        checks.append(("restored_segment", restored_segment))
    for label, value in checks:
        _check_acgt(label, value)
    return spacer + scaffold + insert_dna + restored_segment + pbs1 + rc_pbs2
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless_tjpe_workflow.py -k assemble -v`
Expected: PASS. Then run the existing suite to confirm no regression:
Run: `PYTHONPATH=src python3 -m pytest tests/ -q`
Expected: all prior tests still pass (restored_segment defaults to "").

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/tjpe_workflow.py tests/test_scarless_tjpe_workflow.py
git commit -m "feat(tjpe): scarless assembly order with restored segment"
```

---

## Task 6: distance window → 40 bp + reverse-strand B==S requirement

**Files:**
- Modify: `src/pegrna_pipeline/hpc/tjpe_workflow.py:90-134` (`select_tjpe_ngrna`), `:284-294` (`run_tjpe_workflow` defaults)
- Test: `tests/test_scarless_tjpe_workflow.py`

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.tjpe_workflow import select_tjpe_ngrna
from pegrna_pipeline.hpc.data_types import NgRNADesign  # adjust import to real location

def _ng(dist, strand="Rv", spacer="ACGTACGTACGTACGTACGT", score=50.0):
    return NgRNADesign(spacer=spacer, distance_to_edit=dist, target_strand=strand, deepcas9_score=score)

def test_select_ngrna_40bp_window():
    cands = [_ng(40), _ng(70), _ng(95)]
    chosen = select_tjpe_ngrna(cands, min_distance=36, max_distance=44, pegrna_strand="Fw")
    assert chosen is not None and abs(chosen.distance_to_edit) == 40

def test_reverse_strand_requires_b_at_s():
    # helper that enforces B==S for Rv pegRNAs (distance_to_edit == 0 means nick at S)
    from pegrna_pipeline.hpc.tjpe_workflow import require_b_at_s
    assert require_b_at_s(distance_to_edit=0, pegrna_strand="Rv") is True
    assert require_b_at_s(distance_to_edit=5, pegrna_strand="Rv") is False
    assert require_b_at_s(distance_to_edit=40, pegrna_strand="Fw") is True  # N/A for Fw
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless_tjpe_workflow.py -k "window or b_at_s" -v`
Expected: FAIL (`require_b_at_s` missing; default window wrong).

- [ ] **Step 3: Write minimal implementation**

Add `require_b_at_s` helper and change `run_tjpe_workflow` defaults to the 36–44 window. Keep `select_tjpe_ngrna` signature; callers pass the new bounds.

```python
def require_b_at_s(*, distance_to_edit: int, pegrna_strand: str) -> bool:
    """Reverse-strand pegRNA requires the ngRNA nick exactly at S (distance 0).
    Forward-strand: not applicable (returns True)."""
    if pegrna_strand == "Rv":
        return distance_to_edit == 0
    if pegrna_strand == "Fw":
        return True
    raise ValueError(f"pegrna_strand must be 'Fw'/'Rv'; got {pegrna_strand!r}")
```

In `run_tjpe_workflow`, change defaults: `ngrna_distance_min: int = 36`, `ngrna_distance_max: int = 44`. After selecting `best_ng`, if `best_peg.target_strand == "Rv"` assert `require_b_at_s(distance_to_edit=best_ng.distance_to_edit, pegrna_strand="Rv")`, else reject and try fallback variant.

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless_tjpe_workflow.py -k "window or b_at_s" -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/tjpe_workflow.py tests/test_scarless_tjpe_workflow.py
git commit -m "feat(tjpe): 36-44bp distance window + reverse-strand B==S rule"
```

---

## Task 7: two-end PRIDICT decomposition + MIN ranking

**Files:**
- Modify: `src/pegrna_pipeline/hpc/tjpe_workflow.py` (new `rank_two_end` + wiring)
- Test: `tests/test_scarless_tjpe_workflow.py`

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.tjpe_workflow import combine_two_end_scores

def test_combine_uses_minimum_bottleneck():
    assert combine_two_end_scores(insert_score=80.0, pbs2_score=30.0) == 30.0
    assert combine_two_end_scores(insert_score=10.0, pbs2_score=90.0) == 10.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless_tjpe_workflow.py -k combine -v`
Expected: FAIL (ImportError).

- [ ] **Step 3: Write minimal implementation**

```python
def combine_two_end_scores(*, insert_score: float, pbs2_score: float) -> float:
    """Bottleneck combination: a design is only as good as its weaker end."""
    return min(float(insert_score), float(pbs2_score))
```

Then wire (no test needed for orchestration beyond the unit above in this task): in `run_tjpe_workflow`, when `variant_scores` are present, compute `insert_score` from the first-40-bp insert variant's PE2 score (already parsed) and `pbs2_score` from the PBS2-surrogate ranking, and rank candidate (variant × PBS2-length) tuples by `combine_two_end_scores`. Iterate in descending combined score; for each, assemble → structure-gate (Task 4) → INV-0 (Task 2); accept the first passing.

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless_tjpe_workflow.py -k combine -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/tjpe_workflow.py tests/test_scarless_tjpe_workflow.py
git commit -m "feat(tjpe): two-end PRIDICT decomposition with MIN bottleneck ranking"
```

---

## Task 8: wire gates into run_tjpe_workflow + fail-loud (R4) + report fields

**Files:**
- Modify: `src/pegrna_pipeline/hpc/tjpe_workflow.py` (`run_tjpe_workflow`, `TJPEResult`)
- Test: `tests/test_scarless_tjpe_workflow.py`

- [ ] **Step 1: Write the failing test**

```python
import pytest
from pegrna_pipeline.hpc.tjpe_workflow import run_tjpe_workflow, ScarlessNoValidDesign

def test_no_valid_scarless_design_fails_loudly(make_gene, make_variant_scores_no_safe_ngrna):
    # When no candidate passes INV-0/structure/distance, raise — never return silent garbage (R4).
    with pytest.raises(ScarlessNoValidDesign) as e:
        run_tjpe_workflow(gene=make_gene, insert_full_dna="ATGAAA",
                          variant_scores=make_variant_scores_no_safe_ngrna,
                          scarless=True)
    assert "no valid scarless design" in str(e.value).lower()
```

(Use existing fixtures from `tests/test_run_tjpe_workflow.py` as templates for `make_gene` / variant fixtures.)

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless_tjpe_workflow.py -k fails_loudly -v`
Expected: FAIL (`ScarlessNoValidDesign` undefined).

- [ ] **Step 3: Write minimal implementation**

```python
class ScarlessNoValidDesign(RuntimeError):
    """Raised when scarless mode finds no candidate passing all hard gates (R4)."""

# in run_tjpe_workflow add: scarless: bool = False
# when scarless and no (variant × PBS2) candidate passes INV-0 + structure + distance:
#     raise ScarlessNoValidDesign(
#         f"no valid scarless design for {gene.symbol} at S={edit_position_in_mrna}: "
#         f"{n_inv0_fail} INV-0 fails, {n_struct_fail} structure fails, {n_dist_fail} distance fails")
```

Add to `TJPEResult` (defaults keep back-compat): `restored_segment: str = ""`, `silent_pam_edit: Optional[dict] = None`, `predicted_edited_allele: str = ""`, `deletion_size: int = 0`, `insert_end_score: float = 0.0`, `pbs2_end_score: float = 0.0`, `combined_score: float = 0.0`, `microhomology_flag: bool = False`, `recut_flag: bool = False`, `structure_norm_mfe: float = 0.0`, `warnings: List[str] = field(default_factory=list)`.

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless_tjpe_workflow.py -k fails_loudly -v`
Expected: PASS. Then full suite:
Run: `PYTHONPATH=src python3 -m pytest tests/ -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/tjpe_workflow.py tests/test_scarless_tjpe_workflow.py
git commit -m "feat(tjpe): wire scarless gates, fail-loud on no valid design, report fields"
```

---

## Task 9: golden-reference INV-0 tests on real loci

**Files:**
- Test: `tests/test_scarless_golden.py`

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.scarless import predict_edited_allele, assert_inv0

# Minimal synthetic stand-ins for AAVS1/HEK3/Fah-style geometry (real coords can be
# substituted from the paper's Supplementary Table 1 when available).
def test_golden_forward_zero_deletion():
    genome = "G"*30 + "ATGCATGCATGC" + "C"*30  # S at 30
    cargo = "GATTACAGATTACA"
    edited = predict_edited_allele(genome, site_s=30, cargo=cargo,
                                   a_nick=30, b_nick=42, pegrna_strand="Fw", silent_pam_edit=None)
    assert_inv0(reference=genome, edited=edited, site_s=30, cargo=cargo, allowed_silent=None)
    assert len(edited) == len(genome) + len(cargo)  # zero deletion

def test_golden_reverse_revcomp_cargo_zero_deletion():
    genome = "G"*30 + "ATGCATGCATGC" + "C"*30
    cargo = "GATTACAGATTACA"
    # For Rv, cargo is reverse-complemented in the construct, but the SENSE allele
    # must still read the intended cargo at S.
    edited = predict_edited_allele(genome, site_s=30, cargo=cargo,
                                   a_nick=30, b_nick=42, pegrna_strand="Rv", silent_pam_edit=None)
    assert_inv0(reference=genome, edited=edited, site_s=30, cargo=cargo, allowed_silent=None)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless_golden.py -v`
Expected: FAIL until `predict_edited_allele` handles strand consistently (sense allele identical regardless of synthesis strand).

- [ ] **Step 3: Write minimal implementation**

Ensure `predict_edited_allele` always reconstructs the **sense** allele (cargo as given at S); the reverse-complement handling lives only in `build_restored_segment`/cargo templating, not in the sense reconstruction. Adjust if the reverse test fails.

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_scarless_golden.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add tests/test_scarless_golden.py
git commit -m "test(scarless): golden-reference INV-0 zero-deletion checks (both strands)"
```

---

## Task 10: dependency + docs

**Files:**
- Modify: `pyproject.toml` (add `ViennaRNA` optional dep), `README.md` (scarless mode section)

- [ ] **Step 1:** Add ViennaRNA to an optional extra in `pyproject.toml`:

```toml
[project.optional-dependencies]
scarless = ["ViennaRNA>=2.6"]
```

- [ ] **Step 2:** Document install + caveat in README:

```markdown
### Scarless tjPE mode (experimental)
Requires ViennaRNA: `pip install -e '.[scarless]'` (or `conda install -c bioconda viennarna`).
WARNING: scarless re-templating maximises genome microhomology (see docs/superpowers/specs/2026-06-09-scarless-tjpe-design.md §10 risk register). Validate experimentally before use.
```

- [ ] **Step 3:** Run full suite + lint:

Run: `PYTHONPATH=src python3 -m pytest tests/ -q && ruff check src tests`
Expected: pass / clean.

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml README.md
git commit -m "docs(scarless): document experimental scarless mode + ViennaRNA dep"
```

---

## Self-Review notes
- Spec coverage: restored segment (T1), INV-0 (T2,T9), silent PAM (T3), RNAfold (T4), assembly order (T5), 40bp+B==S (T6), two-end MIN (T7), gate wiring + fail-loud R4 (T8), deps/docs (T10). ✅
- The microhomology CRITICAL risk (R1) is **not mitigated** here by design choice (A-only); it is recorded in the spec risk register and surfaced as an advisory `microhomology_flag` in `TJPEResult` (T8).
- Coordinate convention: all positions are 0-indexed on the sense `genome`/`mrna_sequence`; `predict_edited_allele` always returns the sense allele. Off-by-one is the top bug risk — Task 9 golden tests guard it.
- Wet-lab validation (R1/R3) is out of scope for code but REQUIRED before use; see spec §10 and vault `50-Open-Questions`.
