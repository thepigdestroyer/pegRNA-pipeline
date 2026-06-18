# Recoded TJ-pegRNA Design Pipeline — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Redesign the TJ-pegRNA design workflow so it enumerates PBS1×synonymous-RTT variants through PRIDICT2 (reverse-targeting), derives the nicking guide + PBS2 from a second PRIDICT2 run over off-target-screened nick sites, recodes the deletion-replacement arm to the minimum edits that destroy microhomology while preserving every codon, applies splice/AlphaFold/protein-silence hard gates, surfaces allele-specific pairings + het-SNP advisories from the diploid line genome (with an opt-in `allele_hunt` mode that may move the *nick* — never the insert site S — to capture a discriminating het, byte-exact in non-coding), and emits the first candidate from a ranked pool that survives a final full-pegRNA (RiNALMo 2D) structure check and a scaffold–nCas9 complex (Boltz-2) docking check — no full restart.

**Architecture:** A new orchestrator `recoded_tjpe_workflow.py` drives a **ranked-pool walk**: build candidate pool once (Stages 1–4), cheap gates prune early, expensive gates (splice, AlphaFold, RiNALMo, Boltz-2) run last on survivors, and failure advances to the next pre-ranked candidate rather than restarting. New pure modules (`microhomology.py`, `recode_arm.py`, `nick_search.py`, `genome_context.py`, `splice_screen.py`, `structure_importance.py`, `full_structure.py`, `complex_dock.py`, `allele_specificity.py`) are injected via Protocols so the core + CI stay light; heavy deps (DnaChisel, SpliceAI/Pangolin, AlphaFold/DSSP, RiNALMo, RNAformer, Boltz-2) live behind optional extras and lazy imports, mirroring the existing `rnafold.py` pattern.

**Tech Stack:** Python 3.11, Biopython, pandas, pytest; reuse `codon_expander.py` (synonymous enumeration), `offtarget_scorer.py` (Cas9 off-target HPC job), `pridict_*` (PRIDICT2 IO); new: DnaChisel, SpliceAI+Pangolin, AlphaFold-DB+DSSP+biotite+fpocket, **RiNALMo** (full-pegRNA 2D structure/accessibility — primary), **RNAformer** (CPU/CI-friendly cross-check), **Boltz-2** (pegRNA–nCas9 complex docking), ViennaRNA (kept as cheap inner-loop pre-screen). Genome: **T2T diploid KOLF2.1J assembly** (NHGRI `marbl/KOLF2.1J-Diploid-Assembly`) as sequence-of-truth + KOLF↔GRCh38 coordinate map for the GRCh38-trained gates.

---

## Scope & sub-project boundary

This plan covers the **design-logic redesign** (Stages 1–7 below) as one coherent pipeline. Two pieces are large enough to be their own sub-projects and are **stubbed with fail-loud guards + Protocol seams here, implemented in follow-on plans**:

- **Genome/KOLF2.1J ingest** (Stage 0) — load the **T2T diploid KOLF2.1J assembly** (both haplotypes) + a **KOLF↔GRCh38 coordinate map** (chain file or minimap2 alignment). This plan builds the `GenomicContext` type (now diploid-aware), the Protocol, and a fixture-backed fake; the real fetch/staging is a separate plan.
- **Heavy model adapters** (SpliceAI/Pangolin, AlphaFold importance, RiNALMo/RNAformer full-structure, Boltz-2 complex) — this plan builds the Protocol + in-repo fake + lazy adapter skeleton + optional extras; HPC weight-staging + Apptainer is a separate plan.

Everything else (enumeration reshape, nick search, ngRNA/PBS2 derivation, microhomology recode, allele-specific pairing search, het-SNP advisory, pool-walk control flow, reporting) is implemented fully and is testable with fakes in CI.

---

## Locked design decisions (from brainstorming sessions — B1/B2/B3 RESOLVED)

0. **AUTHORITATIVE tjPE pegRNA structure (B1 RESOLVED).** Full assembled pegRNA, 5′→3′:
   `Spacer · Scaffold · PBS2 · RTT(Insert + Deletion-Replacement-Arm) · PBS1 · PolyU(2–3 U)`
   - **PBS2 occupies the Homology-Arm slot** — in tjPE the HA is *replaced* by PBS2, immediately 3′ of the scaffold and 5′ of the RTT.
   - **PolyU (2–3 uridines)** is the 3′ free tail that binds the **La protein of the PE7 editor** (pegRNA 3′-end stabilisation / La recruitment).
   - PBS1 is now **internal** (just 5′ of PolyU), no longer the 3′ terminus.
   - Mechanism: PBS1 (3′-proximal) primes first-strand synthesis; RT copies the RTT 3′→5′; **PBS2 at the 5′ of the extension is copied last → becomes the 3′ jumping end** that hybridises to the second nick for the template jump + second-strand synthesis.
   - **IMPLEMENTATION-VERIFY (silent-catastrophic class):** at the post-scaffold PBS2 slot, confirm literal-PBS2 vs reverse-complement convention by unit assertion (the position-flip from the historical 3′ "RC-PBS2" likely flips RC→literal). Do not assume.
1. **Enumerate PBS1-length × synonymous-RTT variants only; spacer is fixed** by insertion site S. PRIDICT2 runs **reverse-targeting** (editing into CDS) and ranks; take best.
2. **Nick search 40–70 bp** into the genome (opposite strand to the pegRNA). All candidate nick sites scored by a **Cas9 spacer/on-target/off-target/PAM** tool first; only off-target-clean, PAM-valid sites survive.
3. **PBS2 + ngRNA derived from a 2nd PRIDICT2 run** (reverse-targeting pegRNA at the nick, fed ~last 40 bp of insert). KEEP scaffold+spacer (→ ngRNA) and PBS2 (→ paste into the PBS2 slot of the pegRNA, per decision 0); DELETE that run's Nick-RTgDNA + insert + homology arm.
4. **Deletion-replacement arm = recode of ALL gDNA between the two nicks**, including the Nick-RTgDNA stretch. Recode = **minimum synonymous edits that drop the longest contiguous perfect-identity run with gDNA below threshold L**, preserving every codon. If no synonymous solution exists → abandon candidate, advance.
5. **Hard gates:** protein-silence INV-0 (≤2 silent edits, decision 11); splice (SpliceAI+Pangolin, incl. symmetric-exon-skip); AlphaFold structural importance; conservation depth (≥50 orthologs else fail-loud). **Advisory:** ELM SLiMs, CAI, codon usage, microhomology score, sub-70 bp distance, het-SNP overlap (decision 12), allele-specific pairing availability (decision 13).
6. **Control flow = ranked-pool walk** (no full restart). Cheap structure pre-screen early (ViennaRNA), expensive **RiNALMo full-pegRNA 2D structure check + Boltz-2 complex docking check last**; first survivor wins; exhausting the pool widens search then fails loud.
7. **Final structure checks (TWO distinct gates, both last):**
   - **(7a) 2D self-pairing — RiNALMo (primary).** Folds the full pegRNA, scores self-pairing only among functional regions (spacer, PBS2, RTT, PBS1), **excluding scaffold** (designed to self-structure). RiNALMo emits an L×L symmetric base-pair-probability matrix → per-base accessibility + cross-region pairing; RoPE positional encoding → handles ≤~700 nt (no fixed length cap). **ViennaRNA** stays the cheap inner-loop pre-screen; **RNAformer** is the CPU/CI-friendly cross-check/fallback. (NuFold and RibonanzaNet were evaluated and REJECTED: RibonanzaNet caps at 177 nt with ±16 relpos clamp; NuFold is 3D-only.)
   - **(7b) Scaffold–nCas9 docking — Boltz-2 (final, heaviest).** Predicts the **full pegRNA + nCas9 complex**; compares apo scaffold–nCas9 vs candidate at the scaffold↔Cas9 interface (ipTM / interface-pLDDT + binding-affinity proxy); flags/rejects if docking degrades. MIT-licensed (commercial-safe). Cache the apo baseline; reuse across candidates. Fallbacks: Chai-1 (restraints) → Protenix → RoseTTAFold2NA.
8. **PBS2 gates:** length 11–17 nt, GC 40–70 %, Tm 37–42 °C (treated jointly — search may collapse to ~14–17 nt), not repetitive, ngRNA near a PAM (RC = NCC).
9. **Genome = T2T diploid KOLF2.1J assembly (B3 RESOLVED = B3-C hybrid).** The line's own **phased diploid assembly** (`marbl/KOLF2.1J-Diploid-Assembly`, NHGRI; UCSC hub `KOLF2.1Jv1.1`) is the **source of truth for all emitted/avoided sequence** (spacer, PBS1, PBS2, RTT/recode, off-target search, private-SNP avoidance) on **both haplotypes**. Each locus is **mapped to GRCh38** only to feed the GRCh38-trained gates (SpliceAI/Pangolin, AlphaFold-DB, annotations). Fail-loud on any locus that won't map. Never hg19. (Supersedes the original GRCh38+VCF-overlay plan.)
10. **Assembly order (B1 RESOLVED — see decision 0):** `Spacer · Scaffold · PBS2 · RTT · PBS1 · PolyU`.
11. **PAM destruction (B2 RESOLVED = B2-ii):** INV-0 permits **≤2 silent endogenous edits** — the **pegRNA PAM** break + the **ngRNA PAM** break — placed anywhere they are provably synonymous (NOT constrained to the recoded arm). Each break that lies **outside** the recoded arm is a new endogenous change and must **independently pass the splice gate AND the AlphaFold importance gate** at that residue; breaks **inside** the arm are absorbed by the recode. Hard-cap at 2 (a 3rd endogenous change = INV-0 fail). Report both edits (position, ref>alt, codon, which PAM).
12. **Het-SNP targeting = ADVISORY, never blocking, and often PREFERABLE.** A spacer/PBS/RTT overlapping a heterozygous SNP enables allele-specific editing → homogeneous all-het edited population (more consistent), vs homozygous-region targeting → mixed het+homozygous population. Flag (never gate) with: which element, **position class (PAM/seed/distal)**, predicted discrimination strength, consequence, and require the SNP be **confidently phased**. (Distinct from clone-validation: a Zeocin-surviving clone may still carry one clean + one damaged allele → allele-resolved long-read genotyping downstream.)
13. **Allele-specific pairing = always searched + surfaced (advisory feature).** When a het SNP sits in a discriminating position of the pegRNA (PBS1/spacer/PAM) **OR** the ngRNA (spacer/PAM), the pipeline **always designs and reports ≥1 single-allele-targeting pegRNA–ngRNA pairing**, even if far below the best-predicted pairing. **Both pegRNA-side and ngRNA-side discrimination are valid** (ngRNA-side: allele B sees only an SSB/nick → faithfully repaired, MMR reverts toward WT, and the large cargo cannot install without the 2nd nick → cleanly allele-specific; residual = low-frequency PE2-like byproduct on allele B, raised if PE4/PE5 MMR-inhibition is used → **report note, not a scar warning**). Strength hierarchy: PAM-destroying SNP > seed SNP > PBS1 SNP > distal (≈none). Constraints: SNPs exploited in both guides must be **same-haplotype (phased)**; the **other allele is scored as an explicit off-target**; specificity **cannot be manufactured** (the 2 silent PAM edits write to both alleles) — it is capped by natural heterozygosity (expect an allele-specific option for only ~5–10 % of targets; **compute the exact panel-wide count from the diploid genome and report it**).
14. **Insert size (near-term target): ~200–500 bp** → full pegRNA ≤ ~700 nt (drives the RiNALMo length choice; >800 bp would require a long-context model such as LoRNA).
15. **Mode-gated pegRNA-nick freedom for allele-specificity hunting (Route 4 — probability-raiser, NOT a guarantee).** The **insert site S NEVER moves** (cargo always lands at S; moving S is FORBIDDEN in every mode). **Only in `allele_hunt` mode**, the pegRNA *nick* may move 5′/3′ of its default to sweep its spacer/PAM/seed/PBS1 across more genomic territory and capture a discriminating het SNP. The extra gDNA between the moved nick and S is **re-templated (filled back) by the RTT**, partitioned by annotation:
    - **CDS bases → may be synonymously recoded** (decision 4 logic) to break microhomology.
    - **Non-coding bases (3′UTR / UTR / intron) → re-templated BYTE-EXACT to the *targeted allele's* sequence. ZERO changes, NO silent edits** (non-coding is base-sensitive: RBP / polyA / splice / miRNA).
    - **Accepted coupled tradeoff (advisory):** the 3′UTR is both where het SNPs are densest *and* where recoding is forbidden, so movement into it re-introduces 100 % microhomology on the locked stretch (the Zheng/Fah imperfect-edit mode). Mitigations: search nick movement **both directions** and report the tradeoff (5′-into-CDS = recodable but SNP-poor; 3′-into-3′UTR = SNP-rich but locked + microhomology); **minimise movement distance** (smallest stretch capturing a discriminating het → smallest locked span + shortest RTT); microhomology on the locked non-coding stretch is an **advisory flag** (wet-lab validate).
    - **B2-ii interaction:** "silent edit" is meaningless in non-coding, so the ≤2 silent PAM edits are **CDS-only**. A PAM in non-coding cannot be silently broken → rely on the het SNP itself disrupting that PAM (a PAM-destroying SNP is the strongest discriminator anyway) or leave it intact + re-cut advisory. A discriminating SNP inside the re-templated span is **locked to the targeted allele's base**.

---

## BLOCKERS — ALL RESOLVED (2026-06-18)

- **B1 — assembly linear order → RESOLVED.** Authoritative structure is `Spacer · Scaffold · PBS2 · RTT · PBS1 · PolyU` (decision 0). Supersedes the earlier insert-first/`RC-PBS2`-at-3′ assumption.
- **B2 — INV-0 edit count → RESOLVED = B2-ii.** ≤2 silent edits (both PAMs), each out-of-arm break independently splice+AlphaFold gated (decision 11).
- **B3 — genome source → RESOLVED = B3-C.** T2T diploid KOLF2.1J as sequence-of-truth, mapped to GRCh38 for ML gates (decision 9).

**Execution-time verifications (NOT design unknowns):** RC-vs-literal PBS2 convention at the new 5′ slot (assert in Task wiring); existence of a ready GRCh38↔KOLF2.1Jv1.1 chain file (else build via minimap2); RiNALMo/RNAformer/Boltz-2 exact max-length + I/O + weight-staging on CX3.

---

## File structure

**New modules (`src/pegrna_pipeline/hpc/`):**
- `microhomology.py` — pure: longest-common-substring vs genome, max-run metric, region-pair self-complementarity scoring. No heavy deps.
- `recode_arm.py` — minimum-edit synonymous recode (DnaChisel adapter + pure fallback), gated by INV-0. Lazy DnaChisel.
- `nick_search.py` — enumerate candidate nick sites in [40,70] bp, opposite strand; assemble off-target/PAM screen inputs; rank survivors.
- `genome_context.py` — **diploid-aware** `GenomicContext` type (both haplotypes) + exon map + KOLF↔GRCh38 coord conversions + `MissingGenomicContextError`. (Real fetch deferred; fixture-backed here.)
- `splice_screen.py` — `SpliceScorer` Protocol + lazy SpliceAI/Pangolin adapter + in-repo fake.
- `structure_importance.py` — `ImportanceScorer` Protocol (AlphaFold pLDDT/RSA/contacts/side-chain rules) + lazy adapter + fake.
- `full_structure.py` — `FullStructureScorer` Protocol + lazy **RiNALMo** adapter (primary) + **RNAformer** adapter (CPU/CI fallback) + fake; ViennaRNA pre-screen reused from `rnafold.py`. Scores cross-region pairing (scaffold excluded).
- `complex_dock.py` — `ComplexDockScorer` Protocol + lazy **Boltz-2** adapter + fake; apo-vs-candidate scaffold–nCas9 interface comparison.
- `allele_specificity.py` — pure: het-SNP overlap detection across functional elements, position-class (PAM/seed/distal), discrimination-strength, phasing check, other-allele off-target framing, allele-specific pairing search, **and `allele_hunt`-mode nick-movement search** (decision 15: shift the nick — never S — partition re-templated span into CDS `recode_span` + non-coding `locked_span`, rank by smallest shift, flag the CDS-vs-3′UTR tradeoff + microhomology advisory).
- `recoded_tjpe_workflow.py` — the ranked-pool-walk orchestrator tying it together.

**Modified:**
- `data_types.py` — add diploid `GenomicContext`, extend `VariantScore`/result types as needed (defaulted, back-compatible).
- `scarless.py` — generalise `assert_inv0` to a span-aware invariant taking **a list of ≤2 `allowed_silent` edits** + `recode_span` + injected `protein_check`; legacy `recode_span=None, allowed_silent=None` behaviour byte-identical.
- `pridict_submitter.py` / `pridict_parser.py` — support the 2nd (nick-side) PRIDICT run inputs/outputs (additive).
- `tjpe_workflow.py` — rewrite `assemble_tj_pegrna` to the authoritative order `Spacer · Scaffold · PBS2 · RTT · PBS1 · PolyU` (see Task 11b) + docstring structure note.
- `report_writer.py` / `annotator.py` — surface new fields (2 silent edits, het-SNP advisories, allele-specific pairings, panel-wide allele-specific count) + provenance manifest.

**Reused as-is:** `codon_expander.py` (`enumerate_variants`, `optimise_bulk`, `_gc_content`, `_has_homopolymer`, `_has_restriction_site`, `_has_internal_stop`), `offtarget_scorer.py` (`submit_offtarget_job`, `parse_offtarget_tsv`, `OffTargetReport`), `rnafold.py` (`fold_metrics`, `passes_structure`).

---

## Pipeline overview (target)

```
S0 ingest T2T-KOLF diploid (both haplotypes) + KOLF↔GRCh38 map ─► <genomic context present?> ─no─► FAIL-LOUD
   │
S1 enumerate PBS1-len × synonymous-RTT (spacer fixed) ─► PRIDICT2 reverse-target ─► rank ─► best pegRNA core
   │
S2 insert = PRIDICT2 first ~40bp + manual post-FLAG tail ─► <insert structure/synthesis/homology>
   │
S3 nick search [40,70]bp opposite strand ─► Cas9 off-target/PAM screen (+other-allele off-target) ─► clean nick sites
   │
S4 2nd PRIDICT2 run per nick (reverse pegRNA) ─► keep scaffold+spacer (=ngRNA) + PBS2 (→PBS2 slot); drop rest ─► rank PBS2
   │
S4b het-SNP scan (both haplotypes) ─► always emit ≥1 allele-specific pairing if available  (advisory)
   │
S5 recode deletion-replacement arm (min edits kill microhomology, preserve codons)
   │   └─ no synonymous solution? ─► abandon candidate, advance pool
   │
S6 HARD GATES: INV-0 (≤2 silent PAM edits) · splice · AlphaFold importance · conservation depth
   │   (each out-of-arm PAM break re-checked vs splice + AlphaFold)   (advisory: ELM/CAI/het-SNP)
   │
S7 assemble pegRNA [Spacer·Scaffold·PBS2·RTT·PBS1·PolyU] ─► ViennaRNA pre-screen
   │   └─ (7a) RiNALMo 2D self-pairing (functional regions, scaffold excluded)
   │   └─ (7b) Boltz-2 full pegRNA+nCas9 docking vs apo baseline
   │   └─ fails either? ─► abandon candidate, advance pool
   │
   ▼ first survivor ─► RANK + REPORT (+ allele-specific options + panel-wide count + provenance manifest)
```

---

## Task 1: Microhomology metric (pure)

**Files:**
- Create: `src/pegrna_pipeline/hpc/microhomology.py`
- Test: `tests/test_microhomology.py`

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.microhomology import longest_common_run, max_identity_run

def test_longest_common_run_finds_shared_substring():
    # 12-nt identical stretch embedded in otherwise different sequences
    a = "AAAACGTACGTACGTAGGGG"
    b = "TTTTCGTACGTACGTATTTT"
    assert longest_common_run(a, b) == len("CGTACGTACGTA")  # 12

def test_max_identity_run_against_self_offset():
    # identical sequences => full-length run
    s = "ACGTACGTACGT"
    assert max_identity_run(s, s) == len(s)

def test_no_shared_run_returns_zero_or_one():
    assert longest_common_run("AAAA", "GGGG") == 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_microhomology.py -v`
Expected: FAIL with "cannot import name 'longest_common_run'"

- [ ] **Step 3: Write minimal implementation**

```python
"""Pure microhomology metrics: longest contiguous identity between two strings.

Used by recode_arm.py to drive the deletion-replacement-arm recode toward the
minimum synonymous edits that break microhomology with the displaced genomic
strand, and by the full-pegRNA self-structure check.
"""
from __future__ import annotations

__all__ = ["longest_common_run", "max_identity_run"]


def longest_common_run(a: str, b: str) -> int:
    """Length of the longest contiguous substring common to ``a`` and ``b``.

    Classic DP (O(len(a)*len(b))). Sequences here are short (<=~120 nt) so this
    is fine; swap for a suffix-automaton only if profiling demands it.
    """
    a, b = a.upper(), b.upper()
    if not a or not b:
        return 0
    prev = [0] * (len(b) + 1)
    best = 0
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best:
                    best = cur[j]
        prev = cur
    return best


def max_identity_run(a: str, b: str) -> int:
    """Alias kept explicit for call-site readability (same as longest_common_run)."""
    return longest_common_run(a, b)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_microhomology.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/microhomology.py tests/test_microhomology.py
git commit -m "feat(microhomology): longest-common-run metric for arm recode"
```

---

## Task 2: Span-aware INV-0 (generalise the gate)

**Files:**
- Modify: `src/pegrna_pipeline/hpc/scarless.py:90-124` (`assert_inv0`)
- Test: `tests/test_scarless.py` (add cases), `tests/test_inv0_recode_span.py` (new)

**B2 RESOLVED (= B2-ii):** `allowed_silent` is a **list of ≤2** declared silent PAM edits (pegRNA PAM + ngRNA PAM). Max **2** changes outside the recode span; each must be single-base, translationally silent (verified by injected `protein_check`), and — if outside the recoded arm — independently pass the splice + AlphaFold gates (enforced by the Task 13 orchestrator, not here). A 3rd endogenous change ⇒ INV-0 fail. The Step-3 implementation below uses a `max_outside = len(allowed_silent or [])` capped at 2; `allowed_silent=None`/`recode_span=None` preserves legacy strict (≤1) behaviour byte-identically.

**Allele-hunt extension (decision 15):** when the nick is moved (`allele_hunt` mode), the re-templated arm partitions into a `recode_span` (CDS — synonymous changes allowed, protein-checked) and a `locked_span` (non-coding — **must equal the targeted allele's sequence byte-exact**; zero diffs, even silent). `assert_inv0` therefore also accepts an optional `locked_span` + `locked_expected` (the targeted-allele bytes) and raises if `edited[locked_span] != locked_expected`. Add a test: a single silent change inside `locked_span` must RAISE (`match="non-coding locked span"`), while the same change inside `recode_span` passes. When `locked_span=None` behaviour is unchanged.

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.scarless import assert_inv0
import pytest

def test_inv0_allows_many_synonymous_changes_inside_recode_span():
    ref   = "ATG" "AAA" "CCC" "GGG" "TAA"   # 15 nt
    cargo = "TTTTTT"
    # cargo inserted at S=12 (before stop); recode span [3,12) recoded synonymously
    edited = "ATG" "AAG" "CCA" "GGA" "TTTTTT" "TAA"  # 3 synonymous changes in [3,12)
    assert_inv0(
        reference=ref, edited=edited, site_s=12, cargo=cargo,
        recode_span=(3, 12), protein_check=lambda r, e: True,  # injected translator
        allowed_silent=None,
    )

def test_inv0_rejects_change_outside_recode_span():
    ref    = "ATGAAACCCGGGTAA"
    cargo  = "TTTTTT"
    edited = "ATC" "AAACCCGGG" "TTTTTT" "TAA"  # change at pos2 is OUTSIDE recode span
    with pytest.raises(ValueError, match="outside"):
        assert_inv0(
            reference=ref, edited=edited, site_s=12, cargo=cargo,
            recode_span=(3, 12), protein_check=lambda r, e: True, allowed_silent=None,
        )

def test_inv0_legacy_behaviour_unchanged_when_no_recode_span():
    ref    = "ATGAAACCCGGGTAA"
    cargo  = "TTT"
    edited = "ATGAAACCC" "TTT" "GGGTAA"  # clean insertion, no endogenous change
    assert_inv0(reference=ref, edited=edited, site_s=9, cargo=cargo,
                recode_span=None, allowed_silent=None)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_inv0_recode_span.py -v`
Expected: FAIL with "unexpected keyword argument 'recode_span'"

- [ ] **Step 3: Write minimal implementation**

Replace `assert_inv0` with the span-aware version:

```python
def assert_inv0(
    *,
    reference: str,
    edited: str,
    site_s: int,
    cargo: str,
    recode_span=None,        # Optional[Tuple[int,int]] half-open, on reference coords
    protein_check=None,      # Optional[Callable[[str,str],bool]]; ref_cds, edited_cds -> equal?
    allowed_silent=None,
) -> None:
    """HARD GATE — zero endogenous deletion + bounded changes.

    Outside ``recode_span`` (and the cargo), ``edited`` must equal ``reference``
    except at most the single declared ``allowed_silent`` PAM base. INSIDE
    ``recode_span`` any number of changes are allowed *provided* they are
    protein-silent — verified by the injected ``protein_check`` (so this module
    stays free of translation/coordinate concerns and is unit-testable).

    ``recode_span=None`` preserves the legacy strict behaviour byte-for-byte.
    """
    expected = reference[:site_s] + cargo + reference[site_s:]
    if len(expected) != len(edited):
        raise ValueError(
            "INV-0: length mismatch — endogenous insertion/deletion detected "
            f"(expected {len(expected)} nt, got {len(edited)} nt)"
        )
    # Map recode_span (reference coords) into edited coords (shift by cargo after S).
    def shift(pos: int) -> int:
        return pos + len(cargo) if pos >= site_s else pos

    span_edited = None
    if recode_span is not None:
        lo, hi = recode_span
        span_edited = (shift(lo), shift(hi))

    diffs = [i for i, (a, b) in enumerate(zip(expected, edited)) if a != b]
    outside = []
    for i in diffs:
        if span_edited is not None and span_edited[0] <= i < span_edited[1]:
            continue  # inside recode span — allowed if protein-silent (checked below)
        outside.append(i)

    # Outside-span: at most one allowed silent PAM base.
    if recode_span is None:
        if expected == edited:
            return
        if allowed_silent is None:
            raise ValueError(
                f"INV-0: {len(diffs)} endogenous change(s) but no silent edit allowed"
            )
        if len(diffs) != 1:
            raise ValueError(
                f"INV-0: {len(diffs)} endogenous changes; at most one silent PAM base allowed"
            )
        return

    max_outside = 1 if allowed_silent is not None else 0
    if len(outside) > max_outside:
        raise ValueError(
            f"INV-0: {len(outside)} change(s) outside recode_span/cargo "
            f"(max {max_outside} silent PAM base allowed)"
        )
    # Inside-span changes must be protein-silent.
    if protein_check is not None and not protein_check(reference, edited):
        raise ValueError("INV-0: recode_span changes are not protein-silent")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `PYTHONPATH=src python3 -m pytest tests/test_inv0_recode_span.py tests/test_scarless.py -v`
Expected: PASS (new + existing scarless tests green; legacy path unchanged)

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/scarless.py tests/test_inv0_recode_span.py
git commit -m "feat(inv0): span-aware invariant allowing synonymous recode-span changes"
```

---

## Task 3: Minimum-edit synonymous recode (pure fallback)

**Files:**
- Create: `src/pegrna_pipeline/hpc/recode_arm.py`
- Test: `tests/test_recode_arm.py`

This task implements the **pure** greedy recoder (no DnaChisel) so logic is testable in CI. Task 4 adds DnaChisel as the optimiser behind the same interface.

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.recode_arm import recode_min_edits, RecodeResult
from pegrna_pipeline.hpc.microhomology import longest_common_run
from Bio.Seq import Seq

def test_recode_breaks_microhomology_preserving_protein():
    # in-frame CDS span; genomic strand identical to it (worst-case microhomology)
    cds_span = "GAAGAAGAAGAAGAAGAA"  # Glu x6, all GAA
    genome   = cds_span
    res = recode_min_edits(cds_span, genome=genome, frame=0, max_run=6)
    # protein preserved
    assert str(Seq(res.recoded).translate()) == str(Seq(cds_span).translate())
    # microhomology broken below threshold
    assert longest_common_run(res.recoded, genome) <= 6
    # minimal-ish: at least one edit, fewer than every base
    assert 0 < res.n_edits < len(cds_span)

def test_recode_returns_none_when_no_synonymous_solution():
    # Met/Trp only => no synonymous alternatives => cannot break a long run
    cds_span = "ATGTGGATGTGG"  # Met-Trp-Met-Trp, zero degeneracy
    res = recode_min_edits(cds_span, genome=cds_span, frame=0, max_run=3)
    assert res is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_recode_arm.py -v`
Expected: FAIL with "cannot import name 'recode_min_edits'"

- [ ] **Step 3: Write minimal implementation**

```python
"""Minimum-edit synonymous recode of the deletion-replacement arm.

Goal: change the FEWEST bases needed to drop the longest contiguous
perfect-identity run with the displaced genomic strand below ``max_run``,
while keeping every codon synonymous (protein-identical). Pure greedy
implementation; recode_arm_dnachisel (Task 4) swaps in a constrained optimiser
behind the same signature.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

from Bio.Seq import Seq

from pegrna_pipeline.hpc.codon_expander import codons_for_aa
from pegrna_pipeline.hpc.microhomology import longest_common_run

__all__ = ["recode_min_edits", "RecodeResult"]


@dataclass(frozen=True)
class RecodeResult:
    recoded: str
    n_edits: int
    max_run_after: int


def _edit_distance(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def recode_min_edits(
    cds_span: str, *, genome: str, frame: int, max_run: int
) -> Optional[RecodeResult]:
    """Return a synonymous recode of ``cds_span`` whose longest common run with
    ``genome`` is <= ``max_run``, using a greedy minimum-edit heuristic, or
    ``None`` if no synonymous choice achieves it.

    ``frame`` is the reading frame offset of ``cds_span`` (0/1/2). Only whole
    codons fully inside the span are eligible to change.
    """
    if frame != 0:
        # Caller is expected to pass codon-aligned spans; document + guard.
        raise ValueError("recode_min_edits requires a codon-aligned span (frame=0)")
    codons = [cds_span[i : i + 3] for i in range(0, len(cds_span) - len(cds_span) % 3, 3)]
    best = list(codons)

    def assembled(cs: List[str]) -> str:
        return "".join(cs) + cds_span[len(cs) * 3 :]

    # Greedy: while a run is too long, find the codon overlapping the run's
    # midpoint and swap to the synonymous codon that most reduces the run with
    # the fewest base changes.
    for _ in range(len(codons)):  # bounded iterations
        cur = assembled(best)
        if longest_common_run(cur, genome) <= max_run:
            return RecodeResult(cur, _edit_distance(cur, cds_span), longest_common_run(cur, genome))
        improved = False
        for ci, codon in enumerate(best):
            aa = str(Seq(codon).translate())
            alts = [c for c in codons_for_aa(aa) if c != codon]
            if not alts:
                continue
            # choose alt minimising (resulting run, then base changes)
            scored = []
            for alt in alts:
                trial = best.copy()
                trial[ci] = alt
                run = longest_common_run(assembled(trial), genome)
                scored.append((run, _edit_distance(alt, codon), alt))
            scored.sort()
            run0, _, alt0 = scored[0]
            cur_run = longest_common_run(assembled(best), genome)
            if run0 < cur_run:
                best[ci] = alt0
                improved = True
                break
        if not improved:
            break

    final = assembled(best)
    run = longest_common_run(final, genome)
    if run <= max_run:
        return RecodeResult(final, _edit_distance(final, cds_span), run)
    return None
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_recode_arm.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/recode_arm.py tests/test_recode_arm.py
git commit -m "feat(recode): minimum-edit synonymous arm recoder (pure greedy)"
```

---

## Task 4: DnaChisel optimiser behind the same interface (optional dep)

**Files:**
- Modify: `src/pegrna_pipeline/hpc/recode_arm.py` (add `recode_min_edits_dnachisel` + dispatch)
- Modify: `pyproject.toml` (add `[recode]` extra: `dnachisel>=3.2`)
- Test: `tests/test_recode_arm_dnachisel.py`

- [ ] **Step 1: Write the failing test** (skips if DnaChisel absent)

```python
import pytest
from pegrna_pipeline.hpc.recode_arm import recode_min_edits_dnachisel
from pegrna_pipeline.hpc.microhomology import longest_common_run
from Bio.Seq import Seq

dnachisel = pytest.importorskip("dnachisel")

def test_dnachisel_recode_preserves_protein_and_breaks_run():
    cds_span = "GAAGAAGAAGAAGAAGAA"
    res = recode_min_edits_dnachisel(cds_span, genome=cds_span, frame=0, max_run=6)
    assert res is not None
    assert str(Seq(res.recoded).translate()) == str(Seq(cds_span).translate())
    assert longest_common_run(res.recoded, cds_span) <= 6
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_recode_arm_dnachisel.py -v`
Expected: FAIL with "cannot import name 'recode_min_edits_dnachisel'" (or SKIP if dnachisel missing — install with `pip install '.[recode]'`)

- [ ] **Step 3: Write minimal implementation**

```python
def recode_min_edits_dnachisel(cds_span, *, genome, frame, max_run):
    """Constrained-optimisation recode using DnaChisel.

    Constraints: EnforceTranslation (protein identity), AvoidMatches against the
    genomic strand longer than ``max_run``. Objective: minimise edits
    (MinimizeNumberOfChanges-equivalent). Returns RecodeResult or None.
    """
    try:
        from dnachisel import (
            DnaOptimizationProblem, EnforceTranslation, AvoidMatches,
            AvoidChanges,  # used with a small boost to bias toward fewer edits
        )
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("DnaChisel not available; install '.[recode]'") from exc
    if frame != 0:
        raise ValueError("recode_min_edits_dnachisel requires frame=0 (codon-aligned)")

    problem = DnaOptimizationProblem(
        sequence=cds_span,
        constraints=[
            EnforceTranslation(),
            AvoidMatches(match_length=max_run + 1, sequence=genome),
        ],
        objectives=[AvoidChanges(boost=1.0)],  # minimise number of edits
    )
    try:
        problem.resolve_constraints()
        problem.optimize()
    except Exception:
        return None
    recoded = problem.sequence
    from pegrna_pipeline.hpc.microhomology import longest_common_run
    run = longest_common_run(recoded, genome)
    if run > max_run:
        return None
    n = sum(1 for a, b in zip(recoded, cds_span) if a != b)
    return RecodeResult(recoded, n, run)
```

Add a `recode_arm(...)` dispatcher that prefers DnaChisel when importable, else the pure greedy `recode_min_edits`.

- [ ] **Step 4: Run test to verify it passes** (after `pip install '.[recode]'`)

Run: `PYTHONPATH=src python3 -m pytest tests/test_recode_arm_dnachisel.py -v`
Expected: PASS or SKIP (never error)

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/recode_arm.py tests/test_recode_arm_dnachisel.py pyproject.toml
git commit -m "feat(recode): DnaChisel constrained optimiser behind recode_arm dispatch"
```

---

## Task 5: Nick search over [40,70] bp opposite strand

**Files:**
- Create: `src/pegrna_pipeline/hpc/nick_search.py`
- Test: `tests/test_nick_search.py`

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.nick_search import enumerate_nick_sites, NickCandidate

def test_enumerate_nick_sites_in_window_opposite_strand():
    # 100-nt genome; pegRNA on Rv => nicks wanted on Fw; S at 30
    genome = "ACGT" * 25
    sites = enumerate_nick_sites(
        genome, site_s=30, pegrna_strand="Rv",
        distance_min=40, distance_max=70,
    )
    assert all(isinstance(s, NickCandidate) for s in sites)
    assert all(40 <= s.distance_to_edit <= 70 for s in sites)
    assert all(s.target_strand == "Fw" for s in sites)  # opposite of Rv pegRNA
    # every candidate must carry a PAM-anchored protospacer
    assert all(len(s.protospacer) == 20 for s in sites)

def test_no_pam_in_window_returns_empty():
    genome = "AAAAAAAA" * 13  # no NGG PAM on either strand
    assert enumerate_nick_sites(genome, site_s=20, pegrna_strand="Fw",
                                distance_min=40, distance_max=70) == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_nick_search.py -v`
Expected: FAIL with "cannot import name 'enumerate_nick_sites'"

- [ ] **Step 3: Write minimal implementation**

```python
"""Enumerate candidate nicking-guide sites in the [distance_min, distance_max]
window on the strand OPPOSITE the pegRNA, each anchored to an NGG PAM with a
20-nt protospacer. Off-target/on-target scoring is applied separately
(offtarget_scorer.submit_offtarget_job) to the protospacers returned here.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List

from Bio.Seq import Seq

__all__ = ["enumerate_nick_sites", "NickCandidate"]


@dataclass(frozen=True)
class NickCandidate:
    protospacer: str       # 20 nt
    pam: str               # 3 nt NGG
    nick_position: int     # 0-indexed on the sense genome
    distance_to_edit: int  # |nick - S|
    target_strand: str     # "Fw" or "Rv" (opposite the pegRNA)


def _opposite(strand: str) -> str:
    return {"Fw": "Rv", "Rv": "Fw"}[strand]


def enumerate_nick_sites(
    genome: str, *, site_s: int, pegrna_strand: str,
    distance_min: int, distance_max: int,
) -> List[NickCandidate]:
    genome = genome.upper()
    want = _opposite(pegrna_strand)
    out: List[NickCandidate] = []
    n = len(genome)

    def scan(seq: str, strand: str, to_sense):
        for i in range(len(seq) - 2):
            if seq[i] == "N":
                continue
            if seq[i + 1:i + 3] == "GG":  # NGG PAM
                # SpCas9 nick ~3 nt 5' of PAM on the protospacer
                proto_end = i  # protospacer is the 20 nt immediately 5' of PAM
                proto_start = proto_end - 20
                if proto_start < 0:
                    continue
                nick = to_sense(proto_end - 3)
                dist = abs(nick - site_s)
                if distance_min <= dist <= distance_max:
                    out.append(NickCandidate(
                        protospacer=seq[proto_start:proto_end],
                        pam=seq[i:i + 3],
                        nick_position=nick,
                        distance_to_edit=dist,
                        target_strand=strand,
                    ))

    if want == "Fw":
        scan(genome, "Fw", lambda x: x)
    else:
        rc = str(Seq(genome).reverse_complement())
        scan(rc, "Rv", lambda x: n - 1 - x)
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_nick_search.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/nick_search.py tests/test_nick_search.py
git commit -m "feat(nick-search): enumerate PAM-anchored nick sites in [40,70]bp window"
```

---

## Task 6: Off-target screen integration for nick candidates

**Files:**
- Modify: `src/pegrna_pipeline/hpc/nick_search.py` (add `screen_nick_candidates`)
- Test: `tests/test_nick_screen.py`

Reuses `offtarget_scorer.parse_offtarget_tsv` / `OffTargetReport`. The SSH/HPC submission is injected so the test uses a fake scorer.

**Note (decision 13):** for any candidate intended to be allele-specific, the **other haplotype of the same locus** is the most dangerous off-target. The screen must explicitly include the other-allele sequence (from the diploid `GenomicContext`, Task 7) as a scored off-target and confirm the het-SNP mismatch suppresses activity there. This is an additional `screen_nick_candidates` input (`other_allele_seqs`), defaulted empty for the non-allele-specific path.

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.nick_search import screen_nick_candidates, NickCandidate

def _fake_scorer(protospacers):
    # returns {protospacer: (ontarget, offtarget_count, has_offtarget)}
    return {p: (0.7, 0, False) if "ACGT" in p else (0.2, 9, True)
            for p in protospacers}

def test_screen_drops_offtarget_and_low_ontarget():
    cands = [
        NickCandidate("ACGT"*5, "AGG", 10, 45, "Fw"),   # clean
        NickCandidate("TTTT"*5, "AGG", 12, 47, "Fw"),   # offtarget-heavy
    ]
    kept = screen_nick_candidates(cands, scorer=_fake_scorer,
                                  min_ontarget=0.5, max_offtargets=0)
    assert len(kept) == 1
    assert kept[0].protospacer == "ACGT"*5
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_nick_screen.py -v`
Expected: FAIL with "cannot import name 'screen_nick_candidates'"

- [ ] **Step 3: Write minimal implementation**

```python
from typing import Callable, Dict, List, Tuple

def screen_nick_candidates(
    candidates: List["NickCandidate"], *,
    scorer: Callable[[List[str]], Dict[str, Tuple[float, int, bool]]],
    min_ontarget: float = 0.5,
    max_offtargets: int = 0,
) -> List["NickCandidate"]:
    """Keep only PAM-valid nick candidates whose protospacer is off-target-clean
    and has acceptable on-target activity. ``scorer`` maps protospacer ->
    (ontarget_score, offtarget_count, has_offtarget); inject the real Cas9
    off-target tool (offtarget_scorer.submit_offtarget_job) or a fake in tests.
    Returns survivors sorted by descending on-target score (deterministic).
    """
    scores = scorer([c.protospacer for c in candidates])
    kept = []
    for c in candidates:
        on, off, has_off = scores.get(c.protospacer, (0.0, 10**9, True))
        if has_off or off > max_offtargets or on < min_ontarget:
            continue
        kept.append((on, c))
    kept.sort(key=lambda t: (-t[0], t[1].distance_to_edit, t[1].protospacer))
    return [c for _, c in kept]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_nick_screen.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/nick_search.py tests/test_nick_screen.py
git commit -m "feat(nick-search): off-target/on-target screen for nick candidates"
```

---

## Task 7: GenomicContext type + coordinate conversions (fixture-backed)

**Files:**
- Create: `src/pegrna_pipeline/hpc/genome_context.py`
- Test: `tests/test_genome_context.py`

**B3 RESOLVED (= B3-C):** the genome source is the **T2T diploid KOLF2.1J assembly** (`marbl/KOLF2.1J-Diploid-Assembly`, NHGRI). This task builds the type, conversions, and fail-loud guard with a fixture; the real fetch + KOLF↔GRCh38 chain is a follow-on plan. **Diploid-aware:** `GenomicContext` carries **both haplotypes** (so guides/PBS are checkable on each allele — feeds het-SNP/allele-specific Task 17 and the zygosity story) and a `to_grch38(pos)` mapping used to feed the GRCh38-trained gates (SpliceAI/AlphaFold). Add `het_snps` (phased) to the fixture. The single-haplotype `mrna_to_genomic` test below still holds; add a dual-haplotype test.

- [ ] **Step 1: Write the failing test**

```python
import pytest
from pegrna_pipeline.hpc.genome_context import (
    GenomicContext, ExonBlock, MissingGenomicContextError, mrna_to_genomic,
)

def _ctx():
    # two exons: mRNA[0:10) <- genome[100:110); mRNA[10:20) <- genome[200:210)
    return GenomicContext(
        genome_build="GRCh38", chrom="chr1", strand=1,
        exons=(ExonBlock(100,110,0,10), ExonBlock(200,210,10,20)),
        genomic_sequence=None,
    )

def test_mrna_to_genomic_maps_across_exons():
    ctx = _ctx()
    assert mrna_to_genomic(ctx, 5) == 105
    assert mrna_to_genomic(ctx, 12) == 202

def test_mrna_to_genomic_out_of_range_raises():
    with pytest.raises(ValueError):
        mrna_to_genomic(_ctx(), 25)

def test_missing_context_guard():
    with pytest.raises(MissingGenomicContextError):
        mrna_to_genomic(None, 5)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_genome_context.py -v`
Expected: FAIL with import error

- [ ] **Step 3: Write minimal implementation**

```python
"""Genomic context for a transcript: exon map + coordinate conversions.

The real fetch (T2T diploid KOLF2.1J assembly + KOLF<->GRCh38 chain, B3-C) is a
separate plan; this module defines the (diploid-aware) type, the mRNA<->genomic
and KOLF<->GRCh38 conversions, and the fail-loud guard that forbids any
recoding/splice work when genomic context is absent.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple


class MissingGenomicContextError(RuntimeError):
    """Raised when a genome-dependent operation is attempted without context."""


@dataclass(frozen=True)
class ExonBlock:
    g_start: int
    g_end: int
    m_start: int
    m_end: int


@dataclass(frozen=True)
class GenomicContext:
    genome_build: str
    chrom: str
    strand: int  # +1 / -1
    exons: Tuple[ExonBlock, ...]
    genomic_sequence: Optional[str] = None


def mrna_to_genomic(ctx: Optional[GenomicContext], m_pos: int) -> int:
    if ctx is None:
        raise MissingGenomicContextError(
            "genomic context required (recoding/splice forbidden on mRNA-only data)"
        )
    for ex in ctx.exons:
        if ex.m_start <= m_pos < ex.m_end:
            offset = m_pos - ex.m_start
            return ex.g_start + offset  # + strand; strand=-1 handled by fetch layer
    raise ValueError(f"mRNA position {m_pos} not within any exon block")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_genome_context.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/genome_context.py tests/test_genome_context.py
git commit -m "feat(genome): GenomicContext type + mRNA->genomic mapping + fail-loud guard"
```

---

## Task 8: Splice scorer Protocol + fake + lazy adapter

**Files:**
- Create: `src/pegrna_pipeline/hpc/splice_screen.py`
- Modify: `pyproject.toml` (`[splice]` extra: `spliceai`, `pangolin` — heavy, HPC only)
- Test: `tests/test_splice_screen.py`

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.splice_screen import passes_splice_gate, SpliceFake

def test_splice_gate_blocks_high_delta():
    scorer = SpliceFake(delta=0.6)  # strong predicted disruption
    ok, detail = passes_splice_gate("ACGT"*30, "ACGA"*30, scorer=scorer, max_delta=0.2)
    assert ok is False
    assert "delta" in detail

def test_splice_gate_passes_low_delta():
    scorer = SpliceFake(delta=0.05)
    ok, _ = passes_splice_gate("ACGT"*30, "ACGT"*30, scorer=scorer, max_delta=0.2)
    assert ok is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_splice_screen.py -v`
Expected: FAIL with import error

- [ ] **Step 3: Write minimal implementation**

```python
"""Splice-disruption hard gate. SpliceScorer is a Protocol so SpliceAI/Pangolin
(heavy, HPC-staged weights) are injected; CI uses SpliceFake. Real adapters
require genomic context with >=5kb flanks (see genome_context) and score REF vs
ALT delta near annotated junctions, including symmetric-exon-skip risk.
"""
from __future__ import annotations

from typing import Protocol, Tuple


class SpliceScorer(Protocol):
    def delta(self, ref_seq: str, alt_seq: str) -> float: ...


class SpliceFake:
    def __init__(self, delta: float = 0.0) -> None:
        self._d = delta

    def delta(self, ref_seq: str, alt_seq: str) -> float:
        return self._d


def passes_splice_gate(ref_seq: str, alt_seq: str, *, scorer: SpliceScorer,
                       max_delta: float = 0.2) -> Tuple[bool, str]:
    d = scorer.delta(ref_seq, alt_seq)
    if d >= max_delta:
        return False, f"splice delta {d:.3f} >= {max_delta} (predicted disruption)"
    return True, f"splice delta {d:.3f} ok"
```

(Lazy real adapter `SpliceAIPangolinScorer` added in the heavy-deps follow-on plan; this commit ships Protocol + fake only.)

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_splice_screen.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/splice_screen.py tests/test_splice_screen.py pyproject.toml
git commit -m "feat(splice): SpliceScorer Protocol + fake + gate (real adapter deferred)"
```

---

## Task 9: AlphaFold structural-importance Protocol + fake + gate

**Files:**
- Create: `src/pegrna_pipeline/hpc/structure_importance.py`
- Modify: `pyproject.toml` (`[importance]` extra: `biotite`, `fpocket` note)
- Test: `tests/test_structure_importance.py`

Implements the hard-gate spec from the AlphaFold adversarial pass (pLDDT≥70 AND RSA<15%; contact-number; side-chain rules; in-frame indel in DSSP H/G/I/E/B; fail-loud on missing/low-identity model).

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.structure_importance import (
    residue_blocked, ResidueMetrics,
)

def test_buried_ordered_residue_blocks():
    m = ResidueMetrics(plddt=85, rsa=0.10, contact_number=20,
                       is_catalytic=False, is_disulfide=False,
                       is_positive_phi_gly=False, is_structural_pro=False)
    assert residue_blocked(m) is True  # ordered + buried core

def test_exposed_loop_residue_passes():
    m = ResidueMetrics(plddt=55, rsa=0.45, contact_number=4,
                       is_catalytic=False, is_disulfide=False,
                       is_positive_phi_gly=False, is_structural_pro=False)
    assert residue_blocked(m) is False

def test_catalytic_always_blocks():
    m = ResidueMetrics(plddt=40, rsa=0.6, contact_number=2,
                       is_catalytic=True, is_disulfide=False,
                       is_positive_phi_gly=False, is_structural_pro=False)
    assert residue_blocked(m) is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_structure_importance.py -v`
Expected: FAIL with import error

- [ ] **Step 3: Write minimal implementation**

```python
"""AlphaFold structural-importance hard gate (per-residue).

OR across hard-block axes (conservative: when in doubt, block). Metrics come
from an AlphaFold model (AFDB pLDDT/PAE + DSSP RSA + biotite contacts); the
ImportanceScorer Protocol lets the heavy extraction be injected/faked. Fail-loud
(handled by caller) when no model or footprint identity < 98%.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol


@dataclass(frozen=True)
class ResidueMetrics:
    plddt: float
    rsa: float
    contact_number: int
    is_catalytic: bool
    is_disulfide: bool
    is_positive_phi_gly: bool
    is_structural_pro: bool


def residue_blocked(m: ResidueMetrics, *, plddt_min: float = 70.0,
                    rsa_buried: float = 0.15, contact_max: int = 18) -> bool:
    if m.is_catalytic or m.is_disulfide or m.is_positive_phi_gly or m.is_structural_pro:
        return True
    if m.plddt >= plddt_min and m.rsa < rsa_buried:
        return True
    if m.contact_number > contact_max:
        return True
    return False


class ImportanceScorer(Protocol):
    def metrics(self, residue_index: int) -> ResidueMetrics: ...


def span_blocked(scorer: ImportanceScorer, residue_indices) -> bool:
    """Block if ANY residue in the footprint is structurally important."""
    return any(residue_blocked(scorer.metrics(i)) for i in residue_indices)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_structure_importance.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/structure_importance.py tests/test_structure_importance.py pyproject.toml
git commit -m "feat(importance): AlphaFold structural hard-gate (Protocol + fake + rules)"
```

---

## Task 10: Full-pegRNA 2D structure check (RiNALMo primary + ViennaRNA pre-screen + RNAformer fallback)

**Files:**
- Create: `src/pegrna_pipeline/hpc/full_structure.py`
- Modify: `pyproject.toml` (`[fullstructure]` extra: `rinalmo` + `rnaformer` notes — HPC/GPU for RiNALMo, conda-forge/CPU for RNAformer)
- Test: `tests/test_full_structure.py`

**Tool decision (locked):** primary = **RiNALMo** — RNA language model whose fine-tuned SS head emits an **L×L symmetric base-pair-probability matrix** (RoPE positional encoding → handles ≤~700 nt with no fixed length cap; best generalization to engineered/low-homology RNA, which is exactly the pegRNA OOD case). Cross-region pairing + per-base accessibility are read off that matrix. **ViennaRNA** is the cheap inner-loop pre-screen (`rnafold.py`); **RNAformer** (conda-forge, CPU-capable) is the CI/fallback scorer. NuFold (3D-only) and RibonanzaNet (177-nt cap, ±16 relpos clamp) were evaluated and REJECTED. The scorer is injected via a Protocol so this task's logic is tested with a fake (no GPU/weights in CI). Scores self-pairing only among functional regions (spacer, **PBS2**, RTT, PBS1), **excluding scaffold self-pairs** (scaffold is designed to self-structure). Note the region set uses the authoritative structure order (decision 0): PBS2 sits 5′ (after scaffold), PBS1 sits 3′ (before PolyU).

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.full_structure import (
    functional_region_pairing, region_accessibility, FullStructureFake,
)

def test_functional_pairing_excludes_scaffold():
    # fake returns an L x L prob matrix; scaffold pairs must be ignored by the scorer
    # pairs (i,j) with prob >= threshold: (2,90) spacer<->rtt ; (5,6) inside scaffold
    fake = FullStructureFake(pairs={(2, 90): 0.9, (5, 6): 0.9})
    regions = {"spacer": (0, 20), "scaffold": (20, 80),
               "pbs2": (80, 95), "rtt": (95, 120), "pbs1": (120, 133)}
    score = functional_region_pairing("A"*140, regions=regions, scorer=fake, threshold=0.5)
    # only (2,90) couples spacer<->pbs2 across functional regions
    assert score.cross_region_pairs == 1

def test_clean_pegrna_passes():
    fake = FullStructureFake(pairs={(25, 30): 0.9})  # both inside scaffold
    regions = {"spacer": (0, 20), "scaffold": (20, 80), "pbs2": (80, 95),
               "rtt": (95, 120), "pbs1": (120, 133)}
    score = functional_region_pairing("A"*140, regions=regions, scorer=fake, threshold=0.5)
    assert score.cross_region_pairs == 0

def test_region_accessibility_from_prob_matrix():
    # base 0 pairs with prob 0.8 => accessibility ~0.2
    fake = FullStructureFake(pairs={(0, 50): 0.8})
    acc = region_accessibility("A"*60, region=(0, 5), scorer=fake)
    assert 0.0 <= acc <= 1.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_full_structure.py -v`
Expected: FAIL with import error

- [ ] **Step 3: Write minimal implementation**

```python
"""Full-pegRNA 2D structure check.

Primary scorer = RiNALMo: its SS head emits an L x L symmetric base-pair
PROBABILITY matrix (RoPE => no fixed length cap, handles <=~700nt). We read off
(a) cross-region pairing that couples DISTINCT functional regions (spacer, PBS2,
RTT, PBS1) and (b) per-base accessibility = 1 - sum of pair probs. Scaffold rows/
cols are EXCLUDED (scaffold is designed to self-structure). FullStructureScorer
is a Protocol so RiNALMo (GPU) / RNAformer (CPU/CI) / a fake are interchangeable;
ViennaRNA (rnafold.py) stays the cheap inner-loop pre-screen. Threshold is
calibratable (RiNALMo val-tunes ~0.01-0.30). RiNALMo/RNAformer weights are
HPC-staged — validate max-length + I/O before relying.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Protocol, Tuple


class FullStructureScorer(Protocol):
    def pair_probs(self, seq: str) -> Dict[Tuple[int, int], float]:
        """Return {(i, j): probability} for i < j (sparse; absent => ~0)."""
        ...


class FullStructureFake:
    def __init__(self, pairs: Dict[Tuple[int, int], float]) -> None:
        self._pairs = pairs

    def pair_probs(self, seq: str) -> Dict[Tuple[int, int], float]:
        return self._pairs


@dataclass(frozen=True)
class StructureScore:
    cross_region_pairs: int


def _region_of(i: int, regions: Dict[str, Tuple[int, int]]):
    for name, (lo, hi) in regions.items():
        if lo <= i < hi:
            return name
    return None


def functional_region_pairing(seq: str, *, regions: Dict[str, Tuple[int, int]],
                              scorer: FullStructureScorer,
                              threshold: float = 0.5) -> StructureScore:
    pairs = scorer.pair_probs(seq)
    count = 0
    for (i, j), p in pairs.items():
        if p < threshold:
            continue
        ri, rj = _region_of(i, regions), _region_of(j, regions)
        if ri is None or rj is None:
            continue
        if ri == "scaffold" or rj == "scaffold":
            continue  # scaffold self/cross pairs ignored by design
        if ri != rj:
            count += 1
    return StructureScore(cross_region_pairs=count)


def region_accessibility(seq: str, *, region: Tuple[int, int],
                         scorer: FullStructureScorer) -> float:
    """Mean per-base unpaired probability over ``region`` = 1 - sum of pair probs."""
    pairs = scorer.pair_probs(seq)
    lo, hi = region
    width = max(1, hi - lo)
    paired = {k: 0.0 for k in range(lo, hi)}
    for (i, j), p in pairs.items():
        if lo <= i < hi:
            paired[i] += p
        if lo <= j < hi:
            paired[j] += p
    return sum(max(0.0, 1.0 - paired[k]) for k in range(lo, hi)) / width
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_full_structure.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/full_structure.py tests/test_full_structure.py pyproject.toml
git commit -m "feat(full-structure): RiNALMo-backed cross-region pairing + accessibility (scaffold-excluded)"
```

---

## Task 11: Reshape enumeration — PBS1-length × synonymous-RTT (spacer fixed)

**Files:**
- Modify: `src/pegrna_pipeline/hpc/tjpe_workflow.py` (or new `recoded_tjpe_workflow.py`) — add `enumerate_pbs_rtt_variants`
- Test: `tests/test_enumerate_pbs_rtt.py`

Reuses `codon_expander.enumerate_variants` for synonymous RTT variation.

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.recoded_tjpe_workflow import enumerate_pbs_rtt_variants

def test_enumeration_fixes_spacer_varies_pbs_and_rtt():
    out = enumerate_pbs_rtt_variants(
        spacer="ACGTACGTACGTACGTACGT",         # fixed
        rtt_aa="MKP",                           # 3 residues -> synonymous variants
        pbs_lengths=(11, 13, 15, 17),
        max_rtt_variants=8,
    )
    assert all(v.spacer == "ACGTACGTACGTACGTACGT" for v in out)
    assert {v.pbs_length for v in out} <= {11, 13, 15, 17}
    # synonymous RTT variants all translate to MKP
    from Bio.Seq import Seq
    assert all(str(Seq(v.rtt_dna).translate()) == "MKP" for v in out)
    # combinatorial but bounded
    assert len(out) <= 8 * 4
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_enumerate_pbs_rtt.py -v`
Expected: FAIL with import error

- [ ] **Step 3: Write minimal implementation**

```python
from dataclasses import dataclass
from typing import List, Tuple

from pegrna_pipeline.hpc.codon_expander import enumerate_variants


@dataclass(frozen=True)
class PbsRttVariant:
    spacer: str
    rtt_dna: str
    pbs_length: int


def enumerate_pbs_rtt_variants(
    *, spacer: str, rtt_aa: str, pbs_lengths: Tuple[int, ...],
    max_rtt_variants: int = 100,
) -> List[PbsRttVariant]:
    """Enumerate (synonymous RTT) x (PBS length); spacer is FIXED by site S.
    These are the candidates handed to PRIDICT2 (reverse-targeting) for ranking.
    """
    rtts = enumerate_variants(rtt_aa, max_variants=max_rtt_variants)
    out: List[PbsRttVariant] = []
    for rv in rtts:
        for L in pbs_lengths:
            out.append(PbsRttVariant(spacer=spacer, rtt_dna=rv.dna, pbs_length=L))
    return out
```

(Confirm the `CodonVariant` attribute name is `.dna`; adjust to the real field from `codon_expander.CodonVariant` if different.)

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_enumerate_pbs_rtt.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/recoded_tjpe_workflow.py tests/test_enumerate_pbs_rtt.py
git commit -m "feat(enumerate): PBS-length x synonymous-RTT variants, spacer fixed"
```

---

## Task 11b: Rewrite `assemble_tj_pegrna` to the authoritative structure (B1)

**Files:**
- Modify: `src/pegrna_pipeline/hpc/tjpe_workflow.py` (`assemble_tj_pegrna`, ~line 218-251)
- Test: `tests/test_assemble_authoritative_order.py` (new); update `tests/test_scarless_tjpe_workflow.py` assembly assertions

Authoritative 5′→3′ order (decision 0): `Spacer · Scaffold · PBS2 · RTT · PBS1 · PolyU`. This is a **breaking change** to the component order (PBS2 moves from 3′ to post-scaffold; PolyU tail added; PBS1 becomes internal). Update the signature + all callers. **Silent-catastrophic verify:** assert the PBS2-slot convention (literal vs reverse-complement) explicitly.

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.tjpe_workflow import assemble_tj_pegrna, STANDARD_SGRNA_SCAFFOLD

def test_authoritative_order_spacer_scaffold_pbs2_rtt_pbs1_polyu():
    out = assemble_tj_pegrna(
        spacer="ACGTACGTACGTACGTACGT",
        pbs2="GCAGTCGATCG",          # post-scaffold slot (HA replacement)
        rtt="ATGAAACCCGGGTTTAAATTT",  # insert + deletion-replacement-arm
        pbs1="CCCAGTAGCTAGC",
        polyu="UUU",
        scaffold=STANDARD_SGRNA_SCAFFOLD,
    )
    expected = "ACGTACGTACGTACGTACGT" + STANDARD_SGRNA_SCAFFOLD + \
        "GCAGTCGATCG" + "ATGAAACCCGGGTTTAAATTT" + "CCCAGTAGCTAGC" + "UUU"
    assert out == expected
    # PBS2 precedes RTT precedes PBS1; PolyU is the 3' terminus
    assert out.index("GCAGTCGATCG") < out.index("ATGAAACCCGGGTTTAAATTT") < out.index("CCCAGTAGCTAGC")
    assert out.endswith("UUU")

def test_polyu_only_uridines():
    import pytest
    with pytest.raises(ValueError):
        assemble_tj_pegrna(spacer="ACGT"*5, pbs2="GCAGTCGATCG",
                           rtt="ATGAAA", pbs1="CCCAGTAGCTAGC", polyu="UUA")  # A not allowed
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_assemble_authoritative_order.py -v`
Expected: FAIL — current signature is `(spacer, insert_dna, pbs1, rc_pbs2, restored_segment, scaffold)` with insert-first order

- [ ] **Step 3: Write minimal implementation**

```python
def assemble_tj_pegrna(
    *,
    spacer: str,
    pbs2: str,
    rtt: str,                 # = insert + deletion-replacement-arm (already recoded)
    pbs1: str,
    polyu: str = "UUU",       # 2-3 uridines; La/PE7 3' tail
    scaffold: str = STANDARD_SGRNA_SCAFFOLD,
) -> str:
    """Assemble the full tjPE pegRNA in the AUTHORITATIVE 5'->3' order.

    Layout (decision 0 / vault 306):

        spacer + scaffold + PBS2 + RTT + PBS1 + PolyU

    * PBS2 occupies the Homology-Arm slot (HA is REPLACED by PBS2 in tjPE),
      immediately 3' of the scaffold and 5' of the RTT. Its RT copy becomes the
      3' jumping end that hybridises to the second nick (template jump).
    * RTT = insert + deletion-replacement-arm (recoded; passed pre-assembled).
    * PBS1 primes first-strand synthesis; it is now INTERNAL (5' of PolyU).
    * PolyU = 2-3 uridines, the 3' free tail bound by the La protein of the PE7
      editor. Only 'U' permitted.

    IMPLEMENTATION NOTE (verify before wet-lab): the PBS2 string passed here must
    be in the orientation that makes its RT copy the correct jumping end. The
    historical 3'-terminal "RC-PBS2" convention likely flips to LITERAL PBS2 at
    this 5' slot — Task 13 asserts the predicted product against the genome to
    catch a wrong convention.
    """
    for label, value in (("spacer", spacer), ("scaffold", scaffold),
                         ("pbs2", pbs2), ("rtt", rtt), ("pbs1", pbs1)):
        _check_acgt(label, value)
    if not polyu or set(polyu.upper()) - {"U"}:
        raise ValueError(f"polyu must be 2-3 uridines (U only); got {polyu!r}")
    if not (2 <= len(polyu) <= 3):
        raise ValueError(f"polyu length must be 2-3; got {len(polyu)}")
    return spacer + scaffold + pbs2 + rtt + pbs1 + polyu
```

- [ ] **Step 4: Run test + update callers**

Run: `PYTHONPATH=src python3 -m pytest tests/test_assemble_authoritative_order.py tests/test_scarless_tjpe_workflow.py -v`
Expected: PASS after updating `run_tjpe_workflow`/`run_recoded_tjpe` call sites to the new kwargs (`pbs2=`, `rtt=`, `pbs1=`, `polyu=`). Region offsets used by Task 10 must follow this order: spacer, scaffold, **pbs2**, rtt, pbs1 (+ polyu).

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/tjpe_workflow.py tests/test_assemble_authoritative_order.py tests/test_scarless_tjpe_workflow.py
git commit -m "feat(assemble)!: authoritative tjPE order spacer-scaffold-PBS2-RTT-PBS1-PolyU (B1)"
```

---

## Task 12: PBS2 joint thermodynamic + positioning gate

**Files:**
- Modify: `src/pegrna_pipeline/hpc/recoded_tjpe_workflow.py` — add `passes_pbs2_gate`
- Test: `tests/test_pbs2_gate.py`

Reuses `codon_expander._gc_content`; Tm via nearest-neighbour (Biopython `Bio.SeqUtils.MeltingTemp`).

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.recoded_tjpe_workflow import passes_pbs2_gate

def test_pbs2_length_gc_tm_jointly():
    # in-range: 15 nt, GC ~53%, Tm in 37-42 band
    ok, _ = passes_pbs2_gate("ACGTGCACGTGCATG", repetitive=False)
    assert ok is True

def test_pbs2_rejects_too_short():
    ok, why = passes_pbs2_gate("ACGTGCAC", repetitive=False)  # 8 nt < 11
    assert ok is False and "length" in why

def test_pbs2_rejects_repetitive():
    ok, why = passes_pbs2_gate("ACGTGCACGTGCATG", repetitive=True)
    assert ok is False and "repetitive" in why
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_pbs2_gate.py -v`
Expected: FAIL with import error

- [ ] **Step 3: Write minimal implementation**

```python
from typing import Tuple
from Bio.SeqUtils import MeltingTemp as mt
from pegrna_pipeline.hpc.codon_expander import _gc_content


def passes_pbs2_gate(
    pbs2: str, *, repetitive: bool,
    min_len: int = 11, max_len: int = 17,
    gc_lo: float = 0.40, gc_hi: float = 0.70,
    tm_lo: float = 37.0, tm_hi: float = 42.0,
) -> Tuple[bool, str]:
    """Joint PBS2 gate: length 11-17, GC 40-70%, Tm 37-42C, not repetitive.
    Length and Tm are coupled — evaluated together (a short PBS2 rarely reaches
    37C unless GC-rich), so the search may collapse to ~14-17 nt; that is expected.
    """
    if not (min_len <= len(pbs2) <= max_len):
        return False, f"length {len(pbs2)} outside [{min_len},{max_len}]"
    gc = _gc_content(pbs2)
    if not (gc_lo <= gc <= gc_hi):
        return False, f"GC {gc:.2f} outside [{gc_lo},{gc_hi}]"
    tm = float(mt.Tm_NN(pbs2))
    if not (tm_lo <= tm <= tm_hi):
        return False, f"Tm {tm:.1f} outside [{tm_lo},{tm_hi}]"
    if repetitive:
        return False, "repetitive region"
    return True, "ok"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_pbs2_gate.py -v`
Expected: PASS (adjust the literal test sequence if its computed Tm lands outside 37–42; pick a real in-band 15-mer)

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/recoded_tjpe_workflow.py tests/test_pbs2_gate.py
git commit -m "feat(pbs2): joint length/GC/Tm + repetitive gate"
```

---

## Task 13: Ranked-pool-walk orchestrator

**Files:**
- Modify: `src/pegrna_pipeline/hpc/recoded_tjpe_workflow.py` — add `run_recoded_tjpe`
- Test: `tests/test_recoded_tjpe_workflow.py`

Ties Tasks 1–12 together with the **pool-walk** control flow: build ranked pool, cheap gates prune, expensive gates last, first survivor wins, exhaust → fail-loud. All heavy scorers injected (fakes in test).

**`allele_hunt` mode (decision 15).** `run_recoded_tjpe(..., allele_hunt: bool = False)`. When `False` (default): the pegRNA nick stays at its default position, **S never moves**, and the re-templated span is empty/CDS-only as today. When `True`: the orchestrator calls `enumerate_nick_movements` (Task 17) to try nick shifts that capture a discriminating het SNP — **S still never moves** — and for each shift partitions the re-templated span into `recode_span` (CDS, recoded) + `locked_span` (non-coding, byte-exact to the targeted allele, passed to `assert_inv0` per Task 2). Moving the insert site S is FORBIDDEN in both modes (assert S unchanged). Add a test that `allele_hunt=True` with a discriminating het returns an allele-specific design whose `locked_span` matches the targeted allele exactly and whose microhomology advisory is set; and that `allele_hunt=False` never moves the nick.

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.recoded_tjpe_workflow import run_recoded_tjpe, RecodedTJPEResult
from pegrna_pipeline.hpc.splice_screen import SpliceFake
from pegrna_pipeline.hpc.full_structure import FullStructureFake
from pegrna_pipeline.hpc.complex_dock import DockFake

def test_pool_walk_returns_first_survivor(monkeypatch_free_fakes):
    # fakes: splice clean, 2D structure clean, importance all-pass, docking preserved
    result = run_recoded_tjpe(
        **minimal_inputs(),
        splice_scorer=SpliceFake(delta=0.0),
        full_structure_scorer=FullStructureFake(pairs={}),     # no cross-region pairs
        importance_scorer=AllPassImportance(),
        dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
        offtarget_scorer=clean_offtarget,
    )
    assert isinstance(result, RecodedTJPEResult)
    assert result.pegRNA_sequence != ""
    assert result.deletion_size == 0
    assert result.candidates_tried >= 1

def test_pool_walk_advances_past_failing_candidate():
    # first candidate fails RiNALMo 2D (cross-region pairs), second passes
    result = run_recoded_tjpe(
        **minimal_inputs(n_variants=2),
        full_structure_scorer=FailFirstThenPass(),
        splice_scorer=SpliceFake(0.0),
        importance_scorer=AllPassImportance(),
        dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
        offtarget_scorer=clean_offtarget,
    )
    assert result.candidates_tried == 2

def test_pool_walk_advances_past_docking_failure():
    # 2D clean but Boltz-2 says the scaffold-Cas9 interface collapses on candidate 1
    result = run_recoded_tjpe(
        **minimal_inputs(n_variants=2),
        full_structure_scorer=FullStructureFake(pairs={}),
        splice_scorer=SpliceFake(0.0),
        importance_scorer=AllPassImportance(),
        dock_scorer=DockFailFirstThenPass(apo_iptm=0.85),
        offtarget_scorer=clean_offtarget,
    )
    assert result.candidates_tried == 2

def test_pool_exhausted_fails_loud():
    import pytest
    from pegrna_pipeline.hpc.recoded_tjpe_workflow import NoValidPegRNA
    with pytest.raises(NoValidPegRNA):
        run_recoded_tjpe(
            **minimal_inputs(),
            full_structure_scorer=FullStructureFake(pairs={(1, 100): 0.9}),  # always fails
            splice_scorer=SpliceFake(0.0),
            importance_scorer=AllPassImportance(),
            dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.80),
            offtarget_scorer=clean_offtarget,
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_recoded_tjpe_workflow.py -v`
Expected: FAIL with import error

- [ ] **Step 3: Write minimal implementation**

Implement `run_recoded_tjpe` with this control flow (pseudocode → real code):

```
apo_baseline = dock_scorer.apo_interface_iptm(nCas9)   # computed ONCE, reused
build ranked pool: for each PbsRttVariant -> PRIDICT2 rank -> sort desc
for candidate in pool:                      # ranked-pool walk
    candidates_tried += 1
    # --- cheap gates first (fail fast, no heavy calls) ---
    recode = recode_arm(arm_span, genome=..., frame=0, max_run=L)
    if recode is None: continue             # abandon: no synonymous solution
    # INV-0 with <=2 silent PAM edits (pegRNA PAM + ngRNA PAM); may raise -> caught -> continue
    assert_inv0(..., recode_span=arm_span, allowed_silent=[peg_pam_edit, ng_pam_edit],
                protein_check=translate_equal)
    if not passes_pbs2_gate(pbs2, repetitive=is_repetitive): continue
    # assemble in AUTHORITATIVE order: spacer+scaffold+pbs2+rtt+pbs1+polyu (Task 11b)
    full = assemble_tj_pegrna(spacer=.., pbs2=.., rtt=insert+recode.recoded, pbs1=.., polyu="UUU")
    # VERIFY predicted product == genome+cargo at S (catches a wrong PBS2 literal/RC convention)
    if not predicted_product_matches(full, genome, cargo, S): continue
    if not rnafold.passes_structure(full): continue   # cheap ViennaRNA pre-screen
    # --- expensive gates last ---
    if not passes_splice_gate(ref, alt, scorer=splice_scorer)[0]: continue
    # each out-of-arm PAM-break residue ALSO re-checked vs splice + importance (decision 11)
    if any(span_blocked(importance_scorer, [r]) for r in out_of_arm_pam_residues): continue
    if span_blocked(importance_scorer, footprint_residues): continue
    regions = {"spacer":.., "scaffold":.., "pbs2":.., "rtt":.., "pbs1":..}  # NB authoritative order
    if functional_region_pairing(full, regions=regions, scorer=full_structure_scorer).cross_region_pairs > 0:
        continue                            # RiNALMo 2D self-pairing (Task 10)
    if not passes_docking_gate(full, nCas9, scorer=dock_scorer, max_drop=0.15)[0]:
        continue                            # Boltz-2 scaffold-Cas9 docking, LAST + heaviest (Task 16)
    return RecodedTJPEResult(...)           # first survivor wins
# pool exhausted
widen search once (more synonymous variants / next nick); if still none:
raise NoValidPegRNA(diagnostic_table)
```

Also: before the pool loop, run the het-SNP / allele-specific search (Task 17) on both haplotypes and attach any allele-specific pairings to the result (advisory). Emit a diagnostic table (per candidate: which gate killed it) on failure — fail loud, never silent empty.

- [ ] **Step 4: Run tests to verify they pass**

Run: `PYTHONPATH=src python3 -m pytest tests/test_recoded_tjpe_workflow.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/recoded_tjpe_workflow.py tests/test_recoded_tjpe_workflow.py
git commit -m "feat(workflow): ranked-pool-walk recoded TJ-pegRNA orchestrator"
```

---

## Task 14: Reporting + provenance manifest

**Files:**
- Modify: `src/pegrna_pipeline/hpc/report_writer.py`, `src/pegrna_pipeline/hpc/annotator.py`
- Test: `tests/test_recoded_report.py`

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.report_writer import render_recoded_report

def test_report_includes_provenance_and_gates():
    text = render_recoded_report(sample_recoded_result(), manifest=sample_manifest())
    for token in ["genome_build", "GRCh38", "deletion_size", "0",
                  "splice", "alphafold", "rinalmo", "boltz2", "candidates_tried",
                  "pegRNA", "ngRNA", "PBS2", "PolyU", "recode_edits",
                  "allele_specific", "het_snp"]:
        assert token in text
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_recoded_report.py -v`
Expected: FAIL with "cannot import name 'render_recoded_report'"

- [ ] **Step 3: Write minimal implementation**

Add `render_recoded_report(result, manifest)` that emits: design (pegRNA = spacer/scaffold/PBS2/RTT/PBS1/PolyU per decision 0, ngRNA, insert, recoded arm), coordinates (S, A, B), deletion_size=0 (asserted), recode edit count + positions, **both silent PAM edits** (position/ref>alt/codon/which PAM), every gate result (pass/fail/value incl. RiNALMo cross-region pairs + Boltz-2 interface ipTM drop), advisory flags, **het-SNP advisories + allele-specific pairing options (with side/position-class/strength/phasing/other-allele off-target/efficiency penalty) + panel-wide allele-specific count**, and a provenance manifest (genome = T2T-KOLF2.1Jv1.1 + GRCh38 chain id, transcript version, PRIDICT version, DnaChisel/SpliceAI/Pangolin/AlphaFold/RiNALMo/RNAformer/Boltz-2 versions, codon table, git commit, seeds).

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_recoded_report.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/report_writer.py src/pegrna_pipeline/hpc/annotator.py tests/test_recoded_report.py
git commit -m "feat(report): recoded TJ-pegRNA report + provenance manifest"
```

---

## Task 15: Full-suite regression + lint + docs

**Files:**
- Modify: `README.md` (recoded mode + gate stack + optional extras)
- Modify: `docs/superpowers/specs/2026-06-09-scarless-tjpe-design.md` (link the recoded successor) or new spec

- [ ] **Step 1: Run the full suite**

Run: `PYTHONPATH=src python3 -m pytest tests/ -q`
Expected: all prior tests still pass (no regression on the 484 baseline) + new tests green

- [ ] **Step 2: Lint**

Run: `ruff check src/pegrna_pipeline/hpc tests`
Expected: All checks passed

- [ ] **Step 3: Update README + spec**

Document: the 7-stage pipeline (S0–S7 + S4b), the authoritative pegRNA structure (decision 0), the hard-gate stack (INV-0 ≤2 silent edits, splice, AlphaFold importance, conservation, RiNALMo 2D, Boltz-2 docking), the pool-walk control flow, the allele-specific feature + het-SNP advisory, the optional extras (`[recode]`, `[splice]`, `[importance]`, `[fullstructure]`, `[complex]`), and the wet-lab-validation caveats (sub-70 bp distance, PRIDICT OOD, RiNALMo/Boltz-2 thresholds-to-calibrate, microhomology threshold L, allele-specific availability ~5–10 %).

- [ ] **Step 4: Commit**

```bash
git add README.md docs/superpowers/specs tests
git commit -m "docs(recoded-tjpe): pipeline overview, gate stack, optional extras, caveats"
```

---

## Task 16: Scaffold–nCas9 docking gate (Boltz-2 Protocol + fake)

**Files:**
- Create: `src/pegrna_pipeline/hpc/complex_dock.py`
- Modify: `pyproject.toml` (`[complex]` extra: `boltz` note — GPU/HPC, MIT)
- Test: `tests/test_complex_dock.py`

Final + heaviest gate (decision 7b). Predicts the **full pegRNA + nCas9 complex** and flags if the scaffold↔Cas9 interface degrades vs an apo baseline. Boltz-2 is MIT (commercial-safe). Scorer injected via Protocol; apo baseline cached. Fallbacks: Chai-1 → Protenix → RoseTTAFold2NA.

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.complex_dock import passes_docking_gate, DockFake

def test_docking_gate_blocks_interface_degradation():
    # apo iptm 0.85; candidate collapses to 0.40 => degradation beyond tolerance
    scorer = DockFake(apo_iptm=0.85, candidate_iptm=0.40)
    ok, detail = passes_docking_gate("PEGRNA", "nCas9", scorer=scorer, max_drop=0.15)
    assert ok is False
    assert "interface" in detail

def test_docking_gate_passes_when_interface_preserved():
    scorer = DockFake(apo_iptm=0.85, candidate_iptm=0.82)
    ok, _ = passes_docking_gate("PEGRNA", "nCas9", scorer=scorer, max_drop=0.15)
    assert ok is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_complex_dock.py -v`
Expected: FAIL with import error

- [ ] **Step 3: Write minimal implementation**

```python
"""Scaffold-nCas9 docking gate. Boltz-2 (MIT) predicts the pegRNA+nCas9 complex;
we compare the scaffold<->Cas9 interface confidence (ipTM / interface-pLDDT, and
optionally a binding-affinity proxy) of the candidate against a cached apo
baseline. ComplexDockScorer is a Protocol so Boltz-2 (GPU/HPC) or a fake are
interchangeable; the apo baseline is computed once and reused across candidates.
Threshold has no published value for pegRNA-Cas9 => starts advisory, calibrate
on real data. Fallback scorers: Chai-1 (supports contact restraints) > Protenix
> RoseTTAFold2NA.
"""
from __future__ import annotations

from typing import Protocol, Tuple


class ComplexDockScorer(Protocol):
    def interface_iptm(self, pegrna: str, cas9: str) -> float: ...
    def apo_interface_iptm(self, cas9: str) -> float: ...


class DockFake:
    def __init__(self, apo_iptm: float, candidate_iptm: float) -> None:
        self._apo = apo_iptm
        self._cand = candidate_iptm

    def apo_interface_iptm(self, cas9: str) -> float:
        return self._apo

    def interface_iptm(self, pegrna: str, cas9: str) -> float:
        return self._cand


def passes_docking_gate(pegrna: str, cas9: str, *, scorer: ComplexDockScorer,
                        max_drop: float = 0.15) -> Tuple[bool, str]:
    apo = scorer.apo_interface_iptm(cas9)
    cand = scorer.interface_iptm(pegrna, cas9)
    drop = apo - cand
    if drop > max_drop:
        return False, (f"scaffold-Cas9 interface ipTM dropped {drop:.2f} "
                       f"(apo {apo:.2f} -> {cand:.2f} > {max_drop})")
    return True, f"interface preserved (apo {apo:.2f} -> {cand:.2f})"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_complex_dock.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/complex_dock.py tests/test_complex_dock.py pyproject.toml
git commit -m "feat(complex-dock): Boltz-2 scaffold-nCas9 interface gate (apo-vs-candidate)"
```

---

## Task 17: Het-SNP advisory + allele-specific pairing search (diploid)

**Files:**
- Create: `src/pegrna_pipeline/hpc/allele_specificity.py`
- Test: `tests/test_allele_specificity.py`

Pure module (decisions 12 + 13 + 15). Detects het-SNP overlap across functional elements with position-class and discrimination strength, and searches for allele-specific pairings. Consumes the diploid `GenomicContext` (Task 7). Advisory only — never gates.

**Allele-hunt nick-movement search (decision 15).** Add `enumerate_nick_movements(s_site, genome, cds_mask, max_shift, het_snps)` that, **only when invoked in `allele_hunt` mode**, yields candidate pegRNA nick positions shifted 5′/3′ of the default (S never moves). For each shift it returns: the discriminating het SNP captured (if any) + position class + strength; the re-templated span S↔nick partitioned into `recode_span` (CDS, via `cds_mask`) and `locked_span` (non-coding, byte-exact to the targeted allele); the **direction tradeoff** (5′-into-CDS recodable/SNP-poor vs 3′-into-3′UTR SNP-rich/locked+microhomology); and a **microhomology advisory** for the locked span. Candidates are ranked by **smallest movement distance** first (minimise locked span + RTT length). Outside `allele_hunt` mode this function is never called — the nick stays at its default and S↔nick is empty/CDS-only as today.

- [ ] **Step 1b: Write the failing test (nick-movement)**

```python
from pegrna_pipeline.hpc.allele_specificity import enumerate_nick_movements

def test_nick_movement_partitions_cds_vs_noncoding_and_locks():
    genome = "ACGT" * 50            # 200 nt
    s_site = 120                    # insert site (e.g. before stop)
    cds_mask = [i < 120 for i in range(200)]   # CDS upstream of S; 3'UTR downstream
    het = [{"pos": 130, "element": "spacer", "offset": 1, "haplotype": 1, "phased": True}]
    moves = enumerate_nick_movements(s_site=s_site, genome=genome,
                                     cds_mask=cds_mask, max_shift=20, het_snps=het)
    # a 3' move into the 3'UTR captures the het but the re-templated span is locked
    m = next(x for x in moves if x.direction == "3prime")
    assert m.locked_span is not None and m.recode_span is None
    assert m.microhomology_advisory is True          # non-coding cannot be recoded
    # a 5' move into CDS is recodable but captures no het here
    m5 = next(x for x in moves if x.direction == "5prime")
    assert m5.recode_span is not None
    # ranked by smallest movement first
    assert moves == sorted(moves, key=lambda x: x.shift_abs)
```

- [ ] **Step 1: Write the failing test**

```python
from pegrna_pipeline.hpc.allele_specificity import (
    classify_snp_position, discrimination_strength, find_allele_specific,
)

def test_position_class_pam_seed_distal():
    # spacer 20nt; PAM at 3' end; seed = PAM-proximal ~10nt
    assert classify_snp_position(element="spacer", offset=19, length=20) == "seed"   # PAM-proximal
    assert classify_snp_position(element="spacer", offset=1,  length=20) == "distal"
    assert classify_snp_position(element="pam", offset=1, length=3) == "pam"

def test_discrimination_strength_hierarchy():
    # PAM-destroying > seed > pbs1 > distal
    assert discrimination_strength("pam", "pegrna") > discrimination_strength("seed", "pegrna")
    assert discrimination_strength("seed", "pegrna") > discrimination_strength("pbs1", "pegrna")
    assert discrimination_strength("distal", "pegrna") == 0.0
    # ngRNA-side is valid but flagged with a byproduct note (not zero, not a scar)
    assert discrimination_strength("pam", "ngrna") > 0.0

def test_find_allele_specific_returns_option_when_discriminating_snp_present():
    het_snps = [{"element": "pegrna_pam", "offset": 1, "haplotype": 1, "phased": True}]
    opts = find_allele_specific(het_snps)
    assert len(opts) >= 1
    assert opts[0].allele_specific is True
    assert opts[0].side in ("pegrna", "ngrna")

def test_find_allele_specific_empty_when_no_discriminating_snp():
    het_snps = [{"element": "spacer", "offset": 2, "haplotype": 1, "phased": True}]  # distal
    assert find_allele_specific(het_snps) == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m pytest tests/test_allele_specificity.py -v`
Expected: FAIL with import error

- [ ] **Step 3: Write minimal implementation**

```python
"""Het-SNP advisory + allele-specific pairing search (diploid, advisory only).

A het SNP overlapping a DISCRIMINATING position lets a pegRNA-ngRNA pairing edit
only one allele => homogeneous all-het population (often preferable). We classify
the SNP's position (pam/seed/pbs1/distal), estimate discrimination strength, and
ALWAYS surface >=1 allele-specific pairing when available (decision 13). Both
pegRNA-side and ngRNA-side are valid; ngRNA-side carries a low-frequency PE2-like
byproduct NOTE on allele B (only an SSB is induced there; MMR reverts it; the
cargo cannot install without the 2nd nick), NOT a scar warning. Specificity can
only come from a pre-existing phased het (the 2 silent PAM edits write to both
alleles). The other allele is scored as an explicit off-target by the Task 6 gate.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

_SEED_NT = 10  # PAM-proximal window length that discriminates


def classify_snp_position(*, element: str, offset: int, length: int) -> str:
    if element in ("pam", "pegrna_pam", "ngrna_pam"):
        return "pam"
    if element in ("pbs1", "pbs2"):
        return "pbs1" if element == "pbs1" else "pbs2"
    # spacer: seed = PAM-proximal _SEED_NT (3' end of a Fw spacer)
    if element == "spacer":
        return "seed" if offset >= length - _SEED_NT else "distal"
    return "distal"


def discrimination_strength(position_class: str, side: str) -> float:
    base = {"pam": 1.0, "seed": 0.7, "pbs1": 0.4, "pbs2": 0.3, "distal": 0.0}.get(
        position_class, 0.0
    )
    # ngRNA-side discrimination is valid but modestly down-weighted (byproduct note)
    return base * (0.85 if side == "ngrna" else 1.0)


@dataclass(frozen=True)
class AlleleSpecificOption:
    allele_specific: bool
    side: str               # "pegrna" | "ngrna"
    position_class: str
    strength: float
    haplotype: int
    note: str


def find_allele_specific(het_snps: List[dict]) -> List[AlleleSpecificOption]:
    """Return allele-specific pairing options (possibly empty) from phased het SNPs."""
    opts: List[AlleleSpecificOption] = []
    for snp in het_snps:
        if not snp.get("phased", False):
            continue
        element = snp["element"]
        side = "ngrna" if element.startswith("ngrna") else "pegrna"
        pos = classify_snp_position(
            element=element, offset=snp["offset"], length=snp.get("length", 20)
        )
        strength = discrimination_strength(pos, side)
        if strength <= 0.0:
            continue  # distal => no real discrimination
        note = ("ngRNA-side: allele B sees only an SSB (nick) -> faithfully repaired, "
                "MMR reverts toward WT, cargo cannot install without 2nd nick; "
                "residual low-freq PE2 byproduct, raised under PE4/PE5") if side == "ngrna" \
               else "pegRNA-side: only allele A is nicked -> clean single-allele edit"
        opts.append(AlleleSpecificOption(
            allele_specific=True, side=side, position_class=pos,
            strength=strength, haplotype=snp["haplotype"], note=note,
        ))
    opts.sort(key=lambda o: -o.strength)
    return opts
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m pytest tests/test_allele_specificity.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add src/pegrna_pipeline/hpc/allele_specificity.py tests/test_allele_specificity.py
git commit -m "feat(allele-specificity): het-SNP advisory + allele-specific pairing search (diploid)"
```

---

## Self-review notes (gaps surfaced, addressed inline)

- **Spec coverage:** enumeration reshape (T11), assemble authoritative order (T11b), nick search + off-target incl. other-allele (T5/T6), ngRNA/PBS2 derivation (within T13 build), recode min-edits (T3/T4), INV-0 span-aware ≤2 edits (T2), splice/AlphaFold/RiNALMo gates (T8/T9/T10), Boltz-2 docking (T16), het-SNP + allele-specific (T17), pool-walk (T13), report+provenance (T14). The **2nd-PRIDICT ngRNA-derivation** sequence-surgery (keep scaffold+spacer+PBS2, drop Nick-RTgDNA+insert+HA) is exercised inside T13 — gets its own focused test during execution.
- **Pool-walk gate ordering (T13):** cheap→expensive — recode + INV-0 + PBS2 gate + ViennaRNA pre-screen FIRST; then splice; then AlphaFold importance (incl. the ≤2 out-of-arm PAM-break re-checks); then RiNALMo 2D (T10); then Boltz-2 docking (T16) LAST (cache apo baseline). First candidate passing all wins.
- **Deferred (own plans):** real T2T-KOLF diploid fetch + KOLF↔GRCh38 chain (T7 is type+fake only); real SpliceAI/Pangolin/AlphaFold/RiNALMo/RNAformer/Boltz-2 adapters (T8/9/10/16 ship Protocol+fake). Keeps CI light — correct seam.
- **Blockers:** B1/B2/B3 ALL RESOLVED (decisions 0/11/9). Remaining are execution-time *verifications*, not design unknowns: RC-vs-literal PBS2 convention (assert in T13); GRCh38↔KOLF chain availability (T7 follow-on); RiNALMo/RNAformer/Boltz-2 max-length + I/O + weight-staging.
- **Type consistency:** confirm `codon_expander.CodonVariant` field name (`.dna` assumed in T11) + `enumerate_variants` signature before T11; confirm `_gc_content` returns a fraction (0–1) not percent before T12; `assemble_tj_pegrna` kwargs are now `spacer/pbs2/rtt/pbs1/polyu/scaffold` (T11b) — all callers (T13) must match; Task 10 region keys must be `spacer/scaffold/pbs2/rtt/pbs1`.
- **Calibration parameters (all wet-lab-validate, surfaced in T15 docs + report):** microhomology threshold L (~10), RiNALMo pairing/accessibility threshold (~0.01–0.30), Boltz-2 interface-drop tolerance (~0.15), sub-70 bp nick distance, PRIDICT OOD note, allele-specific availability (~5–10 % of panel — compute exactly from the diploid genome).

---

## Execution handoff

Plan complete; **B1/B2/B3 resolved** (decisions 0, 11, 9). Tasks 1, 3, 5, 7 are pure/fixture-driven and can start immediately. Heavy adapters (T8/9/10/16) ship Protocol+fake here; real model + genome staging are follow-on plans.

Two execution options:
1. **Subagent-Driven (recommended)** — fresh subagent per task, two-stage review between tasks.
2. **Inline Execution** — batch with checkpoints via executing-plans.
