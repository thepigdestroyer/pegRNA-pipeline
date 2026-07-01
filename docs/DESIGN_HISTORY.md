# TJ-PE Recoded pegRNA Design Pipeline — Design & Decision History

> **Read this first.** This is the canonical, self-contained explanation of the
> recoded template-jumping prime-editing (TJ-PE) pegRNA designer: what it does,
> why each piece exists, why it was necessary, and how it is implemented. It also
> carries the full decision/rationale history so any human or agent can reconstruct
> *why the design is the way it is* without re-deriving it.
>
> If you only read one document about this pipeline, read this one. The code
> docstrings are the next layer of detail; the ADRs under `docs/decisions/` and the
> plans under `docs/superpowers/` are the primary sources this synthesises.

**Status (2026-06-18):** Feature-complete on branch `feat/recoded-tjpe` (PR #4).
Full unit suite green in CI (py3.10/3.11/3.12). Heavy ML gates have real adapters
that are validated on the Imperial CX3 HPC, not on a laptop. **Not yet
production-ready:** gate thresholds are literature-derived defaults, not wet-lab
calibrated, and the CX3 GPU smoke test is at 2/5 real gates passing (SpliceAI and
AlphaFold-importance confirmed live; RiNALMo and Boltz-2 blocked on environment).
See [§12 Status & limitations](#12-status--limitations).

---

## Table of contents

1. [The biological problem](#1-the-biological-problem)
2. [The core insight: recode, do not delete](#2-the-core-insight-recode-do-not-delete)
3. [The pegRNA architecture and strand convention](#3-the-pegrna-architecture-and-strand-convention)
4. [The end-to-end algorithm](#4-the-end-to-end-algorithm)
5. [The gate stack](#5-the-gate-stack)
6. [Genome model](#6-genome-model)
7. [The selection cassette](#7-the-selection-cassette)
8. [Tool-selection journey](#8-tool-selection-journey)
9. [Decision log](#9-decision-log)
10. [Adversarial reviews and the fixes they forced](#10-adversarial-reviews-and-the-fixes-they-forced)
11. [HPC deployment (Imperial CX3)](#11-hpc-deployment-imperial-cx3)
12. [Status & limitations](#12-status--limitations)
13. [Module / file map](#13-module--file-map)
14. [How to run](#14-how-to-run)
15. [Glossary](#15-glossary)
16. [Provenance](#16-provenance)

---

## 1. The biological problem

The goal is to insert a cargo (an epitope tag, a fluorophore, a selection cassette,
or a larger payload) at a precise site `S` in a human gene in induced pluripotent
stem cells (the lab uses KOLF2.1J), using **template-jumping prime editing
(TJ-PE / tjPE)**. TJ-PE is a prime-editing variant (Zheng et al. 2023, *Nat
Commun* 14:3369; the scarless tjPE basis is the s41467-023-39137-6 paper) that
uses a pegRNA whose 3' extension templates a new sequence and "jumps" to prime a
second nick made by an ngRNA on the opposite strand. The two nicks plus the
templated flap let you write an insert that is larger than a classic PE3 edit.

Two problems make naive TJ-PE unreliable:

1. **Microhomology-driven deletions.** When the templated flap shares a short
   perfect repeat (microhomology) with the genomic sequence near the second nick,
   the cell's MMEJ/SSA repair can collapse the intended product into a deletion.
   This is the dominant failure mode for larger inserts.

2. **Silent collateral damage.** An insert or its flanking edits can disrupt a
   splice site, fold into a structure that occludes the spacer or scaffold, break
   the pegRNA's ability to dock nCas9, or land in a structurally critical protein
   residue. None of these show up as an obvious "the edit didn't happen"; they show
   up as a dead cell line, a mis-spliced transcript, or a non-functional protein.

The pipeline exists to design a pegRNA + ngRNA pair that performs the insert
**and** is defended against every one of those failure modes before any DNA is
ordered.

The recoded designer is built on top of an existing pegRNA pipeline (PRIDICT2
submission, PE7, prime-amplification workflows, off-target scoring, report
writing). The recoded layer is the `*_recoded*`, gate, and orientation modules
under `src/pegrna_pipeline/hpc/`.

---

## 2. The core insight: recode, do not delete

The naming "scarless" originally implied deleting the microhomology. The redesign
replaced that with a stronger idea:

**Instead of deleting the microhomology-prone genomic window between the insert
site `S` and the ngRNA nick, re-template that window with a synonymously *recoded*
copy of itself.** The recoded arm:

- encodes the **same protein** (every coding change is silent), and
- **breaks the perfect microhomology run** that drives the deletion, by spending
  the minimum number of synonymous codon swaps needed to get the longest perfect
  run below the efficient-homology onset.

This is why `RecodedTJPEResult.deletion_size == 0`: the arm is not deleted, it is
length-preserved and recoded. The net edit is "insert cargo at `S` + a few silent
codon changes in the immediately adjacent CDS," which is invisible at the protein
level and resistant to MMEJ collapse.

The arm can sit 3' of `S` (default, `arm_side="3prime"`) or 5' of `S`
(`arm_side="5prime"`, the C-terminal "recode-into-CDS" geometry used for
C-terminal tags). Where the arm crosses a coding/non-coding boundary, only the CDS
run is recoded; any non-coding run is re-templated **byte-exact** to the targeted
allele (Decision 15), because synonymous recoding is undefined outside a CDS.

This insight is implemented in `recode_arm.py` (DnaChisel-backed synonymous
recoder) and `recoded_tjpe_workflow._partition_arm` (splits the arm into CDS /
non-coding runs and applies recode-or-lock per run).

---

## 3. The pegRNA architecture and strand convention

This section took the most iteration to get right and is the most common source of
silent, catastrophic bugs. It is **LOCKED** (Decision B1 and the strand
convention). See ADR `0004-pegrna-architecture-and-strand-convention.md`.

### 3.1 Authoritative pegRNA structure

A TJ-pegRNA, written **5'→3'**, is:

```
Spacer · Scaffold · PBS2 · RTT · PBS1 · PolyU
```

| Segment   | Role                                                                                   |
|-----------|----------------------------------------------------------------------------------------|
| Spacer    | Targets and nicks at site `S`. **Fixed** by `S`; never enumerated.                     |
| Scaffold  | Standard sgRNA scaffold (`STANDARD_SGRNA_SCAFFOLD`). Binds nCas9.                       |
| PBS2      | Occupies the homology-arm slot. The captured-strand window 5' of the **ngRNA** nick.   |
| RTT       | Reverse-transcription template: the templated flap (cargo + recoded/locked arm).       |
| PBS1      | Primer-binding site for the pegRNA's own nick: primes off the 3'-OH the spacer creates.|
| PolyU     | 2-3 U 3' tail (La/PE7 motif) that protects the extension.                              |

Assembled by `assemble_tj_pegrna_v2(spacer, pbs2, rtt, pbs1, polyu, scaffold)` in
`tjpe_workflow.py`. PolyU defaults to `"UUU"`.

**Why PBS2 is in the homology-arm slot and not where intuition puts it:** TJ-PE
template-jumps. The pegRNA first primes its own nick (PBS1), reverse-transcribes
the RTT, then the new 3' flap must anneal to the second (ngRNA) nick. PBS2 is the
sequence that captures that second strand. Putting it before the RTT in the 5'→3'
layout is what makes the jump work.

### 3.2 Strand convention (the bookkeeping that bites)

- **Everything is stored 5'→3'.** Always.
- Antiparallel relationships are expressed with `reverse_complement`, never with a
  bare `reverse` or a bare `complement`. Those are three distinct operations and
  conflating them is the #1 latent bug class here.
- The **RTT tail** on the pegRNA is the **reverse complement of the sense flap**
  for a reverse-targeting pegRNA (`orient_rtt_tail`).
- **PBS2** is the **literal/identity** copy of the captured strand for a
  reverse-targeting pegRNA, not its reverse complement (`orient_pbs2`). This is the
  "complement of a complement collapses to identity" reasoning. Two independent
  reviewers proposed RC here; it was deliberately kept as identity-for-Rv. See
  [§10](#10-adversarial-reviews-and-the-fixes-they-forced).
- **PBS1** is derived from the **genome** at the pegRNA nick:
  `reverse_complement(genome[S-L : S])` for a reverse-targeting pegRNA, identity for
  forward (`_derive_pbs1`). It is **not** a slice of the spacer. Deriving PBS1 from
  the spacer produced a pegRNA that could not prime, a dead design that looked fine
  on paper. This was a high-severity review fix.
- The default `pegrna_strand="Rv"` (reverse-targeting) because the C-terminal
  recode-into-CDS geometry needs the pegRNA on the strand that lets it recode the
  codons 5' of `S`.

Every assembled candidate is run through `simulate_tprt` (`tprt_sim.py`), which
performs the reverse-transcription in silico from the *oriented* RTT tail and
asserts it reconstructs the intended sense allele. If any
complement/reverse/reverse_complement is wrong, `TPRTMismatch` fires loudly and the
candidate is rejected. This is the safety net that makes the strand convention
auditable rather than faith-based. (It is a self-consistency check, not a
wet-lab-validated strand proof; see limitations.)

---

## 4. The end-to-end algorithm

Implemented in `recoded_tjpe_workflow.run_recoded_tjpe`. The shape is a
**ranked-pool-walk**: build a deterministic candidate pool once, then walk it
running gates cheap→expensive, and the **first candidate that passes every gate
wins**. The insert site `S` **never moves**. If the pool is exhausted,
`NoValidPegRNA` is raised carrying a per-candidate diagnostic of which gate killed
each candidate (fail-loud, never fail-silent or return a junk design).

Order of operations:

1. **Allele options (advisory).** If `allele_hunt=True`, `find_allele_specific`
   surfaces het-SNP-based allele-specific pairing options. Advisory only; wrapped in
   a try/except so a malformed SNP record can never crash a design run.
2. **Nick enumeration + off-target screen.** `enumerate_nick_sites` proposes ngRNA
   nicks on the opposite strand within `[ngrna_distance_min, ngrna_distance_max]` =
   40-70 bp of `S` (smaller preferred; **not** auto-calibrated). `screen_nick_candidates`
   removes off-target-unsafe nicks (including against the other allele).
3. **Candidate enumeration.** `enumerate_pbs_rtt_variants` = (synonymous RTT
   variants) × (PBS lengths, default 13/15/17). The spacer is fixed.
4. **Apo docking baseline.** Warm the Boltz-2 scaffold+nCas9 apo baseline once
   (cached for reuse). Guarded: if it cannot compute, raise `NoValidPegRNA` before
   wasting the walk.

Then for each (nick × variant), in this exact order:

5. **Arm partition + recode** (`_partition_arm`): split `genome[lo:hi)` into CDS /
   non-coding runs; recode CDS runs synonymously with longest perfect run ≤ 8 bp;
   keep non-coding runs byte-exact. Reject if a CDS run has no synonymous solution
   or the arm is non-monotone (crosses more than one coding boundary).
6. **Assemble the edit** (`_assemble_edit`): cargo always at `S`; arm on the chosen
   side; produce the sense flap and the predicted edited allele.
7. **Orient** the RTT tail and PBS2 per strand.
8. **PBS2 gate** (`passes_pbs2_gate`): length 11-17, GC 40-70%, Tm 37-42 °C
   (Biopython nearest-neighbour), not repetitive. Length and Tm are coupled, so the
   feasible PBS2 set tends to collapse toward 14-17 nt; that is expected.
9. **INV-0** (`assert_inv0`): the only thing that may differ between reference and
   product is the cargo at `S` plus protein-silent recode-span edits plus byte-exact
   locked-span re-templating. The CDS recode is verified protein-silent
   (`_arm_protein_equal`); the non-coding locked span is verified byte-exact.
10. **TPRT strand check** (`simulate_tprt`): reconstruct the allele from the oriented
    RTT; reject on `TPRTMismatch`.
11. **Derive PBS1** from the genome (not the spacer).
12. **Assemble** the full pegRNA.
13. **Splice gate** (`passes_splice_gate`, SpliceAI/Pangolin): reject if predicted
    splice Δ > 0.2. Guarded + fail-loud.
14. **AlphaFold importance gate** (`span_blocked`): only runs if the recode footprint
    maps to protein residues; rejects if any recoded residue is structurally critical
    (buried/confident/high-contact). Guarded.
15. **RiNALMo full-pegRNA 2D structure gate** (`functional_region_pairing`): reject if
    any functional region (spacer/pbs2/rtt/pbs1) cross-pairs with another. Guarded.
16. **Boltz-2 docking gate LAST** (`passes_docking_gate`): reject if the
    scaffold+nCas9 interface ipTM drops too far vs the apo baseline. Most expensive,
    so it runs last. Guarded.

First survivor returns a `RecodedTJPEResult` with the full pegRNA, its segments, the
ngRNA protospacer, `S`, the nick position, recode edit count, candidates tried, and
any advisory warnings (e.g. residual microhomology above target, or a non-coding
locked span that could not be recoded).

**Why this order:** every gate is a filter, and the cheap deterministic filters
(recode, INV-0, PBS2, strand) eliminate most of the pool before the expensive GPU
models run. The expensive models are ordered by cost (SpliceAI < AlphaFold lookup <
RiNALMo < Boltz-2). Every expensive gate is individually guarded so one scorer
crashing skips that candidate with a diagnostic instead of aborting the whole walk.

---

## 5. The gate stack

Each gate is described as: **what it checks · why it is necessary · how it is
implemented · tool · threshold + citation.** Thresholds live in
`gate_thresholds.py` (`THRESHOLDS`) with a citation string per value. They are
literature-derived defaults, overridable by a wet-lab calibration library.

### 5.1 Microhomology recode gate
- **What:** longest perfect repeat between the templated arm and flanking genome.
- **Why:** perfect microhomology runs drive MMEJ/SSA deletion of the insert.
- **How:** `microhomology.py` measures runs; `recode_arm.py` (DnaChisel) spends the
  minimum synonymous edits to push the longest run below the ceiling, frame-aware.
- **Tool:** DnaChisel.
- **Threshold:** `microhomology_max_run = 8` bp. Rationale: RecA/RAD51 minimal
  efficient processing segment ≈ 8 bp (Shen & Huang 1986; Hsieh et al. 1992); MMEJ
  efficiency onset ≈ 6-8 bp (Sfeir & Symington 2015). The recoder minimises the
  actual run below this ceiling; smaller preferred.

### 5.2 PBS2 gate
- **What:** PBS2 length / GC / Tm / repetitiveness.
- **Why:** PBS2 must anneal the second strand stably but not so stably it mis-primes;
  repetitive PBS2 mis-targets.
- **How:** `passes_pbs2_gate`; Tm via Biopython `Tm_NN`.
- **Threshold:** length 11-17 nt, GC 40-70%, Tm 37-42 °C. Length and Tm coupling
  collapses the feasible set toward 14-17 nt (expected, not a bug).

### 5.3 INV-0 (invariant zero) gate
- **What:** the product differs from reference by exactly: cargo at `S` +
  protein-silent recode + byte-exact locked re-template. At most two silent PAM
  edits (pegRNA PAM + ngRNA PAM) are allowed elsewhere.
- **Why:** this is the design's correctness contract. Anything else is collateral
  damage.
- **How:** `scarless.assert_inv0`, span-aware: `recode_span` changes checked
  protein-silent via a real translate (`_arm_protein_equal`, not a stub);
  `locked_span` checked byte-exact to the targeted allele.
- **Decision:** B2-ii. Each out-of-arm break is itself re-checked against splice +
  AlphaFold. See ADR `0001`.

### 5.4 Splice gate
- **What:** does the edit create/destroy a splice site near the locus?
- **Why:** a silent-at-protein-level edit can still mis-splice the transcript.
- **How:** `splice_screen.passes_splice_gate`; fail-loud (a scorer error rejects the
  candidate, it does not pass it through).
- **Tool:** SpliceAI (Pangolin alternative).
- **Threshold:** `splice_max_delta = 0.2`. SpliceAI Δ cutoffs are 0.2 (high recall),
  0.5 (recommended), 0.8 (high precision) — Jaganathan et al. 2019, *Cell* 176:535.
  As a **safety** gate it uses 0.2 to maximise detection.

### 5.5 AlphaFold structural-importance gate
- **What:** is a recoded residue structurally critical (well-folded, buried,
  densely packed)?
- **Why:** synonymous DNA edits are protein-silent by construction, so this gate is
  really about the *insert's* protein-level footprint and about not placing a
  C-terminal tag into a buried/critical residue.
- **How:** `structure_importance.span_blocked` over the recode→residue footprint,
  using AlphaFold-DB structures parsed with biotite (per-residue pLDDT, relative
  solvent accessibility, heavy-atom contact number). Only runs when the footprint
  maps to residues; a purely non-coding arm changes no residue and skips it.
- **Tool:** AlphaFold Protein Structure Database + biotite. (Fetches `*.cif`; the
  AFDB API serves v6 `.bcif`/`.cif`; the v4→v6 resolution and a biotite `b_factor`
  bug were fixed on CX3.)
- **Thresholds:** pLDDT ≥ 70 (lower edge of "confident", Jumper et al. 2021); RSA
  buried < 0.25 (Levy 2010; Tien et al. 2013 max-ASA normalisation); > 18 heavy-atom
  contacts within 5 Å marks a packed core residue (Yeh et al. 2014; Shih et al.
  2012).

### 5.6 RiNALMo full-pegRNA 2D-structure gate
- **What:** does the assembled pegRNA mis-fold so a functional region (spacer, PBS2,
  RTT, PBS1) base-pairs into another region?
- **Why:** a spacer occluded by intramolecular pairing will not load nCas9; an
  occluded PBS will not prime.
- **How:** `full_structure.functional_region_pairing` computes the base-pair
  probability matrix over the whole pegRNA and rejects any cross-region pairing among
  named functional regions (scaffold excluded by name from cross-region scoring).
- **Tool:** RiNALMo (Penić et al. 2024) as primary; RNAformer as CPU fallback;
  ViennaRNA as the cheap pre-screen.
- **Threshold:** base-pair-probability decision boundary 0.5; min accessibility 0.5
  (a region is flagged occluded when mean unpaired probability < 0.5). 0.5 chosen to
  avoid over-blocking on low-probability pairs.

### 5.7 Boltz-2 scaffold–nCas9 docking gate
- **What:** does the pegRNA extension degrade the scaffold's ability to dock nCas9?
- **Why:** the RTT/PBS extension can fold back and disrupt the scaffold–protein
  interface; if nCas9 cannot bind, nothing happens.
- **How:** `complex_dock.passes_docking_gate` compares the candidate
  scaffold+extension+nCas9 interface ipTM against the **apo** scaffold+nCas9 baseline
  (Decision: apo = scaffold + Cas9, not the scaffold monomer — a review fix).
- **Tool:** Boltz-2 (AlphaFold-Multimer-class ipTM).
- **Thresholds:** ipTM floor 0.6 (below = likely-incorrect interface, Evans et al.
  2022 / Abramson et al. 2024); max drop 0.2 vs apo (one full confidence band).

---

## 6. Genome model

Decision B3-C, refined by ADR 0007. The truth is the KOLF2.1J iPSC line (male/XY), and
the ML gates trained on GRCh38 (SpliceAI, AlphaFold) operate in GRCh38 coordinates. The
pipeline stays in **GRCh38 coordinates** and materialises the real diploid
**coordinate-safely (scope A)**: KOLF2.1J phased substitutions are applied onto the
GRCh38 region sequence per haplotype (`variant_apply.apply_phased_variants`, wired
through `genome_ingest.build_genomic_context`). Same-length substitutions preserve the
coordinate frame, so `to_grch38` stays a valid identity and the GRCh38-trained gates are
untouched; a ref-base mismatch is fail-loud; indels are advisory only.

A clean diploid KOLF2.1J assembly now exists (NCBI GCA_058280895.1 `pri` +
GCA_058280805.1 `alt`, NHGRI, 2026-06-17), so the old availability blocker is gone.
Making it the *truth* genome is **scope B** (separate follow-on): it needs a validated
KOLF<->GRCh38 liftover chain and a GENCODE->KOLF annotation projection, because feeding
KOLF coordinates to the GRCh38-trained gates without a chain silently corrupts scores.
Scope B is not a drop-in data change.

The genome context provides the local sequence window, the CDS mask (which bases are
coding, used to partition the arm), and splice-site context. Het-SNPs are surfaced as
advisory allele-specific options and, in `allele_hunt` mode, can move the **nick**
(never the insert site `S`). Non-coding bases are re-templated byte-exact to the
**targeted allele** so the design does not silently overwrite the other allele.

---

## 7. The selection cassette

For selectable inserts the cargo carries a **2A-peptide–Zeocin-resistance** cassette
downstream of the gene of interest (the CTCF C-terminal design used
`GGGGS-FLAG-P2A-ZeoR-STOP`; the general principle uses a 2A-ZeoR cassette). The
mechanism is a biological self-selection filter: any product that frameshifts the GOI
also frameshifts ZeoR, so it dies under Zeocin. What survives selection is therefore
limited to in-frame products. The residual error classes the in-silico gates still
have to catch are in-frame missense, in-frame indels, and symmetric exon-skips,
which is exactly what the splice, AlphaFold, and INV-0 gates target. Promoter
elements (`T7` / `T7-rev`) are insertable via `insert_elements.py`.

---

## 8. Tool-selection journey

This records *why* each external tool was chosen, including the ones that were
rejected, so the choices are not re-litigated from scratch. See ADR
`0005-gate-and-tool-stack.md`.

- **RNA 2D structure → RiNALMo.** Candidates considered: ViennaRNA (fast, but a
  thermodynamic MFE model, not learned and weak on long multi-domain RNAs like a full
  pegRNA), NuFold and RibonanzaNet (3D / chemical-probing oriented, heavier and less
  suited to a pairwise base-pair-probability read), RNAformer. RiNALMo (a
  language-model-based RNA structure predictor, Penić et al. 2024) won for full-pegRNA
  base-pair-probability matrices. **RNAformer** is kept as a CPU fallback for when no
  GPU is available; **ViennaRNA** is the cheap pre-screen.
- **Recode → DnaChisel.** Purpose-built synonymous-recoding optimiser with codon and
  constraint support; no reason to hand-roll codon tables.
- **Splice → SpliceAI (Pangolin alternative).** SpliceAI is the de-facto standard with
  published, well-characterised Δ cutoffs, which is what let the threshold be set
  defensibly rather than arbitrarily.
- **Protein importance → AlphaFold-DB + biotite.** Pre-computed structures avoid
  folding on the fly; biotite gives per-residue pLDDT/RSA/contacts cheaply.
- **Scaffold–nCas9 docking → Boltz-2.** Open, runnable complex-structure predictor
  giving an AF-Multimer-class ipTM for the interface comparison.

---

## 9. Decision log

The locked decisions. These are the constraints the implementation must satisfy;
violating one is a bug, not a tradeoff. Primary sources: the plan
`docs/superpowers/plans/2026-06-18-recoded-tjpe-pegrna-redesign.md`, ADRs `0001`,
`0003`, `0004`, `0005`, `0006`, and the session memory log.

| #  | Decision | Rationale |
|----|----------|-----------|
| B1 | Authoritative pegRNA structure `Spacer·Scaffold·PBS2·RTT·PBS1·PolyU` (5'→3'); PBS2 in homology-arm slot; PolyU = La/PE7 tail. | Makes the template-jump mechanically correct; locks segment order for assembly + scoring. |
| —  | Strand convention: store everything 5'→3'; antiparallel = `reverse_complement`; RTT tail = RC of sense flap; PBS2 = identity for Rv; PBS1 = RC of genome window (not spacer slice). | `complement`/`reverse`/`reverse_complement` confusion is the #1 silent-bug class; `simulate_tprt` enforces it. |
| B2-ii | INV-0 allows cargo at `S` + protein-silent recode + byte-exact lock + ≤2 silent PAM edits; each out-of-arm break re-checked vs splice+AlphaFold. | Defines the correctness contract precisely. |
| B3-C | Genome = T2T diploid KOLF2.1J as truth, GRCh38 coords for ML gates (wired as Ensembl-GRCh38 + KOLF SNP overlay). | KOLF2.1J is the lab line; ML models are GRCh38-trained; T2T FASTA swap is a later data change. |
| 1  | Insert site `S` is fixed; the pegRNA spacer is fixed; only the nick + RTT variants + PBS length are enumerated. | `S` is the experimental intent; moving it changes the experiment. |
| —  | ngRNA nick distance 40-70 bp, smaller preferred, **not** auto-calibrated. | The one gate left as an experimental-design knob, not a threshold. |
| —  | Recode = minimum synonymous edits to break microhomology, frame-aware. | Fewer edits = less risk, easier synthesis. |
| 8  | PBS2 gate length 11-17 / GC 40-70 / Tm 37-42; coupling collapses to ~14-17 nt. | Stable but specific second-strand capture. |
| 15 | Arm partitioned by CDS mask: recode CDS runs, keep non-coding runs byte-exact to the targeted allele; reject non-monotone (>2-run) arms. | Synonymous recoding is undefined outside a CDS; do not overwrite the other allele. |
| —  | Selection = 2A-ZeoR cassette downstream of GOI; frameshifts self-select out. | Biological pre-filter; in-silico gates handle the in-frame residue. |
| —  | Gate thresholds literature-derived, each with a citation; wet-lab override allowed. | Defensible defaults, not arbitrary placeholders. |
| —  | `allele_hunt` advisory; het-SNPs may move the nick, never `S`. | Allele-specific editing without changing the experimental target. |
| —  | Fail-loud: exhausting the pool raises `NoValidPegRNA` with diagnostics; never return a junk design. | A wrong "successful" pegRNA is worse than an honest failure. |
| —  | Heavy ML scorers are **injected** (fakes in CI, real adapters on CX3); the workflow imports none at module load. | Keeps unit tests fast/deterministic; isolates GPU deps to the cluster. |

---

## 10. Adversarial reviews and the fixes they forced

The design and implementation were each put through adversarial review, which is why
several non-obvious choices look the way they do.

**Design-stage:** multi-model adversarial panels (and a multi-model consensus reject
on an early plan) drove the recode-not-delete insight, the INV-0 contract, the
genome-model choice, and the strand convention being written down explicitly.

**Implementation-stage:** a 5-subagent parallel review (code, security, adversarial,
refactor, user-stand-in) produced ~25 findings; **13 high-severity were fixed and
committed (CI stayed green).** The ones worth knowing about:

1. **PBS1 was a slice of the spacer** → re-derived from the genome at the pegRNA nick
   (`_derive_pbs1`). The old version produced a pegRNA that could not prime: a dead
   design that passed every other check.
2. **INV-0 protein check was `lambda: True`** (a stub that always passed) → replaced
   with a real translate comparison (`_arm_protein_equal`).
3. **Minus-strand coordinates were not mirror-flipped** → fixed coordinate handling
   for reverse-strand targets.
4. **Recode ignored reading frame** → frame-aware recoding (`_frame_at`, `frame=`
   threaded through `recode_arm`).
5. **Splice gate failed open** (a scorer error let the candidate pass) → fail-loud
   (error rejects the candidate).
6. **Boltz apo baseline was the scaffold monomer** → apo = scaffold + Cas9, the
   correct reference interface.
7. **RiNALMo accessibility used a partial matrix** → full base-pair-probability matrix.
8. **Contacts were not heavy-atom** → heavy-atom contacts within 5 Å.
9. **Residue-footprint mapping** (genomic recode span → protein residue indices) added
   so the AlphaFold gate scores the right residues.
10. **Per-gate crash guards** so one scorer failure skips a candidate with a diagnostic
    instead of aborting the walk.
11. **RNAformer checkpoint resolution** by env/glob, not a fixed filename (CX3 smoke
    caught this).

**Deliberately kept against review pressure:** `orient_pbs2` stays
identity-for-Rv (two reviewers suggested reverse-complement); the
complement-of-complement reasoning plus the `simulate_tprt` self-consistency check
justify it. `simulate_tprt` was honestly downgraded in its claims to a
self-consistency check, because a true strand proof needs an experimentally validated
golden pegRNA, which the project does not yet have.

---

## 11. HPC deployment (Imperial CX3)

The heavy ML gates cannot run on the dev laptop (no GPU), so they are validated on
the Imperial CX3 HPC.

- **Connect:** `source ~/.cx3-credentials; sshpass -p "$HPC_PASS" ssh -o
  PreferredAuthentications=password -o PubkeyAuthentication=no -o
  StrictHostKeyChecking=accept-new mj1824@login.cx3.hpc.imperial.ac.uk`. Use direct
  sshpass; the hpc-connect wrapper's `exec` drops non-interactive channels. Login node
  `login-bi`; GPU queue `v1_gpu72`; `qsub` at `/opt/pbs/bin/qsub`.
- **Env:** conda env `pegrna-recoded` (py3.11) at `~/miniforge3/envs/pegrna-recoded`
  with torch / boltz / spliceai / tensorflow-cpu. `scripts/hpc/environment.yml`.
- **Compute nodes are OFFLINE.** Everything (model weights, AFDB structures) must be
  pre-staged on the **login** node. `scripts/hpc/stage_models.sh` downloads RiNALMo
  (~7.3 GB, Zenodo 15043668), RNAformer (~131 MB), and AFDB cif files.
- **`$EPHEMERAL` is empty in a non-login SSH channel** → use `bash -lc` or hardcode
  `/rds/general/user/mj1824/ephemeral/pegrna-data`, or you get false "missing file"
  failures.
- **Validated live so far:** full unit suite 595 passed on-cluster; end-to-end
  fake-mode run works; **2 real gates confirmed (SpliceAI, AlphaFold-importance)** with
  real metrics. GPU smoke (`scripts/hpc/smoke_models.py`, job in `v1_gpu72`) is at 2/5.
- **Blocked:** RiNALMo `flash-attn` will not build on the GPU node (`nvcc` missing,
  needs a CUDA-devel module or a prebuilt wheel); Boltz-2 produces no confidence JSON
  (weights/offline). These are the remaining two of five smoke gates.

Runbook: `scripts/hpc/README.md`.

---

## 12. Status & limitations

**What is done:** all 18 redesign tasks implemented and committed on
`feat/recoded-tjpe` (PR #4); full unit suite green in CI across py3.10/3.11/3.12; real
adapters written for every heavy gate; 13 high-severity adversarial findings fixed;
literature thresholds set and cited; CX3 env + staging done; 2/5 real gates validated
live.

**What is NOT done / do not over-trust:**

- **Thresholds are literature-derived defaults, not wet-lab calibrated.** A
  calibration library can and should override them. They are defensible, not proven.
- **RiNALMo and Boltz-2 gates are not yet validated live** (2/5 smoke). Until the
  GPU environment is fixed, those two gates run as fakes off-cluster.
- **`simulate_tprt` is a self-consistency check, not an experimental strand proof.**
  It catches orientation bugs; it does not prove the design edits in a cell.
- **Genome is GRCh38-coordinate with a coordinate-safe materialised KOLF2.1J diploid**
  (scope A, substitutions only; ADR 0007). Indels are advisory, and full KOLF-as-truth
  (assembly FASTA + KOLF↔GRCh38 chain + annotation projection) is the pending scope B.
- The pipeline is **not production-ready**: do not order DNA off its output without
  wet-lab threshold calibration and a 5/5 smoke pass. This was stated explicitly to
  the user; no false green light.

---

## 13. Module / file map

Recoded layer (`src/pegrna_pipeline/hpc/`):

| File | Responsibility |
|------|----------------|
| `recoded_tjpe_workflow.py` | Orchestrator `run_recoded_tjpe`, candidate enumeration, gate order, ranked-pool-walk, `NoValidPegRNA`. |
| `recoded_driver.py` | CLI driver for `pegrna design-recoded`. |
| `recode_arm.py` | DnaChisel synonymous recoder (frame-aware, microhomology-minimising). |
| `microhomology.py` | Longest-perfect-run measurement. |
| `nick_search.py` | ngRNA nick enumeration (40-70 bp) + off-target screen. |
| `genome_ingest.py` / `genome_context.py` | Ensembl-GRCh38 ingest + KOLF SNP overlay; CDS mask; splice context. |
| `splice_screen.py` | SpliceAI/Pangolin splice gate (fail-loud). |
| `structure_importance.py` | AlphaFold-DB + biotite residue importance (pLDDT/RSA/contacts), `span_blocked`. |
| `full_structure.py` | RiNALMo full-pegRNA 2D base-pair-probability gate, `functional_region_pairing`. |
| `complex_dock.py` | Boltz-2 scaffold–nCas9 docking gate, apo baseline. |
| `allele_specificity.py` | Het-SNP allele-specific options (`allele_hunt`). |
| `strand_orientation.py` | `orient_rtt_tail`, `orient_pbs2` (the strand convention in code). |
| `tprt_sim.py` | `simulate_tprt`, `TPRTMismatch` (strand self-consistency net). |
| `insert_elements.py` | T7 / T7-rev promoter insertion. |
| `scarless.py` | `assert_inv0` (span-aware INV-0 contract). |
| `gate_thresholds.py` | `THRESHOLDS` + per-value citations. |
| `scorer_factory.py` | `build_scorers`: real-adapter-with-fake-fallback, per gate. |
| `tjpe_workflow.py` | `assemble_tj_pegrna_v2`, `STANDARD_SGRNA_SCAFFOLD` (shared base). |

Existing base pipeline (selected): `pridict_submitter.py` / `pridict_parser.py`
(PRIDICT2), `pe7_workflow.py`, `pa_workflow.py`, `offtarget_scorer.py`,
`pair_scorer.py`, `pbs.py`, `rnafold.py`, `report_writer.py`, `sequence_fetcher.py`,
`driver.py`, `annotator.py`, `codon_expander.py`.

Docs: this file; `docs/decisions/*.md` (ADRs); `docs/superpowers/plans/*.md`;
`docs/superpowers/specs/*.md`. HPC: `scripts/hpc/`.

---

## 14. How to run

```bash
# Unit tests (fast, fakes for all heavy gates)
ruff check src tests && black --check src tests   # CI runs both, line-length 100
pytest -q

# Design CLI (recoded TJ-PE)
pegrna design-recoded --help
```

`run_recoded_tjpe` takes the spacer, the RTT amino-acid intent, the cargo, the genome
window, `site_s`, the CDS mask, PBS lengths, nick distance bounds, `pegrna_strand`,
`arm_side`, and the injected scorers. In CI those scorers are fakes; on CX3
`scorer_factory.build_scorers` wires the real adapters (with a fake fallback per gate
when a model is unavailable). For a real run, do it on CX3 after a 5/5 smoke pass.

---

## 15. Glossary

- **TJ-PE / tjPE** — template-jumping prime editing.
- **pegRNA** — prime-editing guide RNA.
- **ngRNA** — nicking guide RNA (second nick, opposite strand).
- **PBS** — primer-binding site. PBS1 = pegRNA's own nick; PBS2 = captures the ngRNA-nicked strand.
- **RTT** — reverse-transcription template (the templated flap).
- **INV-0** — the design's "invariant zero": the exact allowed difference between reference and product.
- **`S`** — the fixed insert site.
- **Arm** — the genomic window between `S` and the nick that is recoded (CDS) or locked (non-coding).
- **Recode** — synonymous codon substitution that preserves protein but breaks microhomology.
- **MMEJ / SSA** — microhomology-mediated end joining / single-strand annealing (the deletion pathways).
- **ipTM** — interface predicted TM-score (docking confidence).
- **RSA** — relative solvent accessibility.

---

## 16. Provenance

This pipeline was designed and implemented across several agent sessions with
adversarial review at the design and implementation stages. The authoritative
sources, in decreasing order of detail, are: this document → the code docstrings →
`docs/decisions/` ADRs → `docs/superpowers/plans/` and `specs/`. The decision/
rationale narrative is also mirrored in the agent memory log
(`memory/observations.json` and the per-session `global.json` entries) so it survives
across sessions and machines. When you change a locked decision, update the relevant
ADR, this file's [§9](#9-decision-log), and add an observation to memory.
