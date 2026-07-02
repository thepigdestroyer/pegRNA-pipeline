# Features

Technical capability + roadmap reference for the recoded template-jumping
prime-editing (TJ-PE) pegRNA design pipeline. Audience: a computational /
ML-biology reader. For the design rationale and decision history see
[`DESIGN_HISTORY.md`](DESIGN_HISTORY.md); for a non-technical walk-through see
[`USER_GUIDE.md`](USER_GUIDE.md).

## What it produces

Given a target locus, an insert site `S`, and a cargo, the pipeline designs a
**recoded TJ-PE pegRNA** (plus its ngRNA) that knocks the cargo in at `S` while
defending the edit against a stack of biological failure modes. Output is a
human-readable report plus a component-breakdown CSV (see
[Output](#output)).

The construct is emitted 5′→3′ as **Spacer · Scaffold · PBS2 · RTT · PBS1 · PolyU**,
where the RTT is decomposed into the **insert** (cargo) and the **replacement gDNA**
(the recoded deletion-replacement arm).

## Core design capabilities

- **Recode, don't delete.** The genomic window between `S` and the ngRNA nick is
  re-templated with a **synonymously recoded** copy of itself (same protein) chosen
  to break the perfect microhomology run that otherwise drives MMEJ/SSA deletion of
  the insert. Length-preserving (`deletion_size == 0`). Among the minimum-edit
  microhomology-breaking options, each **changed** codon is set to the most frequent
  human synonymous codon that still keeps the run below threshold, so the forced
  edits stay translation-competent while unchanged positions remain byte-exact
  (human Codon Adaptation Index tiebreak, ADR 0009).
- **Ranked-pool-walk orchestrator.** A deterministic candidate pool
  (synonymous RTT variants × PBS lengths × screened nicks) is walked
  cheap→expensive; the first candidate passing every gate wins. Exhaustion is
  fail-loud (`NoValidPegRNA`) with a per-candidate diagnostic naming the gate that
  killed each candidate. The insert site `S` never moves.
- **PBS2 length selection.** PBS2 is genome-fixed at the ngRNA nick, so its **length
  is scanned per nick** (10–18 nt) to hit the Tm band, decoupled from PBS1. (This
  replaced a fixed-length gate that collapsed whole candidate pools.)
- **INV-0 correctness contract.** The product differs from the reference by exactly:
  the cargo at `S`, protein-silent recode-span edits, byte-exact re-templating of any
  non-coding span, and ≤2 silent PAM edits. Enforced with a real translation check.
- **Strand convention + TPRT self-check.** All sequences are stored 5′→3′;
  antiparallel relationships use `reverse_complement`. Every candidate is run through
  an in-silico template-jump reconstruction (`simulate_tprt`) that must rebuild the
  intended edited allele, catching any orientation error loudly.
- **Coordinate-safe diploid genome.** KOLF2.1J is materialised as a real diploid on
  **GRCh38 coordinates** by applying phased substitutions per haplotype, so the
  GRCh38-trained gates stay valid (no liftover).
- **Frame-safe T7 promoters.** `T7-G` / `rcT7-G` insertable elements: reverse
  orientation is the reverse complement (directionality preserved), padded with the
  minimum guanines to a length multiple of 3 so the promoter never frameshifts a
  downstream CDS while keeping the G-rich transcription start.

## The gate stack

Each gate is a `Protocol` with an in-repo fake and a real adapter (see
[Modes](#modes-fake-vs-real)). Thresholds are literature-derived defaults in
`gate_thresholds.py`, each with a citation, and are overridable by a wet-lab
calibration set.

| Gate | Tool / model | Checks | Default threshold |
|---|---|---|---|
| Microhomology | DnaChisel (synonymous recoder) | longest perfect repeat between arm and flank | max run **8 bp** (recoder minimises below) |
| PBS2 | Biopython `Tm_NN` | second-strand primer stability | Tm **37–42 °C**, GC 40–70 %, length scan 10–18 nt |
| Cryptic splice (primary) | **MaxEntScan** (Yeo & Burge) | edit creates/strengthens a donor/acceptor motif in the recode span + cassette + junctions | MaxEnt < **8** bits (context-free) |
| Splice (cross-check) | **SpliceAI** (advisory) | native-context splice disruption, when ≥10 kb genomic context is supplied | Δ ≤ **0.2**; skipped on short context |
| Structural importance | **AlphaFold-DB** + biotite | recoded residue is buried / confident / packed | pLDDT ≥ 70, RSA < 0.25, contacts > 18 |
| Full-pegRNA 2D structure | **RiNALMo** (GPU) → **RNAformer** (CPU) fallback | functional regions (spacer/PBS/RTT) mis-pair | base-pair prob boundary **0.5** |
| Docking | **Boltz-2** (AF3-class ipTM) | scaffold–nCas9 interface degraded by the extension | ipTM floor 0.6, max drop 0.2 vs apo |
| Off-target | **Cas-OFFinder** (+ genome) | ngRNA nick specificity | 0 off-targets beyond the intended site |

Gate order (cheap→expensive): recode/microhomology → PBS2 → INV-0 → TPRT →
cryptic-splice (MaxEntScan; SpliceAI advisory cross-check) → AlphaFold importance →
RiNALMo/RNAformer 2D → Boltz-2 docking. Each heavy gate is
guarded so a scorer failure skips a candidate (with a diagnostic) rather than
aborting the walk.

## Modes: fake vs real

- **`--mode fake`** — the full design logic runs with deterministic stand-in scorers.
  No GPU, no weights; produces a real report. Used for CI, development, and dry-runs.
- **`--mode real`** — the real adapters run (SpliceAI, AlphaFold-DB, RiNALMo/RNAformer,
  Boltz-2, Cas-OFFinder). Requires the HPC environment and staged weights; the
  structure gate falls back RiNALMo→RNAformer (both real models) when FlashAttention
  is unavailable.

Scorers are **injected** — the workflow imports no heavy model at module load, so the
unit suite is fast and deterministic.

## Output

Two files per run: `recoded_tjpe_report.txt` (human-readable) and
`recoded_tjpe_components.csv` (`index,component,role,part_of,sequence,length`, one row
per labelled component). The RTT is split into `insert` and `replacement_gDNA`, and
`part_of` distinguishes the tjPE template-jumping portion (PBS2·RTT·PBS1) from the
spacer/scaffold. See [`DESIGN_HISTORY.md`](DESIGN_HISTORY.md) §5 for the gate detail.

## Validation status

- Full unit suite (~670 tests) is green and heavy-model-free.
- On the HPC, the login-node/CPU-capable gates are validated live (SpliceAI,
  AlphaFold-importance, structure via RNAformer); the off-target gate runs on CPU.
  The Boltz-2 and RiNALMo GPU paths are validated via the GPU smoke job.
- **Thresholds are literature-derived defaults, not yet wet-lab calibrated.** Treat
  outputs as engineering-grade until a calibration set is applied and the GPU smoke
  passes end-to-end.

## Roadmap (planned / not yet implemented)

Candidate additions, roughly ordered by biological impact:

1. **Intron/exon-boundary protection** — reject edits that fall on an intron–exon
   junction. Needs an exon-junction annotation source wired into the genome context
   (the current CDS mask distinguishes coding vs non-coding but not intron vs UTR).
2. **Strand-targetability validation** — confirm whether a locus is targetable on the
   intended strand (and flag when only the template strand is usable).
3. **Non-coding (UTR) microhomology handling** — the current design re-templates
   non-coding spans byte-exact (a locked decision); an alternative that perturbs a few
   non-functional bases to break microhomology in UTRs is a biological-judgment change,
   not a mechanical one, and should be gated behind explicit configuration.
4. **Harder comparison gates** — additional cross-checks to catch anomalous constructs
   beyond the current per-gate filters.
5. **Wet-lab threshold calibration** — replace the literature-derived defaults with
   values fit to an experimental calibration set.

## Testing

```bash
pip install -e '.[dev]'
ruff check src tests && black --check --line-length 100 src tests
pytest -q            # ~670 tests, heavy-model-free
```

Coverage spans the design logic (enumeration, ranked-pool-walk, INV-0, strand/TPRT,
arm partition, nick search), each gate's pure logic (PBS2 selection, splice
indel-delta, structure region-pairing, docking thresholds, off-target aggregation),
the scorer adapters (fakes + fail-loud paths, incl. the Boltz offline YAML), the
output decomposition + CSV, and an end-to-end example run. Real models are validated
separately on the HPC via `scripts/hpc/gpu_smoke.pbs`.
