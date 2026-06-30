# 0005. Gate stack, external tools, and literature-derived thresholds

Date: 2026-06-18
Status: Accepted

## Context

The recoded TJ-PE design must be defended against several independent failure modes
(microhomology deletion, mis-splicing, structurally critical protein residues, pegRNA
mis-folding, loss of nCas9 docking). Each needs a gate, a tool, and a threshold.
Earlier placeholders were "wet-lab-pending," i.e. arbitrary. Thresholds had to become
defensible, and the heavy ML tools had to be chosen deliberately and isolated from the
fast unit suite.

## Decision

Gates run cheap→expensive in `run_recoded_tjpe`; the first candidate passing all wins;
exhaustion raises `NoValidPegRNA` with diagnostics. Heavy scorers are injected (fakes
in CI, real adapters on CX3 via `scorer_factory.build_scorers`), each guarded so a
scorer crash skips a candidate rather than aborting the walk.

Tools: DnaChisel (synonymous recode); SpliceAI, Pangolin alternative (splice);
AlphaFold-DB + biotite (residue importance); RiNALMo (full-pegRNA 2D structure) with
RNAformer as CPU fallback and ViennaRNA as cheap pre-screen; Boltz-2 (scaffold–nCas9
docking).

Thresholds (`gate_thresholds.py`, each with a citation): splice Δ ≤ 0.2 (Jaganathan
2019); pLDDT ≥ 70 (Jumper 2021); RSA buried < 0.25 (Levy 2010; Tien 2013); > 18
heavy-atom contacts within 5 Å = packed core (Yeh 2014; Shih 2012); docking ipTM floor
0.6, max drop 0.2 vs apo (Evans 2022; Abramson 2024); microhomology max run 8 bp (Shen
& Huang 1986; Sfeir & Symington 2015); RiNALMo pair-probability boundary 0.5. Nick
distance (40-70 bp, smaller preferred) is intentionally excluded as an experimental-
design knob, not a threshold. Thresholds are defaults overridable by a wet-lab
calibration library.

## Alternatives considered

| Alternative | Why rejected |
|---|---|
| ViennaRNA as the primary 2D structure gate | Thermodynamic MFE model, weak on long multi-domain pegRNAs; kept only as cheap pre-screen. |
| NuFold / RibonanzaNet for 2D structure | 3D / chemical-probing oriented, heavier, less suited to a pairwise base-pair-probability read than RiNALMo. |
| Hand-rolled codon recoder | DnaChisel already solves constrained synonymous recoding. |
| Keep arbitrary "wet-lab-pending" thresholds | Not defensible; replaced with cited literature defaults. |
| Import heavy scorers in the workflow module | Would make unit tests slow/GPU-bound; scorers are injected instead. |

## Consequences

- Positive: thresholds are defensible and individually citable; tool choices are
  recorded so they are not re-litigated.
- Positive: unit suite stays fast and deterministic; GPU deps isolated to CX3.
- Negative: RiNALMo and Boltz-2 are not yet validated live (CX3 smoke 2/5); off-cluster
  they run as fakes. Thresholds remain uncalibrated against wet-lab outcomes.

## References

- `src/pegrna_pipeline/hpc/gate_thresholds.py` (citations inline).
- `splice_screen.py`, `structure_importance.py`, `full_structure.py`, `complex_dock.py`,
  `recode_arm.py`, `microhomology.py`, `scorer_factory.py`.
- `docs/DESIGN_HISTORY.md` §5, §8.
