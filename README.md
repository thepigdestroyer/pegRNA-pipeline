# pegrna_pipeline

Automated pegRNA design pipeline for prime editing C-terminal protein tagging.

## Modes

The `pegrna run --mode <mode>` command supports three pipeline modes, picked
according to insert size and the editing strategy needed:

| Mode   | Insert size       | Strategy                                                                    |
|--------|-------------------|-----------------------------------------------------------------------------|
| `pe7`  | ≤ 40 bp           | Standard PE7 (PE3-strategy): one pegRNA + one opposite-strand ngRNA pair.   |
| `tjpe` | ≤ ~500 bp         | Tethered-junction PE: codon-enumerated first 8 aa + ngRNA-driven RC-PBS2.   |
| `pa`   | > 500 bp / >1 kb  | Prime Assembly: two pegRNAs install upstream + downstream flaps + gBlock.   |

Each mode dispatches in `pegrna_pipeline.hpc.driver.run_driver` and writes a
text + CSV report via `pegrna_pipeline.hpc.report_writer.write_reports`.

Foundations layer (agent + codon expander) in active development.

## Scarless tjPE mode (experimental)

`tjpe` normally **deletes** the genomic span between the pegRNA nick and the
nicking-sgRNA nick. The scarless mode (`run_tjpe_workflow(..., scarless=True)`)
re-templates that span so the product is the reference allele plus the cargo at
site S, with **zero endogenous deletion** (enforced by the `INV-0` hard gate;
at most one silent nicking-guide-PAM base may change). Reverse-strand pegRNAs
require the nicking nick exactly at S and template the cargo as its reverse
complement.

Install the optional structure dependency:

```bash
pip install -e '.[scarless]'      # or: conda install -c bioconda viennarna
```

> **WARNING (experimental, validate before use).** Scarless re-templating makes
> the restored segment 100% homologous to the genome the flap must displace —
> the exact microhomology condition Zheng 2023 (Fah) engineered *away* because
> it caused imperfect edits. The ~40 bp nick distance is also below the 70 bp
> minimum shown to work. These risks are surfaced as advisory warnings on the
> result, not blocked. See the risk register in
> `docs/superpowers/specs/2026-06-09-scarless-tjpe-design.md` (vault note 306 §10)
> and the implementation plan in
> `docs/superpowers/plans/2026-06-09-scarless-tjpe.md`. Wet-lab validation is
> required before relying on this mode.

## Recoded TJ-pegRNA pipeline (design-logic layer)

The recoded pipeline (`pegrna_pipeline.hpc.recoded_tjpe_workflow.run_recoded_tjpe`)
re-templates the deleted span with **synonymous codons** (CDS) to break
microhomology while preserving the protein, and runs a stacked gate set on a
**ranked-pool walk** (build the pool once; cheap gates prune early; expensive
models run last; the first candidate passing every gate wins; exhausting the
pool fails loud).

Authoritative pegRNA structure (5'→3'):

```
Spacer · Scaffold · PBS2 · RTT(Insert + Deletion-Replacement-Arm) · PBS1 · PolyU(2-3 U)
```

(PBS2 occupies the homology-arm slot; PolyU is the La/PE7 3' tail.)

**Gate stack:** INV-0 (≤2 silent PAM edits + non-coding `locked_span` byte-exact)
· splice (SpliceAI/Pangolin) · AlphaFold structural importance · RiNALMo 2D
self-pairing (RNAformer CPU fallback; ViennaRNA cheap pre-screen) · Boltz-2
scaffold–nCas9 docking. Genome source: the **T2T diploid KOLF2.1J assembly**
(both haplotypes), mapped to GRCh38 for the GRCh38-trained gates. An opt-in
`allele_hunt` mode may move the *nick* (never the insert site S) to capture a
discriminating het SNP; non-coding bases in the re-templated span are
re-templated byte-exact.

Optional dependency extras:

```bash
pip install -e '.[recode]'        # DnaChisel constrained recoder
pip install -e '.[splice]'        # SpliceAI + TensorFlow            (HPC/GPU)
pip install -e '.[importance]'    # biotite (+ a DSSP binary)        for AlphaFold metrics
pip install -e '.[fullstructure]' # torch (RiNALMo / RNAformer)      (RiNALMo: Ampere+ GPU)
pip install -e '.[complex]'       # Boltz-2                          (large GPU + MSA)
```

### Strand convention (locked)

Everything is stored **5′→3′**; antiparallel/template relationships are
`reverse_complement`. The manually-designed RTT tail in the emitted oligo is the
reverse complement of the desired sense flap; PBS2 equals the captured-strand
segment (complement-of-complement). `complement`, `reverse`, and
`reverse_complement` are three distinct operations — `tprt_sim.simulate_tprt`
reconstructs the integrated **sense** allele from the designed extension and
fails loudly on any wrong operation (its tests prove it rejects plain
`complement` and plain `reverse`).

### Insertable cargo elements: T7 / T7-rev

`insert_elements.assemble_insert` builds the cargo from named segments, including
two directional T7 promoters for linear RNA amplification:

- `T7` (`TAATACGACTCACTATAGGG`) transcribes the segments **downstream** (3′).
- `T7_rev` (reverse complement) transcribes the segments **upstream** (5′).

e.g. `FLAG-GGGGS-T7-P2A-ZeoR` → T7 transcribes `P2A-ZeoR`; swap in `T7_rev` →
transcribes `GGGGS-FLAG`.

### Running it

```bash
# Dry-run anywhere (full design logic, stand-in scorers, no GPU):
pegrna design-recoded --config design_config.json --mode fake --output-dir results/dry

# Real end-to-end with all gate models on Imperial CX3:
#   see scripts/hpc/README.md  (env, model staging, smoke validation, PBS submit)
```

### Gate thresholds (literature-derived defaults)

Defaults now come from published cutoffs (`hpc/gate_thresholds.py`, each with a
citation), not arbitrary placeholders. A wet-lab calibration library can still
override them, but they are rationally defensible out of the box:

| Gate | Default | Basis |
|---|---|---|
| splice Δ | **0.2** | SpliceAI high-recall cutoff (Jaganathan 2019 Cell); safety gate maximises detection |
| AlphaFold pLDDT | **≥70** | "confident" band (Jumper 2021; AFDB) — gate-enabler for burial calls |
| AlphaFold RSA buried | **<0.25** | interior/core rASA (Levy 2010); Tien 2013 max-ASA normalisation |
| AlphaFold contact no. | **>18** | densely-packed core; CN ~ mutational sensitivity (Yeh 2014) |
| docking ipTM floor | **0.6** | AF-Multimer/AF3 "likely incorrect" below 0.6 (Evans 2022; Abramson 2024) |
| docking ipTM drop | **0.2** | one full confidence band vs apo scaffold–Cas9 baseline |
| microhomology max-run | **8 bp** | RecA MEPS ~8 bp / MMEJ onset ~6–8 bp; recoder minimises below this |
| RiNALMo pair prob | **0.5** | sigmoid pair-call boundary (RiNALMo val-tunes 0.01–0.30 for F1) |
| nick distance | **40–70 bp** | experimental design; smaller preferred (NOT auto-calibrated) |

These are *rational priors*; final tuning still benefits from a wet-lab
calibration set, but the pipeline no longer ships unjustified numbers.

> **STATUS — design-logic layer is implemented and unit-tested; real model
> execution runs only on GPU/HPC.** Each heavy gate (SpliceAI, AlphaFold,
> RiNALMo, RNAformer, Boltz-2) ships a Protocol + in-repo fake + a real adapter
> written against the tool's documented API. Validated on Imperial CX3: the full
> suite runs, the end-to-end pipeline runs (fake mode), and the **ViennaRNA and
> AlphaFold-importance gates run live**; `scripts/hpc/setup_and_smoke.pbs` brings
> up the GPU gates (SpliceAI, RiNALMo, RNAformer, Boltz-2) and `smoke_models.py`
> must pass 5/5 before a production design. Genome ingest is real Ensembl-GRCh38
> + KOLF variant overlay (`hpc/genome_ingest.py`); the full T2T-KOLF diploid
> FASTA swaps in later behind the same interface. Gate thresholds are now
> literature-derived (above); the sub-70 bp nick distance and final threshold
> tuning still benefit from wet-lab calibration. See
> `docs/superpowers/plans/2026-06-18-recoded-tjpe-pegrna-redesign.md` and
> `scripts/hpc/README.md`.
> **Production-ready only after `smoke_models.py` passes 5/5 on CX3. Until then,
> outputs are engineering previews, not bench-ready designs.**
