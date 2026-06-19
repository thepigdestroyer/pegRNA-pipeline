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

> **STATUS — design-logic layer is implemented and unit-tested; real model
> execution is NOT yet wired and runs only on GPU/HPC.** Each heavy gate
> (SpliceAI, AlphaFold, RiNALMo, Boltz-2) ships a Protocol + an in-repo fake +
> a real adapter whose inference body is staged for Imperial CX3 (this dev box
> lacks an Ampere+ GPU). Stage-0 genome ingest (real T2T-KOLF diploid fetch +
> KOLF↔GRCh38 chain) is fixture-backed here. The microhomology threshold L,
> RiNALMo/Boltz-2 cutoffs, and the sub-70 bp nick distance are **calibration
> parameters requiring wet-lab validation**. See
> `docs/superpowers/plans/2026-06-18-recoded-tjpe-pegrna-redesign.md`.
> **Not production-ready: do not deploy to end users until the heavy-model
> follow-on plan is executed and validated on CX3.**
