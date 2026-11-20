# pegrna_pipeline
pegrna.ukdri.org

Automated pegRNA design pipeline for prime editing: C-terminal protein tagging with
template-jumping PE (tj-PE), and cassettes or whole genes with Prime Assembly. A live
web app is at <https://pegrna.ukdri.org>.

> **New here? Read [`docs/DESIGN_HISTORY.md`](docs/DESIGN_HISTORY.md) first.** It is
> the canonical, self-contained explanation of the recoded TJ-PE designer: what it
> does, why each gate exists, how it is implemented, the full locked-decision /
> rationale history, the adversarial-review fixes, the HPC deployment, and the
> current status and limitations. The ADRs under [`docs/decisions/`](docs/decisions/)
> are the primary sources it synthesises.

## Modes

The `pegrna run --mode <mode>` command (the cluster path) supports two pipeline modes, picked
according to insert size and the editing strategy needed:

| Mode   | Insert size       | Strategy                                                                    |
|--------|-------------------|-----------------------------------------------------------------------------|
| `pe7`  | ≤ 40 bp           | Standard PE7 (PE3-strategy): one pegRNA + one opposite-strand ngRNA pair.   |
| `tjpe` | ≤ ~500 bp         | Template-jumping PE: codon-enumerated first 8 aa + ngRNA-driven RC-PBS2.    |

Each mode dispatches in `pegrna_pipeline.hpc.driver.run_driver` and writes a
text + CSV report via `pegrna_pipeline.hpc.report_writer.write_reports`.

**Prime Assembly is `pegrna design-pa`, not a `run` mode.** For a cassette or a gene the insert is
better carried by a separate donor DNA than by a pegRNA that grows with it. Prime Assembly
(Liu et al., *Nature* 2026) does that: two pegRNAs on opposite strands write ~35 nt flaps at two
nicks, and the ends of a linear donor match them. Here the donor also restores the DNA between the
nicks, so nothing is deleted. `design-pa` chooses the guides itself, runs the same gates (folding,
docking, splicing, off-target) with real models, and writes both pegRNAs, the donor and a list of
every nucleotide that differs from your input. It needs an insert of at least two flaps (70 nt at the
default flap of 35), and it is the web page's second way to write the insert ("As a separate donor DNA"). **No donor that restores the
DNA between the nicks has been tested at the bench**, and every result says so. The design is
specified in `docs/superpowers/specs/2026-11-04-prime-assembly-design.md`, drawn in
[docs/STRUCTURE.md](docs/STRUCTURE.md) §2, and its real-model evidence and open points are in
[ADR 0030](docs/decisions/0030-prime-assembly-design-pa.md). The earlier prototype (`run --mode pa`)
could never pair real PRIDICT output ([ADR 0027](docs/decisions/0027-prime-assembly-status-and-what-a-designer-needs.md))
and was removed; `pegrna run --mode pa` now refuses, with a pointer to `design-pa`.

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
Spacer · Scaffold · rcPBS2 · RTT(Insert + Deletion-Replacement-Arm) · PBS1 · PolyU(2-3 U)
```

(rcPBS2, the reverse complement of PBS2, sits 5' of the RTT; PBS2 itself is synthesised onto the new
flap and captures the strand the ngRNA nicked. PolyU is the La/PE7 3' tail.
[docs/STRUCTURE.md](docs/STRUCTURE.md) draws the whole structure, with numbers from a real run.)

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

### Insertable cargo elements: T7-G / rcT7-G (frame-safe T7 promoters)

`insert_elements.assemble_insert` builds the cargo from named segments, including
two directional, frame-safe T7 promoters for linear RNA amplification:

- `T7-G` (`TAATACGACTCACTATAGGGG`, 21 nt) transcribes the segments **downstream** (3′).
- `rcT7-G` (`CCCCTATAGTGAGTCGTATTA`, 21 nt) transcribes the segments **upstream** (5′).

Two correctness points are baked in:
- **Reverse orientation is the reverse *complement* of the promoter**, not the
  reverse. Inserting `rcT7-G` puts the functional promoter 5′→3′ on the *non-targeted*
  strand, preserving its directionality. (Reversing the letters would place it 3′→5′
  and break it.)
- **Frame-safe with a guanine start.** T7 RNA polymerase initiates at a **G** and the
  start is G-rich (consensus `TAATACGACTCACTATAGGGAGA`, transcription begins at the
  first G; Rong et al. 1998 *PNAS* 95:515). The element is the coding strand, so the
  initiator is G there. The 20-nt core isn't a multiple of 3, so mid-cargo it would
  frameshift the downstream CDS (e.g. ZeoR); we append the minimum **G** bases (one
  for the 20-mer → 21 nt) to fix the frame while keeping the start G-rich. (C would be
  wrong — cytosines on the coding strand are not the initiator.)

e.g. `FLAG-GGGGS-T7-G-P2A-ZeoR` → `T7-G` transcribes `P2A-ZeoR`; swap in `rcT7-G` →
transcribes `GGGGS-FLAG`.

### Installation & running on the HPC (one command)

Clone the repo on an **Imperial CX3 login node** (it needs internet to build the env
and download weights), then run the installer:

```bash
git clone https://github.com/thepigdestroyer/pegRNA-pipeline && cd pegRNA-pipeline
export PEGRNA_DATA=$EPHEMERAL/pegrna-data          # where weights/genome are staged
bash scripts/hpc/install.sh                         # env + package + weights + login-node smoke
```

`install.sh` is idempotent. It (1) creates the `pegrna-recoded` conda env from
`scripts/hpc/environment.yml`, (2) `pip install -e '.[recode,splice,importance,fullstructure,complex]'`,
(3) stages all model weights + the GRCh38 reference (`stage_models.sh`), and (4) runs
the login-node gate smoke (SpliceAI, AlphaFold-importance, structure via RNAformer).

Then submit the GPU jobs:

```bash
# All 5 gate models on a GPU node (target 5/5):
qsub -v PEGRNA_DATA=$PEGRNA_DATA scripts/hpc/gpu_smoke.pbs

# End-to-end recoded TJ-PE design with real gates (ships a real CTCF example):
qsub -v PEGRNA_DATA=$PEGRNA_DATA,CONFIG=examples/ctcf_cterm/config.json \
     scripts/hpc/run_recoded_pipeline.pbs
```

### Running a design (any machine, no GPU)

```bash
# Dry-run: full design logic with stand-in scorers (produces a real report):
pegrna design-recoded --config examples/ctcf_cterm/config.json --mode fake --output-dir results/dry

# Prime Assembly: two pegRNAs and a donor for a cassette (the example is the 499 nt CTCF cassette):
pegrna design-pa --config examples/ctcf_cterm/config_pa.json --mode fake --output-dir results/dry_pa
```

The config is a JSON of the design inputs (`rtt_aa`, `cargo`, `genome`,
`site_s`, `cds_mask`, `arm_side`, ...); keys starting with `_` or provenance keys
(`gene`, `transcript`, `uniprot_acc`) are routed into the report manifest. The
`spacer` is **optional** — omit it and PRIDICT2 auto-designs the pegRNA#1 spacer
from `site_s`; include `spacer` only to pin a fixed guide.

**What may be recoded.** Only bases that code for protein are ever changed in the cargo; a promoter (T7, either
orientation) is never touched. Say more with two optional keys: `cargo_parts`, a list of
`{"name", "dna", "locked"}` whose DNA joins to `cargo` (an unmarked part is locked; `pegrna validate "FLAG-GGGGS-..."`
prints a paste-ready list), and `cargo_locked`, `[first, last]` nucleotide ranges (1-based, inclusive) that must
never change. Without them only the ORF up to the first stop codon counts. The report lists every nucleotide that
differs from what you gave (`recoded_tjpe_changes.csv`; ADR 0029). See
`examples/ctcf_cterm/config.json` for a complete, runnable CTCF C-terminal knock-in
(explicit spacer) and `examples/ctcf_cterm/config_autospacer.json` for the same
design with the spacer auto-designed.

### How it runs (no AI agent required)

The pipeline is a plain command-line tool — you do **not** need an AI/LLM agent to run
it. (The Python package is named `...agent...`, but that layer is just an SSH job
launcher, not an AI.) Three ways to run, in increasing automation:

1. **Directly (recommended for most use).** Run the CLI yourself:
   - dry run anywhere: `pegrna design-recoded --config <cfg> --mode fake --output-dir <out>`
   - full run on the HPC as a batch job: `qsub -v PEGRNA_DATA=$PEGRNA_DATA,CONFIG=<cfg> scripts/hpc/run_recoded_pipeline.pbs`
2. **Via the built-in launcher (optional automation).** If you'd rather not SSH to the
   HPC by hand, the `pegrna` CLI can orchestrate a run from your laptop — it validates
   the construct locally, submits the HPC job, polls for completion, and downloads the
   reports:
   ```bash
   pegrna setup                      # once: write ~/.pegrna/config.yaml (HPC host, user)
   pegrna run --gene CTCF --construct "<construct>" --mode tjpe --output-dir results
   pegrna status        # list runs        pegrna tail <run-id>   # stream the driver log
   pegrna list          # recorded runs    pegrna cancel <run-id> # stop a run
   ```
   The launcher currently drives the `tjpe` / `pe7` workflows; the recoded design
   is run via the batch script or the `design-recoded` CLI above, and Prime Assembly via `design-pa`.
3. **From your own agent/automation.** Because everything is a normal CLI with JSON
   config in and text/CSV out, any external agent or workflow system can call it — it
   just shells out to `pegrna design-recoded` (or `qsub` the PBS script) and reads the
   files from `--output-dir`. No special integration is needed.

The only hard requirement for a *full* (real-gate) run is the **HPC environment**
(GPU + staged model weights), not an agent; `--mode fake` runs anywhere.

### Output

Each `design-recoded` (tj-PE) run writes these files to `--output-dir`:

- `recoded_tjpe_report.txt` — human-readable design + coordinates + gate results.
- `recoded_tjpe_components.csv` — one row per labelled construct component
  (`index,component,role,part_of,sequence,length`) for table/list consumption.
- `recoded_tjpe_changes.csv` — every nucleotide that differs from what you specified, with its
  part and codon (ADR 0029).
- `recoded_tjpe_gates.json` — the result of every gate, for scripts and records.
- `recoded_tjpe_verification.docx` — a colour-annotated manual cross-check (ADR 0012), written
  only when the optional `python-docx` package is installed.

A `design-pa` (Prime Assembly) run writes five: `recoded_pa_report.txt`, `recoded_pa_components.csv`,
`recoded_pa_changes.csv`, `recoded_pa_donor.fasta` (pegRNA A, pegRNA B and the donor, the three
sequences to order) and `recoded_pa_gates.json`.

The tj-PE construct is reported 5′→3′ as **Spacer · Scaffold · rcPBS2 · RTT · PBS1 · PolyU**
(rcPBS2 is the reverse complement of PBS2, ADR 0014).
The **RTT is decomposed** into its labelled sub-parts rather than a single opaque
block: `insert` (the cargo payload) and `replacement_gDNA` (the synonymously recoded
deletion-replacement arm) — these reconcatenate to the RTT exactly. `part_of`
distinguishes the **tjPE template-jumping portion** (rcPBS2 · RTT · PBS1) from the
spacer/scaffold that complete the full construct.

Full install/usage/troubleshooting: [`scripts/hpc/README.md`](scripts/hpc/README.md)
and [`docs/TROUBLESHOOTING.md`](docs/TROUBLESHOOTING.md).

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
| docking (interface ipTM) | **midway** between the length-matched *neutral* and *non-binder* references | Measured on Boltz-2, 5 samples averaged. ipTM falls with RNA length regardless of design quality, so the AF-Multimer floor of 0.6 is unreachable for tjPE-sized pegRNAs; a scaffold-sequestering extension scores 0.2–0.3 below neutral, so the gate still catches it. **Judged only up to 232 nt**: past that the model's scores for a normal RNA and a non-binder converge, so a longer design is reported against its own length and flagged, not accepted or rejected. [ADR 0024](docs/decisions/0024-calibrated-docking-and-folding-gates.md) |
| microhomology max-run | **8 bp** | RecA MEPS ~8 bp / MMEJ onset ~6–8 bp; recoder minimises below this |
| pegRNA folding, RiNALMo pair prob | **0.5** | sigmoid pair-call boundary (RiNALMo val-tunes 0.01–0.30 for F1) |
| pegRNA folding, ViennaRNA pair prob | **0.99** | calibrated against 300 random pegRNAs per size: at 0.5 only 9–12% of random tjPE-sized designs pass, at 0.99 80–84% ([ADR 0024](docs/decisions/0024-calibrated-docking-and-folding-gates.md)) |
| nick distance | **40–70 bp** | experimental design; smaller preferred (NOT auto-calibrated) |

These are *rational priors*; final tuning still benefits from a wet-lab
calibration set, but the pipeline no longer ships unjustified numbers.

Two of them started life as values borrowed from a different tool or a different
molecule size, and were replaced once measured (the docking floor and the ViennaRNA
folding cutoff). A threshold that rejects nearly every *random* design is not
screening out bad ones, and a threshold that sits inside the model's own sampling noise
is not a threshold. `scripts/analysis/` holds the scripts that produce the numbers.

## Testing

```bash
pip install -e '.[dev]'
ruff check src tests && black --check --line-length 100 src tests   # lint (CI gates on both)
pytest -q                                                            # ~660 tests
```

The suite (≈660 tests, 29 skipped) is heavy-model-free — every real gate is a
`Protocol` with an in-repo fake, so the whole design logic is exercised on CPU with
no weights. Coverage by area:

- **Design logic**: candidate enumeration, the ranked-pool-walk, INV-0, strand
  convention (`simulate_tprt`), arm partition/recode, nick search.
- **Gates (pure logic)**: PBS2 length-selection (`test_pbs2_selection.py`), splice
  indel-delta, structure region-pairing, docking thresholds, off-target aggregation.
- **Scorer adapters (fakes + fail-loud paths)**: `test_scorer_fixes.py`,
  `test_complex_dock.py` (incl. the Boltz offline single-sequence YAML).
- **End-to-end**: `test_ctcf_example.py` runs the shipped CTCF example config through
  the driver in fake mode and asserts a report is produced.
- **Real models**: validated separately on CX3 via `scripts/hpc/gpu_smoke.pbs`
  (`smoke_models.py`), not in the unit suite.

## Documentation map

| Doc | What |
|---|---|
| [`docs/USER_GUIDE.docx`](docs/USER_GUIDE.docx) | **Non-technical manual for biologists** — what it does, what you provide, how to run it, how to read the output. |
| [`docs/FEATURES.md`](docs/FEATURES.md) | Technical capability + roadmap reference (ML-biology audience): models, gate stack, algorithm, planned features. |
| [`docs/DESIGN_HISTORY.md`](docs/DESIGN_HISTORY.md) | **How it works + full decision/rationale history**: architecture, the gate stack, algorithm, genome model, every locked decision, adversarial-review fixes. |
| [`docs/decisions/`](docs/decisions/) | ADRs (0001-0007) — the primary sources DESIGN_HISTORY synthesises. |
| [`docs/TROUBLESHOOTING.md`](docs/TROUBLESHOOTING.md) | Install/run failure modes on CX3 and their fixes. |
| [`scripts/hpc/README.md`](scripts/hpc/README.md) | HPC runbook: install, staging, smoke, PBS submit. |
| [`examples/ctcf_cterm/`](examples/ctcf_cterm/) | A complete, runnable CTCF C-terminal knock-in design config. |

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

## Licence and citation

MIT licence, see [`LICENSE`](LICENSE). Third-party code and terms are listed in
[`NOTICE`](NOTICE). To cite pegRNA, use [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this
repository" button reads it) or the Zenodo record for the release.
