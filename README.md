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
