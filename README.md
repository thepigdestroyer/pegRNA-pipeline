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
