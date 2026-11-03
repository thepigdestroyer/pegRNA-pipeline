# Calibration analyses

Scripts that justify (or challenge) the gate thresholds. Every threshold in
`src/pegrna_pipeline/hpc/gate_thresholds.py` should be defensible with numbers from
here, not by assertion.

| Script | Question it answers |
|---|---|
| `offtarget_table.py` | Why did the nick screen reject (or keep) each candidate ngRNA? |
| `scan_structure_gate.py` | Which region pairs does the folding gate actually trip on, and at what probability? |
| `calibrate_structure_gate.py` | What fraction of *random* pegRNAs pass the folding gate at each threshold? |
| `calibrate_docking_gate.py` | **Measures** Boltz-2 ipTM for neutral, non-binder and scaffold-sequestering RNAs, and appends to `data/docking_calibration.jsonl` (spends NVIDIA credits; resumable) |
| `docking_reference.py` | **Derives** `DOCKING_REFERENCE` from that file, and tests the gate on it (leave-one-out, a held-out spacer, the sequestering sweep) |

A gate that rejects ~all random designs is not screening out bad ones; a threshold that sits
inside the model's own sampling noise is not a threshold. These scripts exist to tell the
two apart.

## The docking data

`data/docking_calibration.jsonl` is one JSON record per RNA: its class, spacer, extension length,
the sequence (where it was recorded), the 5 per-sample interface ipTMs, and how many scaffold
nucleotides ViennaRNA pairs with the extension. It is the evidence behind the docking reference
table, and `tests/test_docking_reference.py` fails if the table in `gate_thresholds.py` stops
matching it.

| source | what it is |
|---|---|
| `legacy-calibration`, `legacy-holdout` | the first measurements (2026-11-03), imported from files that lived in `/tmp`; sequences were not recorded except for the pipeline's own designs |
| `calibration` | the extended grid from `calibrate_docking_gate.py`, with sequences |
| `spacer2` | a second spacer, **never** used to build the table: it tests it |

```bash
python scripts/analysis/calibrate_docking_gate.py --dry-run          # what would be measured
python scripts/analysis/calibrate_docking_gate.py --per-length 6 --controls 4 --sweep
python scripts/analysis/docking_reference.py                         # tables and validation
python scripts/analysis/docking_reference.py --emit                  # the constant to paste
python scripts/analysis/docking_reference.py --check                 # exit 1 if gate_thresholds is stale
```

A run skips any RNA already in the file, so re-running after a stall or a rate limit continues
where it stopped. Use `--seed` to draw different random sequences.
