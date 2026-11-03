# Calibration analyses

Scripts that justify (or challenge) the gate thresholds. Every threshold in
`src/pegrna_pipeline/hpc/gate_thresholds.py` should be defensible with numbers from
here, not by assertion.

| Script | Question it answers |
|---|---|
| `offtarget_table.py` | Why did the nick screen reject (or keep) each candidate ngRNA? |
| `scan_structure_gate.py` | Which region pairs does the folding gate actually trip on, and at what probability? |
| `calibrate_structure_gate.py` | What fraction of *random* pegRNAs pass the folding gate at each threshold? |
| `calibrate_docking_gate.py` | What does Boltz-2 ipTM report for known controls, and how noisy is it? (spends NVIDIA credits) |

A gate that rejects ~all random designs is not screening out bad ones; a threshold that sits
inside the model's own sampling noise is not a threshold. These scripts exist to tell the
two apart.
