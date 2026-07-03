"""The shipped CTCF example config must stay valid and produce a report end-to-end.

Runs examples/ctcf_cterm/config.json through the real design driver in fake mode
(no GPU/weights). This is the regression guard for the reproducible CTCF example
the README points users at, and it also pins the PBS2-fix behaviour (the same config
previously failed with no_valid_pegrna under the old fixed-length PBS2 gate).
"""

from __future__ import annotations

import json
from pathlib import Path

from pegrna_pipeline.hpc.recoded_driver import design_recoded

_CONFIG = Path(__file__).resolve().parents[1] / "examples" / "ctcf_cterm" / "config.json"


def _load_design_params():
    cfg = json.loads(_CONFIG.read_text())
    # Drop provenance/metadata keys (the CLI routes these to the manifest).
    import inspect

    valid = set(inspect.signature(design_recoded).parameters)
    return {k: v for k, v in cfg.items() if not k.startswith("_") and k in valid}


def test_ctcf_example_config_exists_and_is_wellformed():
    cfg = json.loads(_CONFIG.read_text())
    assert cfg["gene"] == "CTCF"
    assert len(cfg["genome"]) == len(cfg["cds_mask"])  # mask aligns with the window
    assert 0 < cfg["site_s"] < len(cfg["genome"])
    assert cfg["arm_side"] == "5prime"  # C-terminal recode-into-CDS


def test_ctcf_example_produces_a_report(tmp_path):
    params = _load_design_params()
    out = design_recoded(mode="fake", output_dir=str(tmp_path), **params)
    # A pegRNA was designed (no NoValidPegRNA) and a report was written.
    assert out.result.pegRNA_sequence
    assert out.result.site_s == 147
    assert out.result.deletion_size == 0  # recode, not delete
    assert out.report_text and "Spacer.Scaffold.rcPBS2.RTT.PBS1.PolyU" in out.report_text
    assert Path(out.report_path).exists()
    # PBS2 length is chosen by the per-nick scan (the fix); it need not be 13/15/17.
    assert 10 <= len(out.result.rcpbs2) <= 18
