"""A FAILED design must still say which gates were real.

"No candidate survived all gates" means something rejected every design. Which
gates were real models and which were placeholders is precisely what the reader
needs at that point: a rejection by a stand-in is not evidence of anything, and
a rejection by a real model is. Previously gates.json was only written on
success, so a failed run disclosed nothing at all.
"""

import json
from pathlib import Path

from click.testing import CliRunner

from pegrna_pipeline.agent.cli import cli


def _broken_config(tmp_path: Path) -> Path:
    """A config that reaches the pipeline and then fails (no spacer, no runner)."""
    cfg = {
        "gene": "TEST",
        "uniprot_acc": "P49711",
        "rtt_aa": "L",
        "cargo": "GACTACAAGGACGACGATGACAAG",
        # Too short for any real design, so the run fails loud.
        "genome": "ACGT" * 20,
        "site_s": 40,
        "cds_mask": [True] * 40 + [False] * 40,
        "pbs_lengths": [13],
        "ngrna_distance_min": 40,
        "ngrna_distance_max": 90,
        "pegrna_strand": "Rv",
        "arm_side": "5prime",
        "max_run": 8,
    }
    p = tmp_path / "config.json"
    p.write_text(json.dumps(cfg))
    return p


def test_failed_fake_design_still_writes_stand_in_gates(tmp_path):
    out = tmp_path / "out"
    result = CliRunner().invoke(
        cli,
        [
            "design-recoded",
            "--config",
            str(_broken_config(tmp_path)),
            "--mode",
            "fake",
            "--output-dir",
            str(out),
        ],
    )
    assert result.exit_code == 1, result.output

    gates_path = out / "recoded_tjpe_gates.json"
    assert gates_path.exists(), "a failed run must still record what it ran"

    gates = json.loads(gates_path.read_text())
    assert gates["failed"] is True
    assert gates["scorer_mode"] == "fake"
    assert gates["error"]
    # fake mode is entirely stand-ins, so the rejection proves nothing.
    assert "docking" in gates["stand_in_gates"]


def test_the_diagnostic_never_masks_the_real_failure(tmp_path, monkeypatch):
    """If writing the diagnostic blows up, the original error still surfaces."""
    import pegrna_pipeline.hpc.scorer_factory as sf

    def _explode(**_kwargs):
        raise RuntimeError("scorer probe failed")

    monkeypatch.setattr(sf, "make_scorers", _explode)

    out = tmp_path / "out"
    result = CliRunner().invoke(
        cli,
        [
            "design-recoded",
            "--config",
            str(_broken_config(tmp_path)),
            "--mode",
            "fake",
            "--output-dir",
            str(out),
        ],
    )
    assert result.exit_code == 1
    assert "design-recoded failed" in result.output
    assert "scorer probe failed" not in result.output
