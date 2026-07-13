"""Stage 3c: the ``design-recoded`` CLI routes to the het-edit dual pass.

Before this stage the CLI only ever called the single-pass ``design_recoded``;
the dual-pass ``design_recoded_with_het_edit`` (antisense primary + an
allele-specific sense-window pass) was unreachable from the command line and the
web. These tests pin the additive, backward-compatible routing entirely through
the CLI in ``--mode fake`` (deterministic stand-in scorers + PridictRunnerFake,
no PRIDICT/network):

* het requested (``het_edit`` truthy) AND a discriminating phased het present ->
  the dual pass runs and BOTH designs are written (primary at the top level for
  the web/downloads, het under ``het_edit/``), and ``het_edit.json`` records the
  produced pair + note.
* het requested but NO discriminating het -> the primary is still produced and
  ``het_edit.json`` carries the driver's explanatory note (no crash, exit 0).
* no het requested -> single pass, unchanged (regression guard: no ``het_edit``
  subdir, no ``het_edit.json``).
* auto-spacer flows to BOTH passes (PRIDICT designs each spacer) when no explicit
  ``spacer`` is supplied.
"""

from __future__ import annotations

import json

from click.testing import CliRunner

from pegrna_pipeline.agent.cli import cli

# The fake scorer bundle's PRIDICT stand-in emits this pegRNA#1 spacer for every
# row (see tests/test_autospacer_wiring.py), so the auto path resolves to it.
_FAKE_BUNDLE_SPACER = "GAATTCGTGTTCCTAAGCGT"

_GENOME = "ACGATGCATGGCATGACTGA" * 16  # 320 nt; the known-passing fake-mode fixture
_SPACER = "ACGTACGTACGTACGTACGT"

# A discriminating phased het (PAM-destroying -> strength 1.0) and a
# non-discriminating one (distal spacer position -> strength 0.0).
_DISCRIMINATING = {"element": "pegrna_pam", "offset": 1, "haplotype": 1, "phased": True}
_NON_DISCRIMINATING = {"element": "spacer", "offset": 0, "haplotype": 1, "phased": True}


def _base_config() -> dict:
    """Antisense primary config (mirrors the known-passing het-edit fixture)."""
    return {
        "spacer": _SPACER,
        "rtt_aa": "MKP",
        "cargo": "ATGAAACCCGGG",
        "genome": _GENOME,
        "site_s": 150,
        "cds_mask": [True] * 320,
        "pbs_lengths": [15],
        "pegrna_strand": "Rv",
        "max_run": 8,
    }


def _run(cfg: dict, tmp_path):
    cfg_path = tmp_path / "config.json"
    cfg_path.write_text(json.dumps(cfg))
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(
        cli,
        [
            "design-recoded",
            "--config",
            str(cfg_path),
            "--mode",
            "fake",
            "--output-dir",
            str(out_dir),
        ],
    )
    return result, out_dir


def test_dual_pass_runs_when_discriminating_het_present(tmp_path):
    cfg = _base_config()
    cfg["het_edit"] = True
    cfg["het_snps"] = [_DISCRIMINATING]

    result, out_dir = _run(cfg, tmp_path)

    assert result.exit_code == 0, result.output
    # Primary artifacts are replicated to the top level so the web + downloads
    # (which read run_dir/recoded_tjpe_*.{txt,csv,json}) keep working unchanged.
    assert (out_dir / "recoded_tjpe_report.txt").exists()
    assert (out_dir / "recoded_tjpe_components.csv").exists()
    assert (out_dir / "recoded_tjpe_gates.json").exists()
    # The dual pass wrote both sub-designs.
    assert (out_dir / "primary" / "recoded_tjpe_report.txt").exists()
    assert (out_dir / "het_edit" / "recoded_tjpe_report.txt").exists()
    # het_edit.json records the produced pair + a het-edit note.
    het = json.loads((out_dir / "het_edit.json").read_text())
    assert het["requested"] is True
    assert het["produced"] is True
    assert "het-edit" in het["note"].lower()
    assert "het-edit" in result.output.lower()


def test_primary_only_when_no_discriminating_het(tmp_path):
    cfg = _base_config()
    cfg["het_edit"] = True
    cfg["het_snps"] = [_NON_DISCRIMINATING]

    result, out_dir = _run(cfg, tmp_path)

    assert result.exit_code == 0, result.output  # advisory het failure never kills the run
    assert (out_dir / "recoded_tjpe_report.txt").exists()
    het = json.loads((out_dir / "het_edit.json").read_text())
    assert het["requested"] is True
    assert het["produced"] is False
    assert "no discriminating" in het["note"].lower()
    # No het pair was designed, so no het_edit sub-design exists.
    assert not (out_dir / "het_edit").exists()


def test_single_pass_unchanged_when_no_het_requested(tmp_path):
    cfg = _base_config()  # no het_edit, no het_snps

    result, out_dir = _run(cfg, tmp_path)

    assert result.exit_code == 0, result.output
    assert (out_dir / "recoded_tjpe_report.txt").exists()
    # Single pass writes top-level artifacts directly: no dual-pass side files.
    assert not (out_dir / "het_edit.json").exists()
    assert not (out_dir / "primary").exists()
    assert not (out_dir / "het_edit").exists()


def test_auto_spacer_flows_to_both_passes(tmp_path):
    cfg = _base_config()
    del cfg["spacer"]  # -> auto_spacer: PRIDICT designs each pass's pegRNA#1 spacer
    cfg["het_edit"] = True
    cfg["het_snps"] = [_DISCRIMINATING]

    result, out_dir = _run(cfg, tmp_path)

    # Exit 0 proves neither pass raised the "no spacer supplied" error, i.e. the
    # het pass also went through PRIDICT auto-spacer design.
    assert result.exit_code == 0, result.output
    het = json.loads((out_dir / "het_edit.json").read_text())
    assert het["produced"] is True
    # The fake PRIDICT bundle spacer is the one both passes adopted.
    assert _FAKE_BUNDLE_SPACER in (out_dir / "recoded_tjpe_components.csv").read_text()
