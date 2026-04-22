"""Tests for `pegrna validate <construct>` CLI command."""

from click.testing import CliRunner

from pegrna_pipeline.agent.cli import cli


def test_validate_known_motif_exits_zero():
    runner = CliRunner()
    result = runner.invoke(cli, ["validate", "FLAG"])
    assert result.exit_code == 0, result.output
    # Output should mention the motif name and its aa/dna.
    assert "FLAG" in result.output
    assert "DYKDDDDK" in result.output


def test_validate_chained_construct():
    runner = CliRunner()
    result = runner.invoke(cli, ["validate", "GGGGS-FLAG"])
    assert result.exit_code == 0, result.output
    assert "GGGGS" in result.output
    assert "FLAG" in result.output


def test_validate_repeated_motif():
    runner = CliRunner()
    result = runner.invoke(cli, ["validate", "GGGGSx3"])
    assert result.exit_code == 0, result.output
    # Should repeat GGGGS three times in the resolved aa.
    assert result.output.count("GGGGS") >= 3


def test_validate_raw_aa():
    runner = CliRunner()
    result = runner.invoke(cli, ["validate", "aa:MK"])
    assert result.exit_code == 0, result.output
    # DNA should be 6 chars (MK codon-optimised).
    assert "MK" in result.output


def test_validate_unknown_motif_non_interactive_exits_nonzero():
    runner = CliRunner()
    result = runner.invoke(cli, ["validate", "TotallyFakeTag"])
    # Without --interactive, unknowns are a hard error.
    assert result.exit_code != 0
    assert "TotallyFakeTag" in result.output or "unknown" in result.output.lower()
