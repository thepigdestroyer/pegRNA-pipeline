"""Tests for `pegrna motif list` and `pegrna motif show`."""

from click.testing import CliRunner

from pegrna_pipeline.agent.cli import cli


def test_motif_list_prints_names_grouped_by_category():
    runner = CliRunner()
    result = runner.invoke(cli, ["motif", "list"])
    assert result.exit_code == 0, result.output
    # A few representative names should appear.
    assert "FLAG" in result.output
    assert "GGGGS" in result.output
    assert "EGFP" in result.output


def test_motif_list_filters_by_category():
    runner = CliRunner()
    result = runner.invoke(cli, ["motif", "list", "--category", "linker"])
    assert result.exit_code == 0, result.output
    assert "GGGGS" in result.output
    # A non-linker motif should not appear.
    assert "EGFP" not in result.output


def test_motif_show_prints_details():
    runner = CliRunner()
    result = runner.invoke(cli, ["motif", "show", "FLAG"])
    assert result.exit_code == 0, result.output
    assert "FLAG" in result.output
    assert "DYKDDDDK" in result.output
    assert "GACTACAAGGACGACGATGACAAG" in result.output


def test_motif_show_unknown_exits_nonzero():
    runner = CliRunner()
    result = runner.invoke(cli, ["motif", "show", "TotallyFakeTag"])
    assert result.exit_code != 0
    assert "TotallyFakeTag" in result.output or "not found" in result.output.lower()


def test_motif_help_mentions_add_deferred():
    runner = CliRunner()
    result = runner.invoke(cli, ["motif", "--help"])
    assert result.exit_code == 0, result.output
    # The help docstring should hint that `motif add` lives in Plan 3.
    assert "Plan 3" in result.output or "add" in result.output.lower()
