"""CLI scaffold tests: --help and --version must work."""

from click.testing import CliRunner

from pegrna_pipeline import __version__
from pegrna_pipeline.agent.cli import cli


def test_cli_help_exits_zero():
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "pegRNA" in result.output or "pegrna" in result.output


def test_cli_version_prints_package_version():
    runner = CliRunner()
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output
