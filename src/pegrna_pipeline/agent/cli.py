"""Click-based CLI entry point (scaffold; real commands wired in Task 2+)."""

import click

from pegrna_pipeline import __version__


@click.group()
@click.version_option(__version__, prog_name="pegrna")
def cli() -> None:
    """pegRNA design pipeline agent."""


if __name__ == "__main__":  # pragma: no cover
    cli()
