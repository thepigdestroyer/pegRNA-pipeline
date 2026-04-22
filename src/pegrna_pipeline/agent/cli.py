"""Click-based CLI entry point."""

from __future__ import annotations

import sys

import click

from pegrna_pipeline import __version__
from pegrna_pipeline.agent.construct_parser import (
    parse_components,
    resolve_construct,
    tokenize,
)
from pegrna_pipeline.agent.interactive_resolver import (
    AbortedByUser,
    build_interactive_resolver,
)
from pegrna_pipeline.agent.motif_library import (
    MotifLibrary,
    MotifNotFoundError,
)


@click.group()
@click.version_option(__version__, prog_name="pegrna")
def cli() -> None:
    """pegRNA design pipeline agent."""


@cli.command("validate")
@click.argument("construct")
@click.option(
    "--interactive/--no-interactive",
    default=False,
    help="Prompt for unknown motif names instead of erroring out.",
)
def validate_cmd(construct: str, interactive: bool) -> None:
    """Parse CONSTRUCT and print its resolved components, DNA, and aa."""
    library = MotifLibrary.load_full()
    raw_components = parse_components(tokenize(construct))

    if interactive:
        resolver = build_interactive_resolver(library)
        resolved = []
        try:
            for rc in raw_components:
                if rc.kind.name == "MOTIF" and rc.value not in library:
                    rc = resolver(rc.value)
                resolved.append(rc)
        except AbortedByUser as exc:
            click.echo(f"Aborted: {exc}", err=True)
            sys.exit(2)
        raw_components = resolved

    try:
        construct_obj = resolve_construct(raw_components, library=library)
    except MotifNotFoundError as exc:
        click.echo(f"unknown motif: {exc.args[0]}", err=True)
        sys.exit(1)

    click.echo(f"Construct: {construct}")
    click.echo(f"Components ({len(construct_obj.components)}):")
    for i, comp in enumerate(construct_obj.components, start=1):
        click.echo(
            f"  {i:>3}. {comp.source:<7} {comp.name:<20} "
            f"aa={comp.aa or '-':<12} dna={comp.dna}"
        )
    click.echo(f"Total aa ({len(construct_obj.aa)}): {construct_obj.aa}")
    click.echo(f"Total dna ({len(construct_obj.dna)}): {construct_obj.dna}")


@cli.group("motif")
def motif_group() -> None:
    """Inspect the motif library.

    The ``motif add`` subcommand for persistent user additions is deferred
    to Plan 3; today's workflow for custom motifs is to hand-edit
    ``~/.pegrna/motifs.yaml``.
    """


@motif_group.command("list")
@click.option("--category", default=None, help="Only show motifs in this category.")
def motif_list(category: str | None) -> None:
    """List motif names grouped by category."""
    library = MotifLibrary.load_full()
    rows: dict[str, list] = {}
    for m in library:
        if category is not None and m.category != category:
            continue
        rows.setdefault(m.category, []).append(m)
    for cat in sorted(rows):
        click.echo(f"[{cat}]")
        for m in sorted(rows[cat], key=lambda x: x.name.lower()):
            click.echo(f"  {m.name:<14} {m.aa:<32} {m.description}")


@motif_group.command("show")
@click.argument("name")
def motif_show(name: str) -> None:
    """Show details of a single motif by name (case-insensitive)."""
    library = MotifLibrary.load_full()
    try:
        m = library.lookup(name)
    except MotifNotFoundError:
        click.echo(f"Motif {name!r} not found.", err=True)
        sys.exit(1)
    click.echo(f"Name:        {m.name}")
    click.echo(f"Category:    {m.category}")
    click.echo(f"Description: {m.description}")
    click.echo(f"AA  ({len(m.aa):>3} aa):  {m.aa}")
    click.echo(f"DNA ({len(m.dna):>3} nt): {m.dna}")


if __name__ == "__main__":  # pragma: no cover
    cli()
