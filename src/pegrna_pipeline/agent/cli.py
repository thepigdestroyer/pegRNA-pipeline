"""Click-based CLI entry point."""

from __future__ import annotations

import getpass
import os
import shlex
import sys
from pathlib import Path
from typing import Callable

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


def _resolve_hpc_password_provider(non_interactive: bool) -> Callable[[], str]:
    """Build a password-provider callable for the HPC SSH client.

    Priority order:

    1. ``CX3_PASSWORD`` environment variable.
    2. ``HPC_PASS=...`` line in ``~/.cx3-credentials`` (chmod 600 file
       written by the ``hpc-connect`` skill / ``connect.sh`` setup).
    3. Interactive ``getpass.getpass`` prompt — disabled when
       ``--non-interactive`` is set, in which case launcher will get an
       empty password and SSH will fail loudly with a clear auth error.
    """
    pw = os.environ.get("CX3_PASSWORD")
    if pw:
        return lambda: pw

    cred_path = Path("~/.cx3-credentials").expanduser()
    if cred_path.exists():
        for line in cred_path.read_text().splitlines():
            if line.startswith("HPC_PASS="):
                rhs = line.split("=", 1)[1]
                # The credentials file is typically written by bash via
                # `printf 'HPC_PASS=%q\n' "$HPC_PASS"`, which shell-quotes
                # specials (e.g. ``a*b`` becomes ``a\*b``). shlex.split
                # round-trips those forms cleanly.
                tokens = shlex.split(rhs)
                value = tokens[0] if tokens else ""
                return lambda: value

    if non_interactive:
        return lambda: ""
    return lambda: getpass.getpass("HPC password: ")


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
            f"  {i:>3}. {comp.source:<7} {comp.name:<20} " f"aa={comp.aa or '-':<12} dna={comp.dna}"
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


@cli.command("run")
@click.option("--gene", required=True)
@click.option("--construct", required=True)
@click.option("--mode", type=click.Choice(["tjpe", "pa", "pe7"]), required=True)
@click.option("--cell-line", default="KOLF2_iPSC")
@click.option("--output-dir", default="results", type=click.Path())
@click.option("--cleanup", is_flag=True)
@click.option("--non-interactive", is_flag=True)
def run_cmd(
    gene: str,
    construct: str,
    mode: str,
    cell_line: str,
    output_dir: str,
    cleanup: bool,
    non_interactive: bool,
) -> None:
    """Submit a pegRNA design run to HPC."""
    from pathlib import Path

    import yaml

    from pegrna_pipeline.agent.db import RunsDB
    from pegrna_pipeline.agent.launcher import HPCLauncher

    config_path = Path("~/.pegrna/config.yaml").expanduser()
    if not config_path.exists():
        click.echo("Run `pegrna setup` first to create ~/.pegrna/config.yaml.", err=True)
        sys.exit(2)
    config = yaml.safe_load(config_path.read_text())
    config.setdefault("cell_line", cell_line)

    # Resolve locally for validation (errors surfaced before HPC contact)
    library = MotifLibrary.load_full()
    raw = parse_components(tokenize(construct))
    if non_interactive:
        resolved_components = raw
    else:
        resolver = build_interactive_resolver(library)
        resolved_components = []
        try:
            for rc in raw:
                if rc.kind.name == "MOTIF" and rc.value not in library:
                    rc = resolver(rc.value)
                resolved_components.append(rc)
        except AbortedByUser as exc:
            click.echo(f"Aborted: {exc}", err=True)
            sys.exit(2)
    try:
        resolve_construct(resolved_components, library=library)
    except MotifNotFoundError as exc:
        click.echo(f"Unresolved motif: {exc.args[0]}", err=True)
        sys.exit(1)
    click.echo(f"Construct resolved locally; submitting {mode} run for {gene} ...")

    launcher = HPCLauncher(
        config=config,
        password_provider=_resolve_hpc_password_provider(non_interactive),
    )
    db = RunsDB.local_default()
    result = launcher.launch_run(
        gene=gene,
        construct=construct,
        mode=mode,
        local_output_dir=output_dir,
        cleanup=cleanup,
    )
    db.record_run(
        run_id=result.run_id,
        gene=gene,
        mode=mode,
        construct=construct,
        driver_job_id=result.driver_job_id,
        remote_run_dir=result.remote_run_dir,
        local_reports_dir=result.local_reports_dir,
        status="success" if result.success else "failed",
    )
    click.echo(
        f"Run {result.run_id} -> reports in {result.local_reports_dir} "
        f"({'success' if result.success else 'FAILED'})"
    )


@cli.command("status")
@click.argument("run_id", required=False)
def status_cmd(run_id: str | None) -> None:
    """Show run status (all recent runs by default, or one run if given)."""
    from pegrna_pipeline.agent.db import RunsDB

    db = RunsDB.local_default()
    if run_id:
        try:
            run = db.get(run_id)
        except KeyError:
            click.echo(f"No run {run_id!r} in local history.", err=True)
            sys.exit(1)
        click.echo(f"{run.run_id}  {run.status}  remote={run.remote_run_dir}")
    else:
        for run in db.list_recent(20):
            click.echo(f"{run.run_id:<40} {run.status:<10} {run.gene}/{run.mode}  {run.created_at}")


@cli.command("cancel")
@click.argument("run_id")
def cancel_cmd(run_id: str) -> None:
    """qdel the driver job for a run + mark cancelled in local DB."""
    from pathlib import Path

    import yaml

    from pegrna_pipeline.agent.db import RunsDB
    from pegrna_pipeline.agent.launcher import HPCLauncher

    config_path = Path("~/.pegrna/config.yaml").expanduser()
    config = yaml.safe_load(config_path.read_text())
    db = RunsDB.local_default()
    run = db.get(run_id)
    launcher = HPCLauncher(config=config)
    with launcher._build_ssh() as ssh:
        ssh.exec(f"/opt/pbs/bin/qdel {run.driver_job_id}")
    db.mark_cancelled(run_id)
    click.echo(f"Cancelled {run_id} (driver job {run.driver_job_id}).")


@cli.command("tail")
@click.argument("run_id")
@click.option("--follow/--no-follow", default=True)
def tail_cmd(run_id: str, follow: bool) -> None:
    """Stream the remote driver log for a run. Ctrl-C detaches (does not cancel)."""
    from pathlib import Path

    import yaml

    from pegrna_pipeline.agent.db import RunsDB
    from pegrna_pipeline.agent.launcher import HPCLauncher

    config_path = Path("~/.pegrna/config.yaml").expanduser()
    config = yaml.safe_load(config_path.read_text())
    db = RunsDB.local_default()
    run = db.get(run_id)
    launcher = HPCLauncher(config=config)
    driver_log = f"{run.remote_run_dir}/driver.log"
    cmd = f"tail -f {driver_log}" if follow else f"cat {driver_log}"
    with launcher._build_ssh() as ssh:
        # paramiko exec_command for streaming (bypass the wrapper's buffering)
        client = getattr(ssh, "_client", None)
        if client is None:
            click.echo("SSH client does not expose a paramiko channel for streaming.", err=True)
            sys.exit(2)
        _stdin, stdout, _stderr = client.exec_command(cmd)
        try:
            for line in iter(stdout.readline, ""):
                click.echo(line, nl=False)
        except KeyboardInterrupt:
            click.echo("\n[detached — run continues remotely]")


@cli.command("list")
def list_cmd() -> None:
    """List every run recorded in the local history."""
    from pegrna_pipeline.agent.db import RunsDB

    for run in RunsDB.local_default().list_all():
        click.echo(f"{run.run_id:<40} {run.status:<10} {run.gene}/{run.mode}  {run.created_at}")


@cli.command("fetch")
@click.argument("run_id")
def fetch_cmd(run_id: str) -> None:
    """Download reports for a past run from HPC."""
    from pathlib import Path

    import yaml

    from pegrna_pipeline.agent.db import RunsDB
    from pegrna_pipeline.agent.launcher import HPCLauncher

    config_path = Path("~/.pegrna/config.yaml").expanduser()
    config = yaml.safe_load(config_path.read_text())
    db = RunsDB.local_default()
    run = db.get(run_id)
    launcher = HPCLauncher(config=config)
    launcher.download_reports(
        remote_run_dir=run.remote_run_dir,
        local_reports_dir=run.local_reports_dir,
    )
    click.echo(f"Downloaded reports to {run.local_reports_dir}")


@cli.command("setup")
def setup_cmd() -> None:
    """Interactive first-run config: writes ~/.pegrna/config.yaml."""
    from pathlib import Path

    import yaml

    config_path = Path("~/.pegrna/config.yaml").expanduser()
    config_path.parent.mkdir(parents=True, exist_ok=True)
    if config_path.exists() and not click.confirm(f"{config_path} exists. Overwrite?"):
        return
    cfg = {
        "hostname": click.prompt("HPC hostname", default="login.cx3.hpc.imperial.ac.uk"),
        "username": click.prompt("HPC username"),
        "repo_dir": click.prompt("Remote pegRNA-pipeline install path"),
        "remote_runs_base": click.prompt("Remote runs directory", default="~/pegrna-runs"),
        "entrez_email": click.prompt("NCBI Entrez email"),
        "pridict_install_path": click.prompt("PRIDICT2 install path on HPC"),
        "rscript_path": click.prompt("score_guides.R path on HPC"),
        "conda_env": click.prompt("PRIDICT2 conda env", default="pridict2"),
        "compute_mode": click.prompt(
            "Driver compute mode (local=compute-node qsub, back_ssh=SSH to login)",
            default="local",
            type=click.Choice(["local", "back_ssh"]),
        ),
        "cell_line": click.prompt("Default cell line", default="KOLF2_iPSC"),
        "score_cell_line": click.prompt(
            "PRIDICT2 PE2 score cell line (k562 recommended for KOLF2/iPSC)",
            default="k562",
            type=click.Choice(["k562", "hek293t"]),
        ),
    }
    config_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    click.echo(f"Wrote {config_path}")


if __name__ == "__main__":  # pragma: no cover
    cli()
