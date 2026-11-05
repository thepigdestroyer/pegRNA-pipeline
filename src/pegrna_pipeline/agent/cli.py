"""Click-based CLI entry point."""

from __future__ import annotations

import getpass
import json
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
            f"  {i:>3}. {comp.source:<7} {comp.name:<20} "
            f"{'LOCKED ' if comp.locked else 'coding '}aa={comp.aa or '-':<12} dna={comp.dna}"
        )
    click.echo(f"Total aa ({len(construct_obj.aa)}): {construct_obj.aa}")
    click.echo(f"Total dna ({len(construct_obj.dna)}): {construct_obj.dna}")
    click.echo(
        "LOCKED parts are never recoded. Paste into a design config as cargo_parts:\n"
        + json.dumps(construct_obj.cargo_parts())
    )


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

    if mode == "pa":
        # Said before anything is queued: a PBS job would otherwise run for a long time and
        # end in "no compatible pairs", with nothing to say why.
        click.echo(
            "WARNING: Prime Assembly (--mode pa) is a prototype. Against real PRIDICT output its "
            "pair filter can never find a pair (it compares a distance-from-nick as if it were a "
            "coordinate), so this run will report no compatible pairs. See issue #29 / ADR 0027. "
            "Use --mode tjpe for the supported design.",
            err=True,
        )

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


def _derive_het_edit_window(design_cfg: dict, *, manifest: dict, override: object) -> dict:
    """Derive the sense-window het-edit pass config from the primary config.

    The het-edit pass is the forward/sense-targeting allele-aware variant of the
    same edit (same ``genome`` / ``site_s`` / ``cargo`` / ``rtt_aa`` / spacer
    policy). The one structural change is ``cds_mask``: it is recomputed to the
    sense/UTR window -- coding 5' of the insert site, non-coding (byte-exact,
    Decision 15) 3' of it -- which is exactly what the web derives from site S
    (``web/static/index.html``) and matches the known-passing het-edit fixture.
    The pass's ``pegrna_strand`` / ``arm_side`` are inherited from the primary;
    supply an explicit ``het_edit_window`` object in the config to override any
    key rather than have the pipeline guess a frame. ``het_snps`` / ``allele_hunt``
    / ``mode`` / ``output_dir`` are NOT set here -- ``design_recoded_with_het_edit``
    supplies those to the het pass itself (setting them here would collide).
    """
    genome = design_cfg.get("genome")
    site_s = design_cfg.get("site_s")
    if not isinstance(genome, str) or site_s is None:
        raise ValueError(
            "het-edit dual-pass needs 'genome' and 'site_s' to derive the sense-window "
            "cds_mask; supply them or an explicit 'het_edit_window' override in the config"
        )
    window = {k: v for k, v in design_cfg.items() if k not in ("het_snps", "allele_hunt")}
    window["cds_mask"] = [i < int(site_s) for i in range(len(genome))]
    window["manifest"] = manifest
    if isinstance(override, dict):
        window.update(override)  # explicit power-user overrides win
    window.pop("het_snps", None)
    window.pop("allele_hunt", None)
    return window


def _write_het_toplevel(out_root: Path, pair: object) -> None:
    """Replicate the primary sub-design's web artifacts to the top level and record
    the het-edit outcome.

    ``design_recoded_with_het_edit`` writes the primary under ``<out>/primary`` and
    the het pass under ``<out>/het_edit``. The web readers and the download routes
    expect ``<out>/recoded_tjpe_*.{txt,csv,json}`` at the top level, so copy the
    primary's three artifacts up (keeping single-pass behaviour) and drop a
    ``het_edit.json`` summarising whether the het pair was produced + the note.
    """
    import json

    out_root.mkdir(parents=True, exist_ok=True)
    prim_dir = out_root / "primary"
    for name in (
        "recoded_tjpe_report.txt",
        "recoded_tjpe_components.csv",
        "recoded_tjpe_gates.json",
    ):
        src = prim_dir / name
        if src.exists():
            (out_root / name).write_text(src.read_text())
    (out_root / "het_edit.json").write_text(
        json.dumps(
            {
                "requested": True,
                "produced": pair.het_edit is not None,
                "note": pair.het_edit_note,
            }
        )
    )


@cli.command("design-recoded")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--mode", type=click.Choice(["fake", "real", "api"]), default="fake")
@click.option("--output-dir", default="results", type=click.Path())
@click.option(
    "--progress",
    is_flag=True,
    default=False,
    help="Emit @@P NDJSON per-gate progress lines to stdout (for the web streamer).",
)
def design_recoded_cmd(config_path: str, mode: str, output_dir: str, progress: bool) -> None:
    """Run the recoded TJ-PE design pipeline from a JSON config.

    The config supplies the design inputs (rtt_aa, cargo, genome, site_s,
    cds_mask, pegrna_strand, pbs_lengths, het_snps, manifest, ...). The
    ``spacer`` key is OPTIONAL: omit it and PRIDICT2 designs the pegRNA#1
    spacer from the insert location (``site_s``); supply ``spacer`` only to
    override with a fixed guide. ``--mode fake`` runs the full design logic
    with stand-in scorers (dry-run / smoke); ``--mode real`` uses the real
    models (CX3, GPU, staged weights).
    """
    import inspect
    import json
    from pathlib import Path

    from pegrna_pipeline.hpc.recoded_driver import design_recoded

    cfg = json.loads(Path(config_path).read_text())
    cfg.setdefault("pegrna_strand", "Rv")
    # Split provenance/metadata from design parameters: anything starting with "_"
    # (e.g. "_comment") or the known provenance keys is routed into the report
    # manifest instead of being passed as a design kwarg (design_recoded only
    # accepts real design parameters). This lets configs be self-documenting.
    valid_params = set(inspect.signature(design_recoded).parameters)
    manifest = dict(cfg.pop("manifest", {}) or {})
    # Het-edit dual-pass routing (Stage 3c): capture the request flag + optional
    # window override BEFORE the filter loop routes them into the manifest (neither
    # is a design_recoded parameter). The dual pass runs only when het-edit is
    # requested AND at least one het SNP is supplied.
    het_edit_requested = bool(cfg.get("het_edit"))
    het_window_override = cfg.get("het_edit_window")
    for key in [k for k in cfg if k.startswith("_") or k not in valid_params]:
        manifest[key.lstrip("_")] = cfg.pop(key)
    if cfg.get("uniprot_acc"):
        manifest.setdefault("uniprot", cfg["uniprot_acc"])

    progress_cb = None
    if progress:

        def progress_cb(
            cand: int, gate_id: str, info: str, passed: bool, metrics=None
        ) -> None:  # noqa: F811
            sys.stdout.write(
                "@@P "
                + json.dumps(
                    {
                        "cand": cand,
                        "gate": gate_id,
                        "info": info,
                        "passed": passed,
                        "metrics": metrics or {},
                    }
                )
                + "\n"
            )
            sys.stdout.flush()

    het_snps = list(cfg.get("het_snps") or [])
    route_het = het_edit_requested and bool(het_snps)

    try:
        if route_het:
            # Dual pass: antisense primary + an allele-specific sense-window het pass.
            # manifest/progress_cb ride inside the per-pass kwargs (design_recoded
            # accepts them); het_snps/mode/output_dir are supplied by the driver.
            # progress_cb streams the primary only, to avoid duplicate cand indices.
            from pegrna_pipeline.hpc.recoded_driver import design_recoded_with_het_edit

            out_root = Path(output_dir)
            primary_kwargs = {k: v for k, v in cfg.items() if k != "het_snps"}
            primary_kwargs["manifest"] = manifest
            primary_kwargs["progress_cb"] = progress_cb
            het_window = _derive_het_edit_window(
                cfg, manifest=manifest, override=het_window_override
            )
            pair = design_recoded_with_het_edit(
                primary=primary_kwargs,
                het_edit_window=het_window,
                het_snps=het_snps,
                mode=mode,
                output_dir=output_dir,
            )
            _write_het_toplevel(out_root, pair)
            click.echo(pair.primary.report_text)
            if pair.het_edit is not None:
                click.echo("\n===== het-edit (allele-specific sense-window) pair =====\n")
                click.echo(pair.het_edit.report_text)
            click.echo(f"\nhet-edit: {pair.het_edit_note}")
            top_report = out_root / "recoded_tjpe_report.txt"
            if top_report.exists():
                click.echo(f"\nReport written to {top_report}")
        else:
            out = design_recoded(
                mode=mode, output_dir=output_dir, manifest=manifest, progress_cb=progress_cb, **cfg
            )
            click.echo(out.report_text)
            if out.report_path:
                click.echo(f"\nReport written to {out.report_path}")
    except Exception as exc:  # fail loud with a clear message
        click.echo(f"design-recoded failed: {type(exc).__name__}: {exc}", err=True)
        # Record which gates were real even though the design failed. "No
        # candidate survived" means something REJECTED the design, so which
        # gates were real models and which were placeholders is exactly what the
        # reader needs -- a rejection by a stand-in is not evidence of anything.
        try:
            from pegrna_pipeline.hpc.scorer_factory import make_scorers  # noqa: PLC0415

            # Same inputs the driver uses to decide which gates can be real:
            # a genome FASTA enables cas-offinder, a key enables Boltz-2.
            probe = (
                {
                    "uniprot_acc": cfg.get("uniprot_acc") or manifest.get("uniprot") or "",
                    "genome": os.environ.get("PEGRNA_GENOME") or None,
                }
                if mode in ("real", "api")
                else {}
            )
            stubs = list(make_scorers(mode=mode, **probe).stub_gates)
            out_dir = Path(output_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "recoded_tjpe_gates.json").write_text(
                json.dumps(
                    {
                        "gate_results": [],
                        "failed": True,
                        "error": f"{type(exc).__name__}: {exc}",
                        "stand_in_gates": stubs,
                        "scorer_mode": mode,
                    }
                )
            )
        except Exception:  # noqa: BLE001 - diagnostics must never mask the real failure
            pass
        sys.exit(1)


if __name__ == "__main__":  # pragma: no cover
    cli()
