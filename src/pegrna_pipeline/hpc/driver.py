"""Compute-node driver that runs end-to-end for one gene + construct.

Invoked by :file:`scripts/driver_pbs_wrapper.sh` on a PBS worker, this
module:

1. Reads ``run_dir/inputs.json`` (gene record + mode + flap parameters).
2. Reads the resolved construct from ``run_dir/construct.yaml``.
3. For TJ-PE: enumerates the first 8 amino acids of the construct,
   formats PRIDICT inputs, submits PRIDICT batches, and runs
   :func:`run_tjpe_workflow`.
4. For PA: enumerates upstream (first 10 aa) AND downstream (last 10
   aa) as two separate PRIDICT batches, optimises the bulk, assembles
   the donor, then runs :func:`run_pa_workflow`.
5. Scores the winning guide(s) against the genome off-target database.
6. Calls :func:`write_reports` (Plan 4 — stubbed here).

All expensive operations go through an ``ssh``-like client returned by
:func:`_local_compute_ssh`, which is either a
:class:`LocalShellClient` (plain subprocess, when the driver runs on a
compute node with direct filesystem access) or an :class:`SSHClient`
pointing back at the login node.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

import yaml

from pegrna_pipeline.hpc.data_types import (
    GeneRecord,
    NgRNADesign,
    PegRNADesign,
    VariantScore,
)
from pegrna_pipeline.hpc.input_formatter import format_bracket_insertion
from pegrna_pipeline.hpc.pa_workflow import PAResult, run_pa_workflow
from pegrna_pipeline.hpc.pridict_parser import parse_ngrna_csv, parse_pegrna_csv
from pegrna_pipeline.hpc.pridict_submitter import (
    BatchResult,
    submit_pridict_batches,
)
from pegrna_pipeline.hpc.tjpe_workflow import TJPEResult, run_tjpe_workflow


# --------------------------------------------------------------- LocalShellClient
class LocalShellClient:
    """Subprocess-backed stand-in for :class:`SSHClient`.

    Used when the driver runs on the same node it will dispatch jobs
    from: everything happens via plain ``subprocess`` + :func:`shutil.copy`.
    The API mirrors :class:`SSHClient` so higher-level helpers
    (``submit_pridict_batches`` etc.) don't care which it gets.
    """

    def __init__(self) -> None:
        self._connected: bool = False

    @property
    def is_connected(self) -> bool:
        return self._connected

    def connect(self) -> None:
        self._connected = True

    def disconnect(self) -> None:
        self._connected = False

    def __enter__(self) -> "LocalShellClient":
        self.connect()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.disconnect()

    def exec(self, command: str):
        proc = subprocess.run(
            ["bash", "-lc", command],
            capture_output=True,
            text=True,
            check=False,
        )
        return proc.stdout, proc.stderr, proc.returncode

    def upload(self, local: str, remote: str) -> None:
        Path(remote).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(local, remote)

    def download(self, remote: str, local: str) -> None:
        Path(local).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(remote, local)


# --------------------------------------------------------------- ssh factory
def _local_compute_ssh(config: Dict[str, Any]):
    """Return either a LocalShellClient or a back-SSH SSHClient.

    ``config["compute_mode"]`` selects:

    * ``"local"`` — :class:`LocalShellClient`.
    * ``"back_ssh"`` — :class:`SSHClient` targeted at
      ``config["hostname"]`` / ``config["username"]``. ``config["password"]``
      provides the password (as the driver runs on a compute node with
      limited prompting ability).
    """
    mode = config.get("compute_mode", "local")
    if mode == "local":
        return LocalShellClient()
    if mode == "back_ssh":
        from pegrna_pipeline.hpc.ssh_client import SSHClient

        pw = config.get("password", "")

        def _provider():
            return pw

        return SSHClient(
            hostname=config["hostname"],
            username=config["username"],
            password_provider=_provider,
        )
    raise ValueError(f"unknown compute_mode {mode!r}")


# --------------------------------------------------------------- batch -> scores
def _build_variant_scores_from_batch(
    batch_result: BatchResult,
    variants: Sequence[Dict[str, Any]],
) -> List[VariantScore]:
    """Walk ``batch_result.output_dirs`` and parse the PRIDICT CSVs.

    The ``variants`` argument is used only for ordering / sanity — one
    :class:`VariantScore` is returned per ``(output_dir, sequence_name)``
    pair found on disk.
    """
    out: List[VariantScore] = []
    for out_dir in batch_result.output_dirs:
        out_dir_path = Path(out_dir)
        if not out_dir_path.exists():
            continue
        for peg_path in sorted(out_dir_path.glob("*_pegRNA_Pridict_full.csv")):
            prefix = peg_path.name.replace("_pegRNA_Pridict_full.csv", "")
            ng_path = out_dir_path / f"{prefix}_nicking_guides.csv"
            peg_designs: List[PegRNADesign] = []
            ng_designs: List[NgRNADesign] = []
            try:
                peg_designs = parse_pegrna_csv(peg_path)
            except Exception:  # pragma: no cover - defensive
                peg_designs = []
            if ng_path.exists():
                try:
                    ng_designs = parse_ngrna_csv(ng_path)
                except Exception:  # pragma: no cover - defensive
                    ng_designs = []
            out.append(
                VariantScore(
                    sequence_name=prefix,
                    codon_variant_dna="",
                    pegRNA_designs=peg_designs,
                    ngRNA_candidates=ng_designs,
                )
            )
    return out


# --------------------------------------------------------------- winners -> guides
def _extract_winning_guides(
    result: Union[TJPEResult, PAResult],
) -> List[Dict[str, str]]:
    """Return a list of ``{"guide_id", "spacer"}`` dicts for off-target scoring."""
    out: List[Dict[str, str]] = []
    if isinstance(result, TJPEResult):
        if result.spacer:
            out.append({"guide_id": "pegrna_spacer", "spacer": result.spacer})
        if result.ngrna_spacer:
            out.append({"guide_id": "ngrna_spacer", "spacer": result.ngrna_spacer})
        return out
    if isinstance(result, PAResult):
        if result.winning_pair is None:
            return []
        out.append(
            {
                "guide_id": "upstream_pegrna_spacer",
                "spacer": result.winning_pair.upstream.spacer,
            }
        )
        out.append(
            {
                "guide_id": "downstream_pegrna_spacer",
                "spacer": result.winning_pair.downstream.spacer,
            }
        )
        return out
    raise TypeError(f"unsupported result type {type(result).__name__}")


# --------------------------------------------------------------- write_reports stub
def write_reports(
    *,
    run_dir: Path,
    gene: GeneRecord,
    result: Union[TJPEResult, PAResult],
    offtarget_reports: Sequence[Any],
) -> Path:
    """Plan 4 will implement this properly — for now drop a text placeholder.

    Returns the path of the primary report file.
    """
    reports_dir = run_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    base_name = gene.gene_name
    report_path = reports_dir / f"{base_name}.txt"
    report_path.write_text(
        f"pegRNA-pipeline report placeholder (Plan 4 will expand this)\n"
        f"Gene: {gene.gene_name} ({gene.refseq_accession})\n"
        f"Result type: {type(result).__name__}\n"
        f"Off-target reports: {len(offtarget_reports)}\n"
    )
    return report_path


# --------------------------------------------------------------- enumerate helpers
def _enumerate_first_aa(aa: str, n: int) -> List[str]:
    """Return codon variants for the first ``n`` aa of ``aa``.

    Thin wrapper around :func:`enumerate_variants`.
    """
    from pegrna_pipeline.hpc.codon_expander import enumerate_variants

    window = aa[:n]
    return [v.dna for v in enumerate_variants(window, max_variants=100)]


def _enumerate_last_aa(aa: str, n: int) -> List[str]:
    from pegrna_pipeline.hpc.codon_expander import enumerate_variants

    window = aa[-n:]
    return [v.dna for v in enumerate_variants(window, max_variants=100)]


# --------------------------------------------------------------- driver core
@dataclass
class DriverConfig:
    """Parsed ``inputs.json`` + ``construct.yaml`` shape."""

    gene: GeneRecord
    mode: str  # "tjpe" | "pa"
    construct_aa: str
    construct_dna: str
    cell_line: str
    pridict_install_path: str
    conda_env: str
    rscript_path: str
    compute_mode: str
    back_ssh: Dict[str, Any]


def _load_inputs(run_dir: Path) -> DriverConfig:
    inputs = json.loads((run_dir / "inputs.json").read_text())
    construct = yaml.safe_load((run_dir / "construct.yaml").read_text())

    gene_dict = inputs["gene"]
    gene = GeneRecord(
        gene_name=gene_dict["gene_name"],
        refseq_accession=gene_dict["refseq_accession"],
        mrna_sequence=gene_dict["mrna_sequence"],
        cds_start=gene_dict["cds_start"],
        cds_end=gene_dict["cds_end"],
        stop_codon_start=gene_dict["stop_codon_start"],
        stop_codon=gene_dict["stop_codon"],
        upstream_flank=gene_dict["upstream_flank"],
        downstream_flank=gene_dict["downstream_flank"],
        flanking_bp=gene_dict["flanking_bp"],
    )
    return DriverConfig(
        gene=gene,
        mode=inputs["mode"],
        construct_aa=construct["aa"],
        construct_dna=construct["dna"],
        cell_line=inputs.get("cell_line", "KOLF2_iPSC"),
        pridict_install_path=inputs["pridict_install_path"],
        conda_env=inputs.get("conda_env", "pridict2"),
        rscript_path=inputs["rscript_path"],
        compute_mode=inputs.get("compute_mode", "local"),
        back_ssh=inputs.get("back_ssh", {}),
    )


def _pridict_variants_for_window(
    gene: GeneRecord,
    dna_variants: Sequence[str],
    *,
    name_prefix: str,
) -> List[Dict[str, str]]:
    """Emit PRIDICT batch rows for every DNA variant."""
    out: List[Dict[str, str]] = []
    for i, dna in enumerate(dna_variants):
        out.append(
            {
                "sequence_name": f"{name_prefix}_{i}",
                "editseq": format_bracket_insertion(gene, dna),
            }
        )
    return out


def run_driver(run_dir: Union[str, Path]) -> Union[TJPEResult, PAResult]:
    """Run the full compute-node pipeline for ``run_dir``.

    Writes a report into ``run_dir/reports/<gene>.txt`` (Plan 4 will
    expand this). Returns the :class:`TJPEResult` or :class:`PAResult`
    so callers (the integration test) can assert on the shape.
    """
    run_path = Path(run_dir)
    cfg = _load_inputs(run_path)
    ssh_cfg = {"compute_mode": cfg.compute_mode, **cfg.back_ssh}

    with _local_compute_ssh(ssh_cfg) as ssh:
        if cfg.mode == "tjpe":
            dna_variants = _enumerate_first_aa(cfg.construct_aa, 8)
            variants = _pridict_variants_for_window(
                cfg.gene, dna_variants, name_prefix="tjpe_var"
            )
            batch = submit_pridict_batches(
                ssh=ssh,
                variants=variants,
                remote_run_dir=str(run_path),
                pridict_install_path=cfg.pridict_install_path,
                conda_env=cfg.conda_env,
            )
            variant_scores = _build_variant_scores_from_batch(batch, variants)
            result: Union[TJPEResult, PAResult] = run_tjpe_workflow(
                gene=cfg.gene,
                insert_full_dna=cfg.construct_dna,
                variant_scores=variant_scores,
            )
        elif cfg.mode == "pa":
            from pegrna_pipeline.hpc.codon_expander import optimise_bulk
            from pegrna_pipeline.hpc.pa_workflow import (
                assemble_linear_donor,
                compute_flap_boundaries,
            )

            up_aa, bulk_aa, down_aa = compute_flap_boundaries(
                cfg.construct_aa, flap_length_aa=10
            )
            bulk_dna = optimise_bulk(bulk_aa)

            up_variants_dna = _enumerate_first_aa(cfg.construct_aa, 10)
            down_variants_dna = _enumerate_last_aa(cfg.construct_aa, 10)

            upstream_variants = _pridict_variants_for_window(
                cfg.gene, up_variants_dna, name_prefix="pa_up"
            )
            downstream_variants = _pridict_variants_for_window(
                cfg.gene, down_variants_dna, name_prefix="pa_down"
            )

            up_batch = submit_pridict_batches(
                ssh=ssh,
                variants=upstream_variants,
                remote_run_dir=str(run_path / "upstream"),
                pridict_install_path=cfg.pridict_install_path,
                conda_env=cfg.conda_env,
            )
            down_batch = submit_pridict_batches(
                ssh=ssh,
                variants=downstream_variants,
                remote_run_dir=str(run_path / "downstream"),
                pridict_install_path=cfg.pridict_install_path,
                conda_env=cfg.conda_env,
            )
            up_scores = _build_variant_scores_from_batch(up_batch, upstream_variants)
            down_scores = _build_variant_scores_from_batch(
                down_batch, downstream_variants
            )

            # Assemble donor from the winning flaps + optimised bulk
            # (we still pass both here so the workflow can validate).
            upstream_flap_dna = up_variants_dna[0] if up_variants_dna else "A"
            downstream_flap_dna = down_variants_dna[0] if down_variants_dna else "A"
            linear_donor = assemble_linear_donor(
                upstream_flap_dna=upstream_flap_dna,
                bulk_dna=bulk_dna,
                downstream_flap_dna=downstream_flap_dna,
            )
            result = run_pa_workflow(
                gene=cfg.gene,
                upstream_variant_scores=up_scores,
                downstream_variant_scores=down_scores,
                upstream_flap_dna=upstream_flap_dna,
                downstream_flap_dna=downstream_flap_dna,
                bulk_dna=bulk_dna,
                donor_dna=linear_donor,
            )
        else:
            raise ValueError(f"unknown mode {cfg.mode!r}")

        # Off-target score the winning guides.
        guides = _extract_winning_guides(result)
        offtarget_reports: List[Any] = []
        if guides:
            try:
                from pegrna_pipeline.hpc.offtarget_scorer import (
                    submit_offtarget_job,
                )

                offtarget_reports = list(
                    submit_offtarget_job(
                        ssh=ssh,
                        guides=guides,
                        remote_run_dir=str(run_path),
                        rscript_path=cfg.rscript_path,
                        conda_env="pegrna-ngrna",
                    )
                )
            except Exception as exc:  # pragma: no cover - integration-only
                print(f"offtarget scoring failed: {exc}", file=sys.stderr)

        write_reports(
            run_dir=run_path,
            gene=cfg.gene,
            result=result,
            offtarget_reports=offtarget_reports,
        )

    # Signal success for the PBS wrapper
    (run_path / "DRIVER_DONE").write_text("ok\n")
    return result


# --------------------------------------------------------------- entry point
def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m pegrna_pipeline.hpc.driver",
        description="Compute-node end-to-end driver.",
    )
    parser.add_argument("run_dir", help="Run directory containing inputs.json + construct.yaml")
    args = parser.parse_args(argv)
    run_driver(args.run_dir)
    return 0


if __name__ == "__main__":  # pragma: no cover — exercised via subprocess/PBS
    raise SystemExit(main())
