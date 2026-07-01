"""End-to-end driver for the recoded TJ-PE pipeline.

Glues the gate-scorer factory, the ranked-pool-walk orchestrator, and the report
writer into one call. This is the entry point the CLI and the HPC job script
invoke. ``mode="fake"`` runs the full design logic with deterministic stand-in
scorers (CI / dry-run / smoke); ``mode="real"`` uses the real models (CX3, GPU,
staged weights) and fails loud if any dependency is missing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Tuple, Union

from pegrna_pipeline.hpc.allele_specificity import find_allele_specific
from pegrna_pipeline.hpc.recoded_tjpe_workflow import (
    NoValidPegRNA,
    RecodedTJPEResult,
    run_recoded_tjpe,
)
from pegrna_pipeline.hpc.report_writer import render_recoded_report
from pegrna_pipeline.hpc.scorer_factory import make_scorers

__all__ = [
    "design_recoded",
    "RecodedDesignOutput",
    "design_recoded_with_het_edit",
    "RecodedPairOutput",
]


@dataclass
class RecodedDesignOutput:
    result: RecodedTJPEResult
    report_text: str
    report_path: Optional[str]


def design_recoded(
    *,
    spacer: str,
    rtt_aa: str,
    cargo: str,
    genome: str,
    site_s: int,
    cds_mask: Sequence[bool],
    pbs_lengths: Tuple[int, ...] = (13, 15, 17),
    ngrna_distance_min: int = 40,
    ngrna_distance_max: int = 70,
    pegrna_strand: str = "Rv",
    arm_side: str = "3prime",
    max_run: int = 8,
    cas9: str = "nCas9",
    allele_hunt: bool = False,
    het_snps: Optional[Sequence[dict]] = None,
    mode: str = "fake",
    output_dir: Optional[Union[str, Path]] = None,
    manifest: Optional[dict] = None,
    uniprot_acc: str = "",
) -> RecodedDesignOutput:
    """Run the recoded TJ-PE design end-to-end and (optionally) write a report.

    Returns the :class:`RecodedTJPEResult`, the rendered report text, and the
    path it was written to (or ``None`` if ``output_dir`` was not supplied).
    Raises :class:`pegrna_pipeline.hpc.recoded_tjpe_workflow.NoValidPegRNA` if no
    candidate survives all gates (fail-loud).
    """
    scorer_kwargs = {"uniprot_acc": uniprot_acc, "cas9": cas9} if mode == "real" else {}
    scorers = make_scorers(mode=mode, **scorer_kwargs)

    result = run_recoded_tjpe(
        spacer=spacer,
        rtt_aa=rtt_aa,
        cargo=cargo,
        genome=genome,
        site_s=site_s,
        cds_mask=cds_mask,
        pbs_lengths=pbs_lengths,
        ngrna_distance_min=ngrna_distance_min,
        ngrna_distance_max=ngrna_distance_max,
        pegrna_strand=pegrna_strand,
        arm_side=arm_side,
        max_run=max_run,
        cas9=cas9,
        allele_hunt=allele_hunt,
        het_snps=het_snps,
        splice_scorer=scorers.splice_scorer,
        full_structure_scorer=scorers.full_structure_scorer,
        importance_scorer=scorers.importance_scorer,
        dock_scorer=scorers.dock_scorer,
        offtarget_scorer=scorers.offtarget_scorer,
    )

    manifest = dict(manifest or {})
    manifest.setdefault("scorer_mode", mode)
    report_text = render_recoded_report(result, manifest=manifest)

    report_path: Optional[str] = None
    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        path = out / "recoded_tjpe_report.txt"
        path.write_text(report_text)
        report_path = str(path)

    return RecodedDesignOutput(result=result, report_text=report_text, report_path=report_path)


@dataclass
class RecodedPairOutput:
    """Primary (antisense recode-into-CDS) pair plus an optional sense-window
    het-edit pair, with a note explaining why the het-edit pair is or is not
    present."""

    primary: RecodedDesignOutput
    het_edit: Optional[RecodedDesignOutput]
    het_edit_note: str


def design_recoded_with_het_edit(
    *,
    primary: dict,
    het_edit_window: Optional[dict] = None,
    het_snps: Optional[Sequence[dict]] = None,
    mode: str = "fake",
    output_dir: Optional[Union[str, Path]] = None,
) -> RecodedPairOutput:
    """Run the antisense primary pass and, when a discriminating phased het is
    present and a sense-window config is supplied, a second het-edit pass.

    ``primary`` and ``het_edit_window`` are kwarg dicts for :func:`design_recoded`;
    each carries its own ``genome``/``cds_mask``/``pegrna_strand``/``cargo`` frame.
    No implicit reverse complement is derived here -- the caller controls strand,
    so a frame error fails loudly in the per-pass TPRT/INV-0 checks rather than
    silently flipping. The het-edit pass uses the byte-exact non-coding
    replacement path (Decision 15) and is allele-specific only because of a
    pre-existing phased het (specificity is capped by natural heterozygosity and
    is never manufactured). Returns both outputs and a note.
    """
    het_snps = list(het_snps or [])
    prim = design_recoded(
        **primary,
        het_snps=het_snps,
        mode=mode,
        output_dir=(Path(output_dir) / "primary") if output_dir else None,
    )
    if het_edit_window is None:
        return RecodedPairOutput(primary=prim, het_edit=None, het_edit_note="no het-edit requested")
    opts = find_allele_specific(het_snps)
    if not opts:
        return RecodedPairOutput(
            primary=prim,
            het_edit=None,
            het_edit_note=(
                "no discriminating phased het -> no allele-specific sense-window pair "
                "(allele specificity is capped by natural heterozygosity, never manufactured)"
            ),
        )
    # The het-edit pass is advisory: a fail-loud here (no surviving sense-window
    # candidate) must NOT kill the successful primary -- report it as a note.
    o = opts[0]
    try:
        het = design_recoded(
            **het_edit_window,
            het_snps=het_snps,
            allele_hunt=True,
            mode=mode,
            output_dir=(Path(output_dir) / "het_edit") if output_dir else None,
        )
    except NoValidPegRNA as exc:
        return RecodedPairOutput(
            primary=prim,
            het_edit=None,
            het_edit_note=(
                f"discriminating {o.side}-side {o.position_class} het present, but no "
                f"sense-window candidate survived the gates: {exc}"
            ),
        )
    note = (
        f"het-edit sense-window pair via {o.side}-side {o.position_class} het "
        f"(discrimination strength {o.strength:.2f}); non-coding re-templated byte-exact"
    )
    return RecodedPairOutput(primary=prim, het_edit=het, het_edit_note=note)
