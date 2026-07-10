"""End-to-end driver for the recoded TJ-PE pipeline.

Glues the gate-scorer factory, the ranked-pool-walk orchestrator, and the report
writer into one call. This is the entry point the CLI and the HPC job script
invoke. ``mode="fake"`` runs the full design logic with deterministic stand-in
scorers (CI / dry-run / smoke); ``mode="real"`` uses the real models (CX3, GPU,
staged weights) and fails loud if any dependency is missing.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Sequence, Tuple, Union

from pegrna_pipeline.hpc.allele_specificity import find_allele_specific
from pegrna_pipeline.hpc.data_types import DEFAULT_SCORE_CELL_LINE, _score_attr
from pegrna_pipeline.hpc.input_formatter import format_tjpe_insertion_editseq
from pegrna_pipeline.hpc.recoded_tjpe_workflow import (
    NoValidPegRNA,
    RecodedTJPEResult,
    run_recoded_tjpe,
)
from pegrna_pipeline.hpc.report_writer import render_recoded_csv, render_recoded_report
from pegrna_pipeline.hpc.scorer_factory import make_scorers
from pegrna_pipeline.hpc.verification import build_verification_report, render_verification_docx

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
    csv_text: str = ""
    csv_path: Optional[str] = None
    #: Path to the colour-annotated manual-verification DOCX (None if python-docx
    #: is not installed in this environment, or if no output_dir was given).
    verification_docx_path: Optional[str] = None


def _real_scorer_kwargs(
    *,
    uniprot_acc: str,
    cas9: str,
    genome_fasta: Optional[str],
    boltz_cache: Optional[str],
    rinalmo_weights: Optional[str],
    rnaformer_model_dir: Optional[str],
    pam: str,
    max_mismatches: int,
) -> dict:
    """Resolve real-scorer resource paths into ``make_scorers`` kwargs (ADR 0008).

    Each resource is taken from its explicit parameter OR the environment variable
    that ``scripts/hpc/smoke_models.py`` and ``run_recoded_pipeline.pbs`` already
    set (explicit wins), so the design run consumes exactly the staged resources the
    5/5 gate self-check validated:

    * ``genome_fasta`` / ``$PEGRNA_GENOME``      -> off-target genome FASTA (cas-offinder)
    * ``boltz_cache`` / ``$BOLTZ_CACHE``          -> pre-staged Boltz-2 weights (offline node)
    * ``rinalmo_weights`` / ``$RINALMO_SS_WEIGHTS``
    * ``rnaformer_model_dir`` / dir of ``$RNAFORMER_CKPT``

    Leaving the genome unset (no param, no env) yields ``genome=None`` so the
    off-target scorer stays fail-loud -- this never fabricates a path.
    """
    # ``or None`` so a hand-exported empty env var degrades to None uniformly
    # (matching the ``if genome:`` fail-loud guard in scorer_factory) rather than
    # reaching a scorer as an empty-string path.
    genome = genome_fasta or os.environ.get("PEGRNA_GENOME") or None
    boltz = boltz_cache or os.environ.get("BOLTZ_CACHE") or None
    rinalmo = rinalmo_weights or os.environ.get("RINALMO_SS_WEIGHTS") or None
    if rnaformer_model_dir is None:
        ckpt = os.environ.get("RNAFORMER_CKPT")
        rnaformer_model_dir = (os.path.dirname(ckpt) or None) if ckpt else None
    return {
        "uniprot_acc": uniprot_acc,
        "cas9": cas9,
        "genome": genome,
        "pam": pam,
        "max_mismatches": max_mismatches,
        "boltz_cache": boltz,
        "rinalmo_weights": rinalmo,
        "rnaformer_model_dir": rnaformer_model_dir,
    }


def _select_top_pegrna(pegrnas: Sequence[Any]) -> Any:
    """Pick the top pegRNA by PRIDICT PE2 (default cell line), else the first.

    Mirrors the ranking policy used elsewhere (:data:`DEFAULT_SCORE_CELL_LINE` +
    :func:`_score_attr`, i.e. ``pe2_k562``). If a design lacks that score attribute
    the first pegRNA (PRIDICT's own top-ranked row) is used.
    """
    attr = _score_attr(DEFAULT_SCORE_CELL_LINE)
    try:
        return max(pegrnas, key=lambda d: float(getattr(d, attr)))
    except (AttributeError, TypeError, ValueError):
        return pegrnas[0]


def _auto_design_spacer(
    *,
    genome: str,
    site_s: int,
    cargo: str,
    auto_spacer: bool,
    min_flank: int,
    runner: Any,
) -> str:
    """Design the pegRNA#1 spacer at the insert site via the injected PRIDICT2 runner.

    Builds the ``(+cargo)`` insertion editseq at ``site_s``, runs one PRIDICT batch,
    and returns the top pegRNA's spacer. Fail-loud (never guesses a spacer): raises if
    auto-design is disabled, no runner is available, or PRIDICT yields no pegRNA. The
    runner is injected (fake in CI, local subprocess in real/api), so this never
    imports or launches real PRIDICT2 itself.
    """
    if not auto_spacer:
        raise ValueError(
            "no spacer supplied and auto_spacer=False: pass spacer= or set auto_spacer=True"
        )
    if runner is None:
        raise RuntimeError(
            "auto_spacer requested but no PRIDICT runner is available: use mode='fake', "
            "configure PRIDICT2_DIR/PRIDICT2_PYTHON for mode 'real'/'api', or pass "
            "pridict_runner="
        )
    editseq = format_tjpe_insertion_editseq(genome, site_s, cargo, min_flank=min_flank)
    batch = runner.run_batch([{"sequence_name": "peg1", "editseq": editseq}])
    entry = batch.get("peg1") or {}
    pegrnas = entry.get("pegrnas") or []
    if not pegrnas:
        raise RuntimeError(
            f"PRIDICT returned no pegRNA for the auto-spacer pegRNA#1 editseq at "
            f"site_s={site_s}; cannot auto-design a spacer"
        )
    return _select_top_pegrna(pegrnas).spacer


def design_recoded(
    *,
    spacer: Optional[str] = None,
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
    auto_spacer: bool = True,
    pridict_min_flank: int = 40,
    pridict_runner: Optional[Any] = None,
    output_dir: Optional[Union[str, Path]] = None,
    manifest: Optional[dict] = None,
    uniprot_acc: str = "",
    genome_fasta: Optional[str] = None,
    boltz_cache: Optional[str] = None,
    rinalmo_weights: Optional[str] = None,
    rnaformer_model_dir: Optional[str] = None,
    pam: str = "NGG",
    max_mismatches: int = 3,
    protein_length: Optional[int] = None,
    structure_threshold: float = 0.5,
    progress_cb=None,
) -> RecodedDesignOutput:
    """Run the recoded TJ-PE design end-to-end and (optionally) write a report.

    The pegRNA ``spacer`` is optional. When omitted and ``auto_spacer`` is ``True``
    (the default), PRIDICT2 designs pegRNA#1 at the insert site and the pipeline
    adopts its spacer (``pridict_min_flank`` sets the required flank each side of the
    editseq); the ngRNA-nick enumeration, recode arm, and gate cascade are unchanged.
    Passing an explicit ``spacer`` keeps the historical fixed-spacer path and never
    invokes PRIDICT. The runner is injected — an explicit ``pridict_runner`` wins,
    else the scorer bundle's (``PridictRunnerFake`` in fake mode, a local-subprocess
    ``PridictLocalRunner`` in real/api). Auto-spacer fails loud if no runner is
    available or PRIDICT returns no pegRNA.

    Returns the :class:`RecodedTJPEResult`, the rendered report text, and the
    path it was written to (or ``None`` if ``output_dir`` was not supplied).
    Raises :class:`pegrna_pipeline.hpc.recoded_tjpe_workflow.NoValidPegRNA` if no
    candidate survives all gates (fail-loud).

    ``genome`` is the LOCAL design-context sequence (for nick search / recode /
    INV-0). In ``mode="real"`` the heavy gate models also need staged resources,
    which are supplied via the distinct parameters ``genome_fasta`` (off-target
    genome FASTA for cas-offinder), ``boltz_cache``, ``rinalmo_weights`` and
    ``rnaformer_model_dir``; each falls back to the environment variable the smoke
    harness / ``run_recoded_pipeline.pbs`` set (``PEGRNA_GENOME`` / ``BOLTZ_CACHE`` /
    ``RINALMO_SS_WEIGHTS`` / ``RNAFORMER_CKPT``) when the parameter is omitted, so a
    real run consumes exactly the resources the 5/5 gate self-check validated
    (ADR 0008). These are ignored in ``mode="fake"``.
    """
    scorer_kwargs = (
        _real_scorer_kwargs(
            uniprot_acc=uniprot_acc,
            cas9=cas9,
            genome_fasta=genome_fasta,
            boltz_cache=boltz_cache,
            rinalmo_weights=rinalmo_weights,
            rnaformer_model_dir=rnaformer_model_dir,
            pam=pam,
            max_mismatches=max_mismatches,
        )
        if mode in ("real", "api")
        else {}
    )
    scorers = make_scorers(mode=mode, **scorer_kwargs)

    # Auto-spacer (wiring spec 2026-07-10): with no explicit spacer, PRIDICT2 designs
    # pegRNA#1 at the insert site and the pipeline adopts its spacer. An explicit
    # spacer keeps the historical fixed-spacer path unchanged (it never touches the
    # runner). The runner is injected — an explicit ``pridict_runner`` wins, else the
    # scorer bundle's (PridictRunnerFake in fake mode) — so CI never runs real PRIDICT.
    if spacer is None:
        spacer = _auto_design_spacer(
            genome=genome,
            site_s=site_s,
            cargo=cargo,
            auto_spacer=auto_spacer,
            min_flank=pridict_min_flank,
            runner=pridict_runner if pridict_runner is not None else scorers.pridict_runner,
        )

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
        maxent_scorer=scorers.maxent_scorer,
        protein_length=protein_length,
        structure_threshold=structure_threshold,
        on_gate=progress_cb,
    )

    manifest = dict(manifest or {})
    manifest.setdefault("scorer_mode", mode)
    report_text = render_recoded_report(result, manifest=manifest)
    csv_text = render_recoded_csv(result, manifest=manifest)

    report_path: Optional[str] = None
    csv_path: Optional[str] = None
    verification_docx_path: Optional[str] = None
    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        path = out / "recoded_tjpe_report.txt"
        path.write_text(report_text)
        report_path = str(path)
        cpath = out / "recoded_tjpe_components.csv"
        cpath.write_text(csv_text)
        csv_path = str(cpath)
        # Machine-readable per-gate info for the web flowchart (each step's score).
        import json as _json

        (out / "recoded_tjpe_gates.json").write_text(
            _json.dumps(
                {
                    "gate_results": result.gate_results,
                    "recode_edits": result.recode_edits,
                    "deletion_size": result.deletion_size,
                    "candidates_tried": result.candidates_tried,
                    "warnings": result.warnings,
                }
            )
        )
        # Colour-annotated manual-verification DOCX (sections 1-4). Guarded: if
        # python-docx is not installed the run still succeeds with text + CSV.
        try:
            vreport = build_verification_report(result, genome=genome, cargo=cargo)
            vpath = out / "recoded_tjpe_verification.docx"
            render_verification_docx(vreport, vpath)
            verification_docx_path = str(vpath)
        except ImportError:
            pass

    return RecodedDesignOutput(
        result=result,
        report_text=report_text,
        report_path=report_path,
        csv_text=csv_text,
        csv_path=csv_path,
        verification_docx_path=verification_docx_path,
    )


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
