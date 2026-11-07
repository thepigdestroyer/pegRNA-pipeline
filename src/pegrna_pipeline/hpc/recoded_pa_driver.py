"""End-to-end driver for the Prime Assembly design: scorers, the pair walk, and the files.

The counterpart of :func:`recoded_driver.design_recoded` for a cargo that arrives as a separate donor
(``docs/superpowers/specs/2026-11-04-prime-assembly-design.md``). The scorer modes are the same: ``fake`` runs the
full design logic with deterministic stand-ins (CI, dry run), ``real`` and ``api`` use the real models, and a gate
that had to be served by a stand-in is declared in the manifest, so a report can never pass a placeholder off as a
verdict.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Optional, Sequence, Tuple, Union

from pegrna_pipeline.hpc.pa_report import (
    pa_gates_payload,
    render_pa_csv,
    render_pa_fasta,
    render_pa_report,
)
from pegrna_pipeline.hpc.recoded_driver import _attach_splice_context, _real_scorer_kwargs
from pegrna_pipeline.hpc.recoded_pa_workflow import RecodedPAResult, run_recoded_pa
from pegrna_pipeline.hpc.report_writer import render_changes_csv
from pegrna_pipeline.hpc.scorer_factory import make_scorers

__all__ = ["design_pa", "RecodedPADesignOutput", "pa_stand_in_gates", "PA_FILES"]

#: The files a design writes, by role. The web API serves exactly these names.
PA_FILES = {
    "report": "recoded_pa_report.txt",
    "components": "recoded_pa_components.csv",
    "fasta": "recoded_pa_donor.fasta",
    "changes": "recoded_pa_changes.csv",
    "gates": "recoded_pa_gates.json",
}


@dataclass
class RecodedPADesignOutput:
    result: RecodedPAResult
    report_text: str
    report_path: Optional[str]
    csv_text: str = ""
    csv_path: Optional[str] = None
    fasta_text: str = ""
    fasta_path: Optional[str] = None
    stand_in_gates: Tuple[str, ...] = ()


def pa_stand_in_gates(scorers: Any) -> Tuple[str, ...]:
    """The gates of this bundle that a Prime Assembly run uses and that a stand-in is serving.

    The fake bundle also lists ``pridict`` (the auto-spacer source of tj-PE). A Prime Assembly run chooses its own
    guides and never calls it, so naming it in a PA report would warn about a gate that does not exist.
    """
    return tuple(g for g in scorers.stub_gates if g != "pridict")


def design_pa(
    *,
    cargo: str,
    genome: str,
    site_s: int,
    cds_mask: Sequence[bool],
    pbs_lengths: Tuple[int, ...] = (13, 15, 17),
    flap_length: int = 35,
    flap_min: int = 30,
    flap_max: int = 50,
    min_span: int = 30,
    max_span: int = 100,
    max_run: int = 8,
    tail: str = "UUU",
    cas9: str = "nCas9",
    max_offtargets: int = 0,
    max_pairs: int = 30,
    cargo_locked: Optional[Sequence[Sequence[int]]] = None,
    cargo_parts: Optional[Sequence[Any]] = None,
    protein_length: Optional[int] = None,
    structure_threshold: Optional[float] = None,
    mode: str = "fake",
    transcript: Optional[str] = None,
    output_dir: Optional[Union[str, Path]] = None,
    manifest: Optional[dict] = None,
    uniprot_acc: str = "",
    genome_fasta: Optional[str] = None,
    boltz_cache: Optional[str] = None,
    rinalmo_weights: Optional[str] = None,
    rnaformer_model_dir: Optional[str] = None,
    pam: str = "NGG",
    max_mismatches: int = 3,
    progress_cb: Optional[Callable[..., None]] = None,
) -> RecodedPADesignOutput:
    """Design a Prime Assembly insertion end to end and (optionally) write its files.

    Two guides on opposite strands each write a flap, a separate linear donor carries the cargo, and the donor's
    arms restore the DNA between the nicks (recoded where coding, copied exactly where not). Returns the result,
    the report, the components CSV and the FASTA; raises
    :class:`~pegrna_pipeline.hpc.recoded_pa_workflow.NoValidPADesign` (with per-pair diagnostics) when no pair
    survives every gate, or when the cargo is shorter than two flaps.

    ``genome`` is the local design window; ``transcript`` (an Ensembl id) lets the SpliceAI cross-check use the real
    genomic flanks, which is what makes its verdict binding (ADR 0028). ``mode`` and the staged-resource parameters
    behave as in ``design_recoded``. ``max_offtargets`` is the number of off-target sites tolerated per guide (0 by
    default; an exact duplicate elsewhere is never tolerated). ``cargo_locked`` / ``cargo_parts`` are validated and
    reported: a donor's cargo is never recoded, so they only change what the report says was kept exactly.
    """
    if pam != "NGG":
        raise ValueError(f"design-pa supports SpCas9 with an NGG PAM only; got pam={pam!r}")
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
    splice_context = _attach_splice_context(
        scorers.splice_scorer, transcript=transcript, genome=genome
    )

    result = run_recoded_pa(
        cargo=cargo,
        genome=genome,
        site_s=site_s,
        cds_mask=cds_mask,
        splice_scorer=scorers.splice_scorer,
        full_structure_scorer=scorers.full_structure_scorer,
        importance_scorer=scorers.importance_scorer,
        dock_scorer=scorers.dock_scorer,
        offtarget_scorer=scorers.offtarget_scorer,
        maxent_scorer=scorers.maxent_scorer,
        pbs_lengths=pbs_lengths,
        flap_length=flap_length,
        flap_min=flap_min,
        flap_max=flap_max,
        min_span=min_span,
        max_span=max_span,
        max_run=max_run,
        tail=tail,
        cas9=cas9,
        structure_threshold=structure_threshold,
        protein_length=protein_length,
        max_offtargets=max_offtargets,
        max_pairs=max_pairs,
        cargo_locked=cargo_locked,
        cargo_parts=cargo_parts,
        on_gate=progress_cb,
    )

    stand_in = pa_stand_in_gates(scorers)
    manifest = dict(manifest or {})
    manifest.setdefault("scorer_mode", mode)
    if splice_context:
        manifest.setdefault("splice_context", splice_context)
    if stand_in:
        manifest.setdefault("stand_in_gates", ", ".join(stand_in))
    report_text = render_pa_report(result, manifest)
    csv_text = render_pa_csv(result, manifest)
    fasta_text = render_pa_fasta(result)

    report_path = csv_path = fasta_path = None
    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / PA_FILES["report"]).write_text(report_text)
        (out / PA_FILES["components"]).write_text(csv_text)
        (out / PA_FILES["fasta"]).write_text(fasta_text)
        # Every nucleotide the donor's arms change from the window, one row each (ADR 0029).
        (out / PA_FILES["changes"]).write_text(render_changes_csv(result))
        (out / PA_FILES["gates"]).write_text(
            json.dumps(pa_gates_payload(result, stand_in_gates=stand_in, mode=mode))
        )
        report_path = str(out / PA_FILES["report"])
        csv_path = str(out / PA_FILES["components"])
        fasta_path = str(out / PA_FILES["fasta"])

    return RecodedPADesignOutput(
        result=result,
        report_text=report_text,
        report_path=report_path,
        csv_text=csv_text,
        csv_path=csv_path,
        fasta_text=fasta_text,
        fasta_path=fasta_path,
        stand_in_gates=stand_in,
    )
