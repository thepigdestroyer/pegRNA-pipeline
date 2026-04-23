"""Text + CSV report writers for TJ-PE and PA results."""

from __future__ import annotations

import csv
import gzip
from datetime import datetime
from pathlib import Path
from textwrap import fill
from typing import Any, Dict, List

from pegrna_pipeline.hpc.annotator import (
    annotate_pa_result,
    annotate_tjpe_result,
)
from pegrna_pipeline.hpc.data_types import GeneRecord
from pegrna_pipeline.hpc.offtarget_scorer import OffTargetReport
from pegrna_pipeline.hpc.pa_workflow import PAResult
from pegrna_pipeline.hpc.tjpe_workflow import TJPEResult

_BAR = "=" * 72
_DASH = "-" * 72

_CAVEATS = """\
- PRIDICT 2.0 is validated for inserts <= 15 bp; our probes (TJ-PE 24 bp,
  PA 30 bp default) are extrapolation. Use the scores as relative rankings
  within a probe size, not absolute efficiencies.
- TJ-PE practical range is <= 500 bp; >500 bp use Prime Assembly mode.
- PA nick-to-nick distance window (30-80 bp) is a house rule; verify
  empirically for your locus.
- TJ-pegRNA architecture is [spacer][scaffold][insert][PBS1][RC-PBS2]
  (RC-PBS2 at the 3' end per Zheng 2023 Fig 1a). Cross-check a published
  TJ-pegRNA from the Zheng 2023 supplement before first synthesis.
- Validate cloning + transfection in a control system before scaling up.
"""


def _wrap(text: str, width: int = 72, indent: str = "") -> str:
    return fill(
        text,
        width=width,
        initial_indent=indent,
        subsequent_indent=indent,
        break_long_words=False,
        break_on_hyphens=False,
    )


# --------------------------------------------------------------------- TJ-PE

def write_tjpe_txt_report(
    result: TJPEResult,
    gene: GeneRecord,
    construct_str: str,
    offtarget_reports: List[OffTargetReport],
    output_path: Path,
) -> None:
    annotations = annotate_tjpe_result(result, gene)
    lines: List[str] = [_BAR, "pegRNA DESIGN REPORT -- TJ-PE mode", _BAR]
    lines.append(f"Generated:      {datetime.utcnow().isoformat(timespec='seconds')}Z")
    lines.append(f"Gene:           {gene.gene_name} ({gene.refseq_accession})")
    lines.append(
        f"Stop codon:     {gene.stop_codon} at mRNA position {gene.stop_codon_start}"
    )
    lines.append(f"Construct:      {construct_str}")
    lines.append(f"Insert length:  {len(result.full_insert_dna)} bp")
    lines.append(f"PE2 HEK293T:    {result.pegrna_pe2_hek293t:.2f}%")
    lines.append(
        f"ngRNA distance: {result.ngrna_distance_to_edit:+d} bp from edit"
    )
    lines.append("")
    lines.append(_BAR)
    lines.append("WINNING TJ-pegRNA + ngRNA")
    lines.append(_BAR)
    lines.append("")
    lines.append(f"FULL TJ-pegRNA OLIGO ({len(result.tjpe_pegRNA_sequence)} nt):")
    lines.append(_wrap(result.tjpe_pegRNA_sequence))
    lines.append("")
    lines.append(
        f"FULL ngRNA (spacer only, append standard scaffold): "
        f"5'-{result.ngrna_spacer}-3'"
    )
    lines.append("")
    lines.append(_DASH)
    lines.append("COMPONENT BREAKDOWN")
    lines.append(_DASH)
    for ann in annotations:
        lines.append("")
        lines.append(
            f"{ann.name.upper()} ({ann.length_bp} nt) "
            f"[source: {ann.source.value}]:"
        )
        if ann.sequence:
            lines.append(_wrap(f"5'-{ann.sequence}-3'"))
        lines.append(_wrap(ann.functional_role, indent="  "))
        if ann.design_notes:
            lines.append(_wrap(f"Notes: {ann.design_notes}", indent="  "))
    lines.append("")
    lines.append(_DASH)
    lines.append("OFF-TARGET ANALYSIS")
    lines.append(_DASH)
    _append_offtarget_block(lines, offtarget_reports)
    lines.append("")
    lines.append(_DASH)
    lines.append("CAVEATS")
    lines.append(_DASH)
    lines.append(_CAVEATS)
    lines.append(_BAR)
    lines.append("END OF REPORT")
    lines.append(_BAR)
    output_path.write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------------------ PA

def write_pa_txt_report(
    result: PAResult,
    gene: GeneRecord,
    construct_str: str,
    offtarget_reports: List[OffTargetReport],
    output_path: Path,
) -> None:
    annotations = annotate_pa_result(result, gene)
    up = result.winning_pair.upstream
    dn = result.winning_pair.downstream
    lines: List[str] = [_BAR, "pegRNA DESIGN REPORT -- Prime Assembly mode", _BAR]
    lines.append(f"Generated:      {datetime.utcnow().isoformat(timespec='seconds')}Z")
    lines.append(f"Gene:           {gene.gene_name} ({gene.refseq_accession})")
    lines.append(
        f"Stop codon:     {gene.stop_codon} at mRNA position {gene.stop_codon_start}"
    )
    lines.append(f"Construct:      {construct_str}")
    lines.append(
        f"Joint score:    {result.winning_pair.joint_score:.2f} "
        f"(geometric mean of up={up.pe2_hek293t:.2f}%, down={dn.pe2_hek293t:.2f}%)"
    )
    lines.append(f"Nick distance:  {result.winning_pair.nick_distance} bp")
    lines.append(f"Donor length:   {len(result.linear_donor_dna)} bp")
    lines.append("")
    lines.append(_BAR)
    lines.append("WINNING pegRNA PAIR + DONOR")
    lines.append(_BAR)
    lines.append("")
    lines.append(
        f"UPSTREAM pegRNA ({len(result.upstream_pegRNA_sequence)} nt):"
    )
    lines.append(_wrap(result.upstream_pegRNA_sequence or "<not assembled>"))
    lines.append("")
    lines.append(
        f"DOWNSTREAM pegRNA ({len(result.downstream_pegRNA_sequence)} nt):"
    )
    lines.append(_wrap(result.downstream_pegRNA_sequence or "<not assembled>"))
    lines.append("")
    lines.append(
        f"LINEAR DONOR DNA ({len(result.linear_donor_dna)} nt) [order as gBlock]:"
    )
    lines.append(_wrap(result.linear_donor_dna))
    lines.append("")
    lines.append(_DASH)
    lines.append("COMPONENT BREAKDOWN")
    lines.append(_DASH)
    for ann in annotations:
        lines.append("")
        lines.append(
            f"{ann.name.upper()} ({ann.length_bp} nt) "
            f"[source: {ann.source.value}]:"
        )
        if ann.sequence:
            lines.append(_wrap(f"5'-{ann.sequence}-3'"))
        lines.append(_wrap(ann.functional_role, indent="  "))
        if ann.design_notes:
            lines.append(_wrap(f"Notes: {ann.design_notes}", indent="  "))
    lines.append("")
    lines.append(_DASH)
    lines.append("FULL INTEGRATED PRODUCT (after successful PA in cells)")
    lines.append(_DASH)
    lines.append(_wrap(result.full_integrated_product or "<pending>"))
    lines.append("")
    lines.append(_DASH)
    lines.append("OFF-TARGET ANALYSIS")
    lines.append(_DASH)
    _append_offtarget_block(lines, offtarget_reports)
    lines.append("")
    lines.append(_DASH)
    lines.append("CAVEATS")
    lines.append(_DASH)
    lines.append(_CAVEATS)
    lines.append(_BAR)
    lines.append("END OF REPORT")
    lines.append(_BAR)
    output_path.write_text("\n".join(lines) + "\n")


def _append_offtarget_block(lines: List[str], reports: List[OffTargetReport]) -> None:
    if not reports:
        lines.append("(no off-target reports generated)")
        return
    for r in reports:
        lines.append("")
        lines.append(f"{r.guide_id}: spacer 5'-{r.spacer}-3'")
        lines.append(f"  DeepSpCas9:    {r.deepspcas9:.3f}")
        lines.append(f"  DeepHF:        {r.deephf:.3f}")
        lines.append(f"  CFD spec:      {r.cfd_specificity:.3f}")
        lines.append(
            f"  0/1/2/3 mismatch hits: "
            f"{r.mm0_count}/{r.mm1_count}/{r.mm2_count}/{r.mm3_count}"
        )


# -------------------------------------------------------------------- CSVs

TJPE_CSV_COLUMNS = [
    "rank",
    "variant_sequence_name",
    "first_8aa_dna",
    "full_insert_dna",
    "pegrna_spacer",
    "pegrna_pbs1",
    "pegrna_rc_pbs2",
    "pegrna_rc_pbs2_length",
    "tjpe_pegRNA_full_oligo",
    "ngrna_spacer",
    "ngrna_distance_to_edit",
    "ngrna_deepcas9_score",
    "pegrna_pe2_hek293t",
    "edit_position_in_mrna",
    "ngrna_nick_in_mrna",
    "cytosine_count",
    "gc_content",
]


def write_tjpe_csv(result: TJPEResult, output_path: Path) -> None:
    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TJPE_CSV_COLUMNS)
        writer.writeheader()
        dna = result.full_insert_dna
        gc = (dna.count("G") + dna.count("C")) / max(len(dna), 1)
        writer.writerow(
            {
                "rank": 1,
                "variant_sequence_name": getattr(
                    result.winning_variant, "sequence_name", ""
                ),
                "first_8aa_dna": getattr(
                    result.winning_variant, "codon_variant_dna", ""
                ),
                "full_insert_dna": result.full_insert_dna,
                "pegrna_spacer": result.spacer,
                "pegrna_pbs1": result.pbs1,
                "pegrna_rc_pbs2": result.rc_pbs2,
                "pegrna_rc_pbs2_length": result.rc_pbs2_length,
                "tjpe_pegRNA_full_oligo": result.tjpe_pegRNA_sequence,
                "ngrna_spacer": result.ngrna_spacer,
                "ngrna_distance_to_edit": result.ngrna_distance_to_edit,
                "ngrna_deepcas9_score": f"{result.ngrna_deepcas9_score:.4f}",
                "pegrna_pe2_hek293t": f"{result.pegrna_pe2_hek293t:.4f}",
                "edit_position_in_mrna": result.edit_position_in_mrna,
                "ngrna_nick_in_mrna": result.ngrna_nick_in_mrna,
                "cytosine_count": dna.count("C"),
                "gc_content": f"{gc:.4f}",
            }
        )


PA_CSV_COLUMNS = [
    "rank",
    "joint_score",
    "nick_distance",
    "upstream_spacer",
    "upstream_pbs",
    "upstream_rtt",
    "upstream_pe2_hek293t",
    "upstream_target_strand",
    "upstream_editing_position",
    "downstream_spacer",
    "downstream_pbs",
    "downstream_rtt",
    "downstream_pe2_hek293t",
    "downstream_target_strand",
    "downstream_editing_position",
    "upstream_flap_dna",
    "downstream_flap_dna",
    "linear_donor_dna",
    "donor_length_bp",
    "full_integrated_product",
]


def write_pa_csv(result: PAResult, output_path: Path) -> None:
    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PA_CSV_COLUMNS)
        writer.writeheader()
        pairs = list(result.top_pairs) if result.top_pairs else [result.winning_pair]
        for rank, pair in enumerate(pairs, start=1):
            up = pair.upstream
            dn = pair.downstream
            writer.writerow(
                {
                    "rank": rank,
                    "joint_score": f"{pair.joint_score:.4f}",
                    "nick_distance": pair.nick_distance,
                    "upstream_spacer": up.spacer,
                    "upstream_pbs": up.pbs,
                    "upstream_rtt": up.rt,
                    "upstream_pe2_hek293t": f"{up.pe2_hek293t:.4f}",
                    "upstream_target_strand": up.target_strand,
                    "upstream_editing_position": up.editing_position,
                    "downstream_spacer": dn.spacer,
                    "downstream_pbs": dn.pbs,
                    "downstream_rtt": dn.rt,
                    "downstream_pe2_hek293t": f"{dn.pe2_hek293t:.4f}",
                    "downstream_target_strand": dn.target_strand,
                    "downstream_editing_position": dn.editing_position,
                    "upstream_flap_dna": result.upstream_flap_dna if rank == 1 else "",
                    "downstream_flap_dna": (
                        result.downstream_flap_dna if rank == 1 else ""
                    ),
                    "linear_donor_dna": result.linear_donor_dna if rank == 1 else "",
                    "donor_length_bp": len(result.linear_donor_dna) if rank == 1 else "",
                    "full_integrated_product": (
                        result.full_integrated_product if rank == 1 else ""
                    ),
                }
            )


OFFTARGET_CSV_COLUMNS = [
    "guide_id",
    "spacer",
    "deepspcas9",
    "deephf",
    "cfd_specificity",
    "mm0_count",
    "mm1_count",
    "mm2_count",
    "mm3_count",
]


def write_offtarget_csv(reports: List[OffTargetReport], output_path: Path) -> None:
    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OFFTARGET_CSV_COLUMNS)
        writer.writeheader()
        for r in reports:
            writer.writerow(
                {
                    "guide_id": r.guide_id,
                    "spacer": r.spacer,
                    "deepspcas9": f"{r.deepspcas9:.4f}",
                    "deephf": f"{r.deephf:.4f}",
                    "cfd_specificity": f"{r.cfd_specificity:.4f}",
                    "mm0_count": r.mm0_count,
                    "mm1_count": r.mm1_count,
                    "mm2_count": r.mm2_count,
                    "mm3_count": r.mm3_count,
                }
            )


ALL_VARIANTS_CSV_COLUMNS = [
    "variant_rank",
    "sequence_name",
    "codon_variant_dna",
    "cytosine_count",
    "gc_content",
    "top_pe2_hek293t",
    "top_spacer",
    "top_pbs",
    "top_rtt",
    "top_target_strand",
]


def write_all_variants_csv_gz(variants: List[Any], output_path: Path) -> None:
    with gzip.open(output_path, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ALL_VARIANTS_CSV_COLUMNS)
        writer.writeheader()
        for rank, v in enumerate(variants, start=1):
            best = getattr(v, "best_pegRNA", None)
            dna = getattr(v, "codon_variant_dna", "") or ""
            gc = (dna.count("G") + dna.count("C")) / max(len(dna), 1)
            writer.writerow(
                {
                    "variant_rank": rank,
                    "sequence_name": getattr(v, "sequence_name", ""),
                    "codon_variant_dna": dna,
                    "cytosine_count": dna.count("C"),
                    "gc_content": f"{gc:.4f}",
                    "top_pe2_hek293t": (
                        f"{best.pe2_hek293t:.4f}" if best else ""
                    ),
                    "top_spacer": best.spacer if best else "",
                    "top_pbs": best.pbs if best else "",
                    "top_rtt": best.rt if best else "",
                    "top_target_strand": best.target_strand if best else "",
                }
            )


# ----------------------------------------------------------------- dispatch

def write_reports(
    *,
    result: Any,
    gene: GeneRecord,
    construct_str: str,
    offtarget_reports: List[OffTargetReport],
    all_evaluated_variants: List[Any],
    reports_dir: Path,
    base_name: str,
) -> Dict[str, Path]:
    reports_dir.mkdir(parents=True, exist_ok=True)
    paths: Dict[str, Path] = {}
    if isinstance(result, TJPEResult):
        txt = reports_dir / f"{base_name}.txt"
        csv_p = reports_dir / f"{base_name}.csv"
        write_tjpe_txt_report(result, gene, construct_str, offtarget_reports, txt)
        write_tjpe_csv(result, csv_p)
    elif isinstance(result, PAResult):
        txt = reports_dir / f"{base_name}.txt"
        csv_p = reports_dir / f"{base_name}.csv"
        write_pa_txt_report(result, gene, construct_str, offtarget_reports, txt)
        write_pa_csv(result, csv_p)
    else:
        raise TypeError(f"unknown result type {type(result).__name__}")
    ot_path = reports_dir / f"{base_name}_offtarget.csv"
    write_offtarget_csv(offtarget_reports, ot_path)
    av_path = reports_dir / f"{base_name}_all_variants.csv.gz"
    write_all_variants_csv_gz(all_evaluated_variants, av_path)
    paths.update(txt=txt, csv=csv_p, offtarget=ot_path, all_variants=av_path)
    return paths
