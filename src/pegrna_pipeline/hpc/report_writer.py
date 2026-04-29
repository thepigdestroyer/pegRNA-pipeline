"""Text + CSV report writers for TJ-PE, PA, and PE7 results."""

from __future__ import annotations

import csv
import gzip
from datetime import datetime, timezone
from pathlib import Path
from textwrap import fill
from typing import Any, Dict, List

from pegrna_pipeline.hpc.annotator import (
    annotate_pa_result,
    annotate_pe7_result,
    annotate_tjpe_result,
)
from pegrna_pipeline.hpc.data_types import GeneRecord
from pegrna_pipeline.hpc.offtarget_scorer import OffTargetReport
from pegrna_pipeline.hpc.pa_workflow import PAResult
from pegrna_pipeline.hpc.pe7_workflow import PE7Result
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD, TJPEResult

_BAR = "=" * 72
_DASH = "-" * 72


def _utc_iso_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S") + "Z"


def _assemble_ngrna(spacer: str) -> str:
    """Return the full, orderable ngRNA oligo: spacer + standard scaffold."""
    return spacer + STANDARD_SGRNA_SCAFFOLD


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
- PRIDICT 2.0's nicking-guide selector occasionally emits 19-nt ngRNA
  spacers (one nt short of the canonical 20-nt SpCas9 protospacer). They
  still bind and nick, but DeepSpCas9 / DeepHF on-target models won't
  score them and CFD-based off-target scoring is skipped — those rows
  show as NA in the off-target CSV. The pegRNA + ngRNA pair itself is
  still synthesisable and orderable.
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
    lines.append(f"Generated:      {_utc_iso_now()}")
    lines.append(f"Gene:           {gene.gene_name} ({gene.refseq_accession})")
    lines.append(f"Stop codon:     {gene.stop_codon} at mRNA position {gene.stop_codon_start}")
    lines.append(f"Construct:      {construct_str}")
    lines.append(f"Insert length:  {len(result.full_insert_dna)} bp")
    lines.append(
        f"PE2 {result.pegrna_pe2_cell_line.upper()}:    " f"{result.pegrna_pe2_score:.2f}%"
    )
    lines.append(f"ngRNA distance: {result.ngrna_distance_to_edit:+d} bp from edit")
    lines.append(
        f"Joint PE3:      {result.joint_pe3_score:.3f} "
        f"(PE2 × DeepCas9 / 100; informational, not used for ranking)"
    )
    lines.append("")
    lines.append(_BAR)
    lines.append("WINNING TJ-pegRNA + ngRNA")
    lines.append(_BAR)
    lines.append("")
    ngrna_full = _assemble_ngrna(result.ngrna_spacer)
    lines.append(f"FULL TJ-pegRNA OLIGO ({len(result.tjpe_pegRNA_sequence)} nt):")
    lines.append(_wrap(result.tjpe_pegRNA_sequence))
    lines.append("")
    lines.append(f"FULL ngRNA OLIGO ({len(ngrna_full)} nt):")
    lines.append(_wrap(ngrna_full))
    lines.append("")
    lines.append(_DASH)
    lines.append("COMPONENT BREAKDOWN")
    lines.append(_DASH)
    for ann in annotations:
        lines.append("")
        lines.append(f"{ann.name.upper()} ({ann.length_bp} nt) " f"[source: {ann.source.value}]:")
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


# ----------------------------------------------------------------------- PE7


def write_pe7_txt_report(
    result: PE7Result,
    gene: GeneRecord,
    construct_str: str,
    offtarget_reports: List[OffTargetReport],
    output_path: Path,
) -> None:
    annotations = annotate_pe7_result(result, gene)
    lines: List[str] = [_BAR, "pegRNA DESIGN REPORT -- PE7 (PE3-strategy) mode", _BAR]
    lines.append(f"Generated:      {_utc_iso_now()}")
    lines.append(f"Gene:           {gene.gene_name} ({gene.refseq_accession})")
    lines.append(f"Stop codon:     {gene.stop_codon} at mRNA position {gene.stop_codon_start}")
    lines.append(f"Construct:      {construct_str}")
    lines.append(f"Insert length:  {len(result.full_insert_dna)} bp")
    cl = result.pegrna_pe2_cell_line.upper() or "K562"
    lines.append(f"PE2 {cl}:    {result.pegrna_pe2_score:.2f}%")
    lines.append(f"ngRNA distance: {result.ngrna_distance_to_edit:+d} bp from edit")
    lines.append(
        f"Joint PE3:      {result.joint_pe3_score:.3f} "
        f"(PE2 × DeepCas9 / 100; informational, not used for ranking)"
    )
    lines.append(f"Total ranked pairs: {len(result.top_pairs)}")
    lines.append("")

    # ---------------- TOP-N RANKED PAIRS TABLE ----------------------------
    if result.top_pairs:
        lines.append(_BAR)
        lines.append(f"TOP {len(result.top_pairs)} RANKED pegRNA + ngRNA PAIRS")
        lines.append(_BAR)
        lines.append("")
        header = (
            f"{'rank':>4}  {'pegrna_spacer':<22}  {'strand':>6}  "
            f"{'PE2_'+cl[:5]:>10}  {'ngrna_spacer':<22}  {'dist':>5}  "
            f"{'DeepCas9':>9}  {'jointPE3':>9}"
        )
        lines.append(header)
        lines.append("-" * len(header))
        for p in result.top_pairs:
            lines.append(
                f"{p.rank:>4}  {p.spacer:<22}  {p.pegrna_strand:>6}  "
                f"{p.pegrna_pe2_score:>10.2f}  {p.ngrna_spacer:<22}  "
                f"{p.ngrna_distance_to_edit:>+5d}  "
                f"{p.ngrna_deepcas9_score:>9.2f}  {p.joint_pe3_score:>9.3f}"
            )
        lines.append("")
        lines.append(_BAR)
        lines.append("ORDERABLE SEQUENCES — ALL RANKED PAIRS")
        lines.append(_BAR)
        lines.append("Copy-paste directly into synthesis order form.")
        lines.append("")
        for p in result.top_pairs:
            ngrna_full = _assemble_ngrna(p.ngrna_spacer)
            lines.append(_DASH)
            lines.append(
                f"RANK {p.rank}  |  pegRNA PE2={p.pegrna_pe2_score:.2f}%  "
                f"ngRNA dist={p.ngrna_distance_to_edit:+d} bp  "
                f"DeepCas9={p.ngrna_deepcas9_score:.2f}"
            )
            lines.append(f"  pegRNA ({len(p.pegRNA_full)} nt):  {p.pegRNA_full}")
            lines.append(f"  ngRNA  ({len(ngrna_full)} nt):  {ngrna_full}")
        lines.append(_DASH)
        lines.append("")

    lines.append(_BAR)
    lines.append("RANK 1 — pegRNA + ngRNA DETAIL")
    lines.append(_BAR)
    lines.append("")
    ngrna_full_r1 = _assemble_ngrna(result.ngrna_spacer) if result.ngrna_spacer else ""
    lines.append(f"FULL pegRNA OLIGO ({len(result.pegRNA_full)} nt):")
    lines.append(_wrap(result.pegRNA_full or "<not assembled>"))
    lines.append("")
    lines.append(f"FULL ngRNA OLIGO ({len(ngrna_full_r1)} nt):")
    lines.append(_wrap(ngrna_full_r1 if ngrna_full_r1 else "<not assembled>"))
    lines.append("")
    lines.append(_DASH)
    lines.append("COMPONENT BREAKDOWN")
    lines.append(_DASH)
    for ann in annotations:
        lines.append("")
        lines.append(f"{ann.name.upper()} ({ann.length_bp} nt) " f"[source: {ann.source.value}]:")
        if ann.sequence:
            lines.append(_wrap(f"5'-{ann.sequence}-3'"))
        lines.append(_wrap(ann.functional_role, indent="  "))
        if ann.design_notes:
            lines.append(_wrap(f"Notes: {ann.design_notes}", indent="  "))
    lines.append("")
    lines.append(_DASH)
    lines.append("OFF-TARGET ANALYSIS (every unique spacer across all ranked pairs)")
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
    score_cell_line: str = "k562",
) -> None:
    from pegrna_pipeline.hpc.data_types import _score_attr

    annotations = annotate_pa_result(result, gene, score_cell_line=score_cell_line)
    up = result.winning_pair.upstream
    dn = result.winning_pair.downstream
    score_attr = _score_attr(score_cell_line)
    cl_label = score_cell_line.upper()
    up_score = float(getattr(up, score_attr))
    dn_score = float(getattr(dn, score_attr))
    lines: List[str] = [_BAR, "pegRNA DESIGN REPORT -- Prime Assembly mode", _BAR]
    lines.append(f"Generated:      {_utc_iso_now()}")
    lines.append(f"Gene:           {gene.gene_name} ({gene.refseq_accession})")
    lines.append(f"Stop codon:     {gene.stop_codon} at mRNA position {gene.stop_codon_start}")
    lines.append(f"Construct:      {construct_str}")
    lines.append(
        f"Joint score:    {result.winning_pair.joint_score:.2f} "
        f"(geometric mean of up={up_score:.2f}%, down={dn_score:.2f}% "
        f"PE2 {cl_label})"
    )
    lines.append(f"Nick distance:  {result.winning_pair.nick_distance} bp")
    lines.append(f"Donor length:   {len(result.linear_donor_dna)} bp")
    lines.append("")
    lines.append(_BAR)
    lines.append("WINNING pegRNA PAIR + DONOR")
    lines.append(_BAR)
    lines.append("")
    lines.append(f"UPSTREAM pegRNA ({len(result.upstream_pegRNA_sequence)} nt):")
    lines.append(_wrap(result.upstream_pegRNA_sequence or "<not assembled>"))
    lines.append("")
    lines.append(f"DOWNSTREAM pegRNA ({len(result.downstream_pegRNA_sequence)} nt):")
    lines.append(_wrap(result.downstream_pegRNA_sequence or "<not assembled>"))
    lines.append("")
    lines.append(f"LINEAR DONOR DNA ({len(result.linear_donor_dna)} nt) [order as gBlock]:")
    lines.append(_wrap(result.linear_donor_dna))
    lines.append("")
    lines.append(_DASH)
    lines.append("COMPONENT BREAKDOWN")
    lines.append(_DASH)
    for ann in annotations:
        lines.append("")
        lines.append(f"{ann.name.upper()} ({ann.length_bp} nt) " f"[source: {ann.source.value}]:")
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
    from pegrna_pipeline.hpc.offtarget_scorer import _fmt_float, _fmt_int

    if not reports:
        lines.append("(no off-target reports generated)")
        return
    for r in reports:
        lines.append("")
        lines.append(f"{r.guide_id}: spacer 5'-{r.spacer}-3'")
        lines.append(f"  DeepSpCas9:    {_fmt_float(r.deepspcas9)}")
        lines.append(f"  DeepHF:        {_fmt_float(r.deephf)}")
        lines.append(f"  CFD spec:      {_fmt_float(r.cfd_specificity)}")
        lines.append(
            "  0/1/2/3 mismatch hits: "
            f"{_fmt_int(r.mm0_count)}/{_fmt_int(r.mm1_count)}/"
            f"{_fmt_int(r.mm2_count)}/{_fmt_int(r.mm3_count)}"
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
    "pegrna_pe2_score",
    "pegrna_pe2_cell_line",
    "joint_pe3_score",
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
                "variant_sequence_name": getattr(result.winning_variant, "sequence_name", ""),
                "first_8aa_dna": getattr(result.winning_variant, "codon_variant_dna", ""),
                "full_insert_dna": result.full_insert_dna,
                "pegrna_spacer": result.spacer,
                "pegrna_pbs1": result.pbs1,
                "pegrna_rc_pbs2": result.rc_pbs2,
                "pegrna_rc_pbs2_length": result.rc_pbs2_length,
                "tjpe_pegRNA_full_oligo": result.tjpe_pegRNA_sequence,
                "ngrna_spacer": result.ngrna_spacer,
                "ngrna_distance_to_edit": result.ngrna_distance_to_edit,
                "ngrna_deepcas9_score": f"{result.ngrna_deepcas9_score:.4f}",
                "pegrna_pe2_score": f"{result.pegrna_pe2_score:.4f}",
                "pegrna_pe2_cell_line": result.pegrna_pe2_cell_line,
                "joint_pe3_score": f"{result.joint_pe3_score:.4f}",
                "edit_position_in_mrna": result.edit_position_in_mrna,
                "ngrna_nick_in_mrna": result.ngrna_nick_in_mrna,
                "cytosine_count": dna.count("C"),
                "gc_content": f"{gc:.4f}",
            }
        )


PE7_CSV_COLUMNS = [
    "rank",
    "variant_sequence_name",
    "full_insert_dna",
    "pegrna_spacer",
    "pegrna_pbs",
    "pegrna_rt",
    "pegrna_strand",
    "pegRNA_full",
    "ngrna_spacer",
    "ngrna_distance_to_edit",
    "ngrna_strand",
    "ngrna_deepcas9_score",
    "pegrna_pe2_score",
    "pegrna_pe2_cell_line",
    "joint_pe3_score",
    "edit_position_in_mrna",
    "ngrna_nick_in_mrna",
    "cytosine_count",
    "gc_content",
]


def write_pe7_csv(result: PE7Result, output_path: Path) -> None:
    """Write one row per ranked (pegRNA, ngRNA) pair in ``result.top_pairs``.

    Falls back to a single rank=1 row built from the legacy top-level
    fields when ``top_pairs`` is empty (e.g. older callers that haven't
    populated it yet).
    """
    pairs = list(result.top_pairs)
    if not pairs and result.spacer:
        # Synthesize a single-row pair from the legacy top-level fields.
        from pegrna_pipeline.hpc.pe7_workflow import PE7DesignPair

        pairs = [
            PE7DesignPair(
                rank=1,
                variant_sequence_name=getattr(result.winning_variant, "sequence_name", ""),
                pegRNA_full=result.pegRNA_full,
                spacer=result.spacer,
                pbs=result.pbs,
                rt=result.rt,
                pegrna_strand=result.pegrna_strand,
                pegrna_pe2_score=result.pegrna_pe2_score,
                pegrna_pe2_cell_line=result.pegrna_pe2_cell_line,
                ngrna_spacer=result.ngrna_spacer,
                ngrna_distance_to_edit=result.ngrna_distance_to_edit,
                ngrna_strand=result.ngrna_strand,
                ngrna_deepcas9_score=result.ngrna_deepcas9_score,
                ngrna_nick_in_mrna=result.ngrna_nick_in_mrna,
                joint_pe3_score=result.joint_pe3_score,
            )
        ]
    dna = result.full_insert_dna
    gc = (dna.count("G") + dna.count("C")) / max(len(dna), 1)
    cytosines = dna.count("C")
    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PE7_CSV_COLUMNS)
        writer.writeheader()
        for p in pairs:
            writer.writerow(
                {
                    "rank": p.rank,
                    "variant_sequence_name": p.variant_sequence_name,
                    "full_insert_dna": result.full_insert_dna,
                    "pegrna_spacer": p.spacer,
                    "pegrna_pbs": p.pbs,
                    "pegrna_rt": p.rt,
                    "pegrna_strand": p.pegrna_strand,
                    "pegRNA_full": p.pegRNA_full,
                    "ngrna_spacer": p.ngrna_spacer,
                    "ngrna_distance_to_edit": p.ngrna_distance_to_edit,
                    "ngrna_strand": p.ngrna_strand,
                    "ngrna_deepcas9_score": f"{p.ngrna_deepcas9_score:.4f}",
                    "pegrna_pe2_score": f"{p.pegrna_pe2_score:.4f}",
                    "pegrna_pe2_cell_line": p.pegrna_pe2_cell_line,
                    "joint_pe3_score": f"{p.joint_pe3_score:.4f}",
                    "edit_position_in_mrna": result.edit_position_in_mrna,
                    "ngrna_nick_in_mrna": p.ngrna_nick_in_mrna,
                    "cytosine_count": cytosines,
                    "gc_content": f"{gc:.4f}",
                }
            )
    return None


PA_CSV_COLUMNS = [
    "rank",
    "joint_score",
    "nick_distance",
    "upstream_spacer",
    "upstream_pbs",
    "upstream_rtt",
    "upstream_pe2_score",
    "upstream_target_strand",
    "upstream_editing_position",
    "downstream_spacer",
    "downstream_pbs",
    "downstream_rtt",
    "downstream_pe2_score",
    "downstream_target_strand",
    "downstream_editing_position",
    "upstream_flap_dna",
    "downstream_flap_dna",
    "linear_donor_dna",
    "pe2_cell_line",
    "donor_length_bp",
    "full_integrated_product",
]


def write_pa_csv(result: PAResult, output_path: Path, score_cell_line: str = "k562") -> None:
    from pegrna_pipeline.hpc.data_types import _score_attr

    score_attr = _score_attr(score_cell_line)
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
                    "upstream_pe2_score": f"{float(getattr(up, score_attr)):.4f}",
                    "upstream_target_strand": up.target_strand,
                    "upstream_editing_position": up.editing_position,
                    "downstream_spacer": dn.spacer,
                    "downstream_pbs": dn.pbs,
                    "downstream_rtt": dn.rt,
                    "downstream_pe2_score": f"{float(getattr(dn, score_attr)):.4f}",
                    "downstream_target_strand": dn.target_strand,
                    "downstream_editing_position": dn.editing_position,
                    "upstream_flap_dna": result.upstream_flap_dna if rank == 1 else "",
                    "downstream_flap_dna": (result.downstream_flap_dna if rank == 1 else ""),
                    "linear_donor_dna": result.linear_donor_dna if rank == 1 else "",
                    "pe2_cell_line": score_cell_line if rank == 1 else "",
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
    from pegrna_pipeline.hpc.offtarget_scorer import _fmt_float, _fmt_int

    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OFFTARGET_CSV_COLUMNS)
        writer.writeheader()
        for r in reports:
            writer.writerow(
                {
                    "guide_id": r.guide_id,
                    "spacer": r.spacer,
                    "deepspcas9": _fmt_float(r.deepspcas9, width=4),
                    "deephf": _fmt_float(r.deephf, width=4),
                    "cfd_specificity": _fmt_float(r.cfd_specificity, width=4),
                    "mm0_count": _fmt_int(r.mm0_count),
                    "mm1_count": _fmt_int(r.mm1_count),
                    "mm2_count": _fmt_int(r.mm2_count),
                    "mm3_count": _fmt_int(r.mm3_count),
                }
            )


ALL_VARIANTS_CSV_COLUMNS = [
    "variant_rank",
    "sequence_name",
    "codon_variant_dna",
    "cytosine_count",
    "gc_content",
    "top_pe2_score",
    "top_pe2_cell_line",
    "top_spacer",
    "top_pbs",
    "top_rtt",
    "top_target_strand",
]


def write_all_variants_csv_gz(
    variants: List[Any], output_path: Path, score_cell_line: str = "k562"
) -> None:
    from pegrna_pipeline.hpc.data_types import _score_attr

    score_attr = _score_attr(score_cell_line)
    with gzip.open(output_path, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ALL_VARIANTS_CSV_COLUMNS)
        writer.writeheader()
        for rank, v in enumerate(variants, start=1):
            best = getattr(v, "best_pegRNA", None)
            dna = getattr(v, "codon_variant_dna", "") or ""
            gc = (dna.count("G") + dna.count("C")) / max(len(dna), 1)
            top_score = f"{float(getattr(best, score_attr)):.4f}" if best else ""
            writer.writerow(
                {
                    "variant_rank": rank,
                    "sequence_name": getattr(v, "sequence_name", ""),
                    "codon_variant_dna": dna,
                    "cytosine_count": dna.count("C"),
                    "gc_content": f"{gc:.4f}",
                    "top_pe2_score": top_score,
                    "top_pe2_cell_line": score_cell_line if best else "",
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
    score_cell_line: str = "k562",
) -> Dict[str, Path]:
    """Write text + CSV reports for any result type.

    ``score_cell_line`` controls which PRIDICT PE2 column is shown in the
    text reports and per-variant CSV. PA results carry the choice through
    its annotator and CSV writer; TJ-PE / PE7 results store it on the
    Result type itself, so the choice is internally consistent.
    """
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
        write_pa_txt_report(
            result,
            gene,
            construct_str,
            offtarget_reports,
            txt,
            score_cell_line=score_cell_line,
        )
        write_pa_csv(result, csv_p, score_cell_line=score_cell_line)
    elif isinstance(result, PE7Result):
        txt = reports_dir / f"{base_name}.txt"
        csv_p = reports_dir / f"{base_name}.csv"
        write_pe7_txt_report(result, gene, construct_str, offtarget_reports, txt)
        write_pe7_csv(result, csv_p)
    else:
        raise TypeError(f"unknown result type {type(result).__name__}")
    ot_path = reports_dir / f"{base_name}_offtarget.csv"
    write_offtarget_csv(offtarget_reports, ot_path)
    av_path = reports_dir / f"{base_name}_all_variants.csv.gz"
    write_all_variants_csv_gz(all_evaluated_variants, av_path, score_cell_line=score_cell_line)
    paths.update(txt=txt, csv=csv_p, offtarget=ot_path, all_variants=av_path)
    return paths
