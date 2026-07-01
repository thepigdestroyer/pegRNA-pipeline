"""Output decomposition + CSV export for recoded TJ-PE designs.

Covers: the RTT split into labelled insert (cargo) / replacement gDNA (recoded arm)
sub-parts that reconcatenate to the RTT (both strands, both arm sides); the text
report showing the decomposition; and the component-breakdown CSV written by the
driver alongside the text report.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path

from pegrna_pipeline.hpc.recoded_driver import design_recoded
from pegrna_pipeline.hpc.recoded_tjpe_workflow import RecodedTJPEResult, _decompose_rtt
from pegrna_pipeline.hpc.report_writer import render_recoded_csv, render_recoded_report
from pegrna_pipeline.hpc.strand_orientation import orient_rtt_tail

_CTCF = Path(__file__).resolve().parents[1] / "examples" / "ctcf_cterm" / "config.json"


# --------------------------------------------------------------- _decompose_rtt
def test_decompose_rtt_forward_nick_3prime():
    # Fw + nick 3' of S: sense flap = cargo + arm, RTT = flap (identity).
    cargo, arm = "AAAACCCC", "GGGGTTTT"
    rtt = orient_rtt_tail(cargo + arm, pegrna_strand="Fw")
    parts = _decompose_rtt(rtt, cargo=cargo, arm_sense=arm, pegrna_strand="Fw", nick_gt_s=True)
    assert parts == [("insert", "AAAACCCC"), ("replacement_gdna", "GGGGTTTT")]
    assert "".join(s for _, s in parts) == rtt


def test_decompose_rtt_reverse_reorders_to_reconstruct():
    # Rv RC's the whole flap, reversing segment order; decomposition must still
    # reconstruct the RTT and keep correct labels.
    cargo, arm = "AAAACCCC", "GGGGTTTT"
    for nick_gt_s in (True, False):
        flap = (cargo + arm) if nick_gt_s else (arm + cargo)
        rtt = orient_rtt_tail(flap, pegrna_strand="Rv")
        parts = _decompose_rtt(
            rtt, cargo=cargo, arm_sense=arm, pegrna_strand="Rv", nick_gt_s=nick_gt_s
        )
        assert "".join(s for _, s in parts) == rtt  # reconstructs exactly
        labels = {lbl for lbl, _ in parts}
        assert labels == {"insert", "replacement_gdna"}
        # the "insert" part is the oriented cargo
        insert = dict(parts)["insert"]
        assert insert == orient_rtt_tail(cargo, pegrna_strand="Rv")


def test_decompose_rtt_empty_arm_omitted():
    cargo = "AAAACCCC"
    rtt = orient_rtt_tail(cargo, pegrna_strand="Rv")
    parts = _decompose_rtt(rtt, cargo=cargo, arm_sense="", pegrna_strand="Rv", nick_gt_s=True)
    assert parts == [("insert", rtt)]  # no replacement_gdna row


# --------------------------------------------------------------- report / CSV
def _demo_result():
    cargo, arm = "AAAACCCC", "GGGGTTTT"
    rtt = orient_rtt_tail(cargo + arm, pegrna_strand="Fw")
    parts = _decompose_rtt(rtt, cargo=cargo, arm_sense=arm, pegrna_strand="Fw", nick_gt_s=True)
    return RecodedTJPEResult(
        pegRNA_sequence="SPACER" + "SCAF" + "PB2" + rtt + "PB1" + "UUU",
        spacer="SPACER",
        pbs2="PB2",
        rtt=rtt,
        pbs1="PB1",
        polyu="UUU",
        ngrna_protospacer="NG",
        site_s=10,
        nick_position=20,
        deletion_size=0,
        recode_edits=2,
        candidates_tried=1,
        rtt_components=parts,
        insert_seq=dict(parts)["insert"],
        replacement_gdna=dict(parts)["replacement_gdna"],
    )


def test_report_shows_insert_and_replacement_decomposition():
    text = render_recoded_report(_demo_result(), manifest={})
    assert "insert (cargo)" in text
    assert "replacement gDNA" in text
    assert "tjPE template-jumping portion" in text  # portion vs full-construct note


def test_csv_has_component_rows_with_split_insert():
    csv_text = render_recoded_csv(_demo_result(), manifest={"gene": "X"})
    # parse the non-comment rows
    rows = [r for r in csv.reader(io.StringIO(csv_text)) if r and not r[0].startswith("#")]
    header, data = rows[0], rows[1:]
    assert header == ["index", "component", "role", "part_of", "sequence", "length"]
    names = [r[1] for r in data]
    assert names == ["spacer", "scaffold", "PBS2", "insert", "replacement_gDNA", "PBS1", "PolyU"]
    # insert and replacement_gDNA are distinct labelled rows (not one opaque RTT)
    assert "RTT" not in names
    # metadata is present as comment rows
    assert "# gene" in csv_text and "# recode_edits" in csv_text


def test_driver_writes_both_report_and_csv(tmp_path):
    import inspect
    import json

    cfg = json.loads(_CTCF.read_text())
    valid = set(inspect.signature(design_recoded).parameters)
    params = {k: v for k, v in cfg.items() if not k.startswith("_") and k in valid}
    out = design_recoded(mode="fake", output_dir=str(tmp_path), **params)
    assert out.report_path and Path(out.report_path).exists()
    assert out.csv_path and Path(out.csv_path).exists()
    assert out.csv_path.endswith("recoded_tjpe_components.csv")
    # the design's RTT decomposes and reconcatenates
    assert out.result.insert_seq and out.result.replacement_gdna
    assert "".join(s for _, s in out.result.rtt_components) == out.result.rtt
