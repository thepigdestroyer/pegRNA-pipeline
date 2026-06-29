"""End-to-end driver for the recoded TJ-PE pipeline (fake mode smoke test)."""

from pathlib import Path

from pegrna_pipeline.hpc.recoded_driver import design_recoded


def _demo_inputs():
    return dict(
        spacer="ACGTACGTACGTACGTACGT",
        rtt_aa="MKP",
        cargo="ATGAAACCCGGG",
        genome="ACGATGCATGGCATGACTGA" * 16,
        site_s=150,
        cds_mask=[True] * 320,
        pbs_lengths=(15,),
        pegrna_strand="Rv",
    )


def test_design_recoded_fake_mode_end_to_end(tmp_path):
    out = design_recoded(
        **_demo_inputs(),
        mode="fake",
        output_dir=tmp_path,
        manifest={"genome_build": "T2T-KOLF2.1Jv1.1", "git_commit": "test"},
    )
    assert out.result.pegRNA_sequence
    assert out.result.deletion_size == 0
    # a report file was written
    report = Path(out.report_path)
    assert report.exists()
    text = report.read_text()
    assert "Recoded TJ-PE" in text
    assert "T2T-KOLF2.1Jv1.1" in text


def test_design_recoded_returns_report_string(tmp_path):
    out = design_recoded(**_demo_inputs(), mode="fake", output_dir=tmp_path)
    assert "pegRNA" in out.report_text
    assert out.result.site_s == 150
