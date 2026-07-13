"""Stage G: the reverse-targeting primary + het-specific sense-window designs run in
parallel when both are needed, and the result is deterministic + identical to the
sequential order (the two designs are independent and write separate output subdirs).
"""

from __future__ import annotations

from pegrna_pipeline.hpc.recoded_driver import design_recoded_with_het_edit

_GENOME = "ACGATGCATGGCATGACTGA" * 16  # 320 nt (the known-passing fake-mode fixture)
_DISCRIMINATING = {"element": "pegrna_pam", "offset": 1, "haplotype": 1, "phased": True}
_NON_DISCRIMINATING = {"element": "spacer", "offset": 0, "haplotype": 1, "phased": True}


def _primary() -> dict:
    return dict(
        spacer="ACGTACGTACGTACGTACGT",
        rtt_aa="MKP",
        cargo="ATGAAACCCGGG",
        genome=_GENOME,
        site_s=150,
        cds_mask=[True] * 320,
        pbs_lengths=(15,),
        pegrna_strand="Rv",
        max_run=8,
    )


def _het_window() -> dict:
    w = _primary()
    w["cds_mask"] = [i < 150 for i in range(len(_GENOME))]  # sense/UTR window
    return w


def test_parallel_het_pair_is_produced(tmp_path):
    out = design_recoded_with_het_edit(
        primary=_primary(),
        het_edit_window=_het_window(),
        het_snps=[_DISCRIMINATING],
        mode="fake",
        output_dir=tmp_path,
    )
    assert out.primary.report_text
    assert out.het_edit is not None
    assert "het-edit" in out.het_edit_note.lower()
    # Both passes ran and wrote their own output subdir (independent, no collision).
    assert (tmp_path / "primary" / "recoded_tjpe_report.txt").exists()
    assert (tmp_path / "het_edit" / "recoded_tjpe_report.txt").exists()


def test_parallel_result_is_deterministic(tmp_path):
    a = design_recoded_with_het_edit(
        primary=_primary(),
        het_edit_window=_het_window(),
        het_snps=[_DISCRIMINATING],
        mode="fake",
        output_dir=tmp_path / "a",
    )
    b = design_recoded_with_het_edit(
        primary=_primary(),
        het_edit_window=_het_window(),
        het_snps=[_DISCRIMINATING],
        mode="fake",
        output_dir=tmp_path / "b",
    )
    assert a.primary.report_text == b.primary.report_text
    assert (a.het_edit is None) == (b.het_edit is None)
    if a.het_edit is not None and b.het_edit is not None:
        assert a.het_edit.report_text == b.het_edit.report_text
    assert a.het_edit_note == b.het_edit_note


def test_primary_only_when_no_het_window(tmp_path):
    out = design_recoded_with_het_edit(
        primary=_primary(), het_edit_window=None, het_snps=[], mode="fake", output_dir=tmp_path
    )
    assert out.primary.report_text
    assert out.het_edit is None
    assert out.het_edit_note == "no het-edit requested"
    assert not (tmp_path / "het_edit").exists()


def test_primary_only_when_no_discriminating_het(tmp_path):
    out = design_recoded_with_het_edit(
        primary=_primary(),
        het_edit_window=_het_window(),
        het_snps=[_NON_DISCRIMINATING],
        mode="fake",
        output_dir=tmp_path,
    )
    assert out.het_edit is None
    assert "no discriminating" in out.het_edit_note.lower()
    # No discriminating het -> the het pass is never dispatched, so no het_edit/ subdir.
    assert not (tmp_path / "het_edit").exists()
