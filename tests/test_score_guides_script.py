"""Static sanity checks for scripts/score_guides.R.

This is a black-box file-content test: we can't execute the R code in
CI without the conda env, but we can ensure it has the right CLI shape
and emits the expected success marker.
"""

from __future__ import annotations

from pathlib import Path


SCRIPT_PATH = Path(__file__).parent.parent / "scripts" / "score_guides.R"


def test_script_exists():
    assert SCRIPT_PATH.exists(), f"missing {SCRIPT_PATH}"


def test_script_uses_optparse():
    text = SCRIPT_PATH.read_text()
    assert "optparse" in text.lower()


def test_script_has_required_cli_flags():
    text = SCRIPT_PATH.read_text()
    for flag in ("--guides", "--nuclease", "--scores", "--genome", "--max-mismatches", "--output"):
        assert flag in text, f"missing {flag}"


def test_script_default_nuclease_is_spcas9_h840a():
    text = SCRIPT_PATH.read_text()
    assert "SpCas9_H840A" in text


def test_script_default_scores_include_required_three():
    text = SCRIPT_PATH.read_text()
    for s in ("DeepSpCas9", "DeepHF", "CFD"):
        assert s in text


def test_script_loads_required_bioconductor_pkgs():
    text = SCRIPT_PATH.read_text()
    for pkg in ("crisprScore", "crisprDesign", "BSgenome.Hsapiens.UCSC.hg38"):
        assert pkg in text, f"missing package {pkg}"


def test_script_emits_success_marker():
    text = SCRIPT_PATH.read_text()
    assert "OFFTARGET_SCORING_DONE" in text


def test_output_tsv_columns():
    text = SCRIPT_PATH.read_text()
    # Columns must be present verbatim somewhere in the header string.
    for col in (
        "guide_id",
        "spacer",
        "deepspcas9",
        "deephf",
        "cfd_specificity",
        "mm0_count",
        "mm1_count",
        "mm2_count",
        "mm3_count",
    ):
        assert col in text, f"missing column {col}"
