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


def test_script_default_nuclease_is_spcas9():
    """Default nuclease is the canonical SpCas9 (predefined in crisprBase).

    Earlier versions used the bare string SpCas9_H840A which crisprDesign
    couldn't resolve — recognition is identical so SpCas9 is the right
    default. The CLI still accepts --nuclease overrides.
    """
    text = SCRIPT_PATH.read_text()
    assert 'default = "SpCas9"' in text


def test_script_loads_required_bioconductor_pkgs():
    text = SCRIPT_PATH.read_text()
    # We dropped crisprDesign / crisprScore in favour of Biostrings-only
    # off-target counting (no precomputed BWA/bowtie index needed); the
    # one mandatory package is BSgenome.Hsapiens.UCSC.hg38.
    for pkg in ("Biostrings", "BSgenome.Hsapiens.UCSC.hg38"):
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
