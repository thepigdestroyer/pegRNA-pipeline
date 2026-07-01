"""Task 14: recoded TJ-pegRNA report + provenance manifest."""

from pegrna_pipeline.hpc.recoded_tjpe_workflow import RecodedTJPEResult
from pegrna_pipeline.hpc.report_writer import render_recoded_report


def _result():
    return RecodedTJPEResult(
        pegRNA_sequence="ACGT" * 20,
        spacer="ACGTACGTACGTACGTACGT",
        pbs2="GCATGACTGAACGAT",
        rtt="ATGAAACCCGGG" + "GCATGCAT",
        pbs1="ACGTACGTACGTA",
        polyu="UUU",
        ngrna_protospacer="CGGACGTACGTACGTACGTA",
        site_s=150,
        nick_position=195,
        deletion_size=0,
        recode_edits=5,
        candidates_tried=3,
        warnings=["microhomology advisory: residual run above target"],
        allele_specific_options=[],
    )


def _manifest():
    return {
        "genome_build": "T2T-KOLF2.1Jv1.1",
        "grch38_chain": "kolf2.1j_to_grch38.chain",
        "pridict_version": "2.0",
        "dnachisel_version": "3.2.16",
        "git_commit": "abc1234",
        "seed": "0",
    }


def test_report_includes_design_gates_and_provenance():
    text = render_recoded_report(_result(), manifest=_manifest())
    for token in [
        "T2T-KOLF2.1Jv1.1",
        "deletion",
        "0",
        "recode",
        "PBS2",
        "PolyU",
        "pegRNA",
        "ngRNA",
        "candidates_tried",
        "dnachisel_version",
        "git_commit",
        "microhomology",
    ]:
        assert token in text, f"missing {token!r} in report"


def test_report_notes_no_allele_specific_option():
    text = render_recoded_report(_result(), manifest=_manifest())
    assert "no allele-specific" in text.lower()


def test_report_lists_allele_specific_options_when_present():
    from pegrna_pipeline.hpc.allele_specificity import AlleleSpecificOption

    res = _result()
    res.allele_specific_options = [
        AlleleSpecificOption(
            allele_specific=True,
            side="pegrna",
            position_class="pam",
            strength=1.0,
            haplotype=1,
            note="pegRNA-side: clean single-allele edit",
        )
    ]
    text = render_recoded_report(res, manifest=_manifest())
    assert "allele-specific" in text.lower()
    assert "pam" in text.lower()
    assert "haplotype" in text.lower()
