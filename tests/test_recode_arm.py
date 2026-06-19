"""Task 3: minimum-edit synonymous recode (pure greedy)."""

from Bio.Seq import Seq

from pegrna_pipeline.hpc.microhomology import longest_common_run
from pegrna_pipeline.hpc.recode_arm import RecodeResult, recode_min_edits


def test_recode_breaks_microhomology_preserving_protein():
    # in-frame CDS span; genomic strand identical to it (worst-case microhomology)
    cds_span = "GAAGAAGAAGAAGAAGAA"  # Glu x6, all GAA
    genome = cds_span
    res = recode_min_edits(cds_span, genome=genome, frame=0, max_run=6)
    assert res is not None
    # protein preserved
    assert str(Seq(res.recoded).translate()) == str(Seq(cds_span).translate())
    # microhomology broken below threshold
    assert longest_common_run(res.recoded, genome) <= 6
    # minimal-ish: at least one edit, fewer than every base
    assert 0 < res.n_edits < len(cds_span)


def test_recode_returns_none_when_no_synonymous_solution():
    # Met/Trp only => zero degeneracy => cannot break a long run
    cds_span = "ATGTGGATGTGG"  # Met-Trp-Met-Trp
    res = recode_min_edits(cds_span, genome=cds_span, frame=0, max_run=3)
    assert res is None


def test_recode_noop_when_already_below_threshold():
    cds_span = "GAAGAAGAA"
    # genome shares nothing -> already below threshold -> zero edits
    res = recode_min_edits(cds_span, genome="TTTTTTTTT", frame=0, max_run=6)
    assert res is not None
    assert res.n_edits == 0
    assert res.recoded == cds_span


def test_recode_result_fields():
    res = recode_min_edits("GAAGAAGAAGAA", genome="GAAGAAGAAGAA", frame=0, max_run=4)
    assert isinstance(res, RecodeResult)
    assert res.max_run_after <= 4
    assert str(Seq(res.recoded).translate()) == "EEEE"


def test_recode_rejects_non_codon_aligned_frame():
    import pytest

    with pytest.raises(ValueError, match="codon-aligned"):
        recode_min_edits("GAAGAA", genome="GAAGAA", frame=1, max_run=3)
