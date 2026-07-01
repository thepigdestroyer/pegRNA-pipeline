"""Human-CAI scoring + opt-in CAI ranking of synonymous variants (ADR 0009).

``human_cai(dna)`` is the Codon Adaptation Index against the human genome-wide
table (geometric mean of relative adaptiveness ``w = f/f_max_in_family``), and
``enumerate_variants(..., ranking="cai")`` reorders survivors by CAI (descending)
with the existing c_count/GC/dna tiebreak. The default ranking is unchanged.
"""

import math

import pytest

from pegrna_pipeline.hpc.codon_expander import (
    CODONS_BY_AA,
    codon_frequency,
    enumerate_variants,
    human_cai,
)

# --- helpers that derive the expectation from the table, not from the SUT -----


def _max_codon(aa: str) -> str:
    return max(CODONS_BY_AA[aa], key=lambda c: codon_frequency(aa, c))


def _min_codon(aa: str) -> str:
    return min(CODONS_BY_AA[aa], key=lambda c: codon_frequency(aa, c))


def _dna(aa_seq, picker):
    return "".join(picker(aa) for aa in aa_seq)


# --- human_cai ---------------------------------------------------------------


def test_human_cai_in_unit_interval():
    v = human_cai("ATGGCCGACGAG")  # M A D E
    assert 0.0 <= v <= 1.0


def test_human_cai_max_codons_is_one():
    # A sequence built from the most-frequent codon of each family is maximally
    # adapted -> every w == 1 -> CAI == 1.0.
    dna = _dna("ADEFGK", _max_codon)
    assert human_cai(dna) == pytest.approx(1.0)


def test_human_cai_min_codons_below_max_and_positive():
    hi = human_cai(_dna("ADEFGK", _max_codon))
    lo = human_cai(_dna("ADEFGK", _min_codon))
    assert 0.0 < lo < hi


def test_human_cai_monotonic_for_single_aa():
    # For one amino acid with multiple codons, the more-frequent codon scores higher.
    hi = human_cai(_max_codon("L"))
    lo = human_cai(_min_codon("L"))
    assert hi >= lo


def test_human_cai_single_codon_families_excluded_gives_one():
    # Met (ATG) and Trp (TGG) are the only codons in their families -> excluded
    # from the geometric mean -> CAI of an all-M/W sequence is 1.0 by convention.
    assert human_cai("ATGTGG") == pytest.approx(1.0)


def test_human_cai_empty_is_zero():
    assert human_cai("") == 0.0


def test_human_cai_case_insensitive():
    assert human_cai("atggccgac") == pytest.approx(human_cai("ATGGCCGAC"))


def test_human_cai_ragged_length_raises():
    with pytest.raises(ValueError):
        human_cai("ATGG")  # 4 nt, not a whole number of codons at frame 0


def test_human_cai_unknown_codon_raises():
    with pytest.raises(ValueError):
        human_cai("ATGNNN")


def test_human_cai_frame_offset():
    seq = _dna("ADEF", _max_codon)
    assert human_cai("G" + seq, frame=1) == pytest.approx(human_cai(seq))


def test_human_cai_never_raises_math_domain_error_on_rare_codons():
    # The 2-decimal table can round a rare codon to 0.00; the w-floor must keep
    # the log finite rather than raising a math domain error.
    dna = _dna("LLLLRRRR", _min_codon)
    val = human_cai(dna)
    assert math.isfinite(val) and val > 0.0


# --- ranking on enumerate_variants ------------------------------------------


def test_default_ranking_is_gc_c_and_unchanged():
    assert enumerate_variants("MKL", max_variants=20) == enumerate_variants(
        "MKL", max_variants=20, ranking="gc_c"
    )


def test_cai_ranking_orders_by_descending_cai():
    variants = enumerate_variants("LL", max_variants=50, ranking="cai")
    cais = [human_cai(v.dna) for v in variants]
    assert cais == sorted(cais, reverse=True)
    # The top variant is at least as adapted as the bottom one.
    assert cais[0] >= cais[-1]


def test_cai_ranking_same_survivor_set_as_default():
    default = {v.dna for v in enumerate_variants("MKL", max_variants=10_000)}
    cai = {v.dna for v in enumerate_variants("MKL", max_variants=10_000, ranking="cai")}
    assert default == cai


def test_cai_ranking_assigns_sequential_ranks():
    variants = enumerate_variants("LL", max_variants=1000, ranking="cai")
    assert [v.rank for v in variants] == list(range(len(variants)))


def test_invalid_ranking_raises():
    with pytest.raises(ValueError):
        enumerate_variants("MK", max_variants=5, ranking="nonsense")


def test_enumerate_pbs_rtt_variants_uses_cai_ranking():
    from pegrna_pipeline.hpc.recoded_tjpe_workflow import enumerate_pbs_rtt_variants

    out = enumerate_pbs_rtt_variants(
        spacer="ACGTACGTACGTACGTACGT", rtt_aa="LL", pbs_lengths=(13,), max_rtt_variants=50
    )
    # One PBS length, so emission order follows the RTT CAI ranking directly.
    cais = [human_cai(v.rtt_dna) for v in out]
    assert cais == sorted(cais, reverse=True)
