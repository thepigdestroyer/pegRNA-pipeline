"""Tests for core codon tables and lookup helpers."""

import pytest

from pegrna_pipeline.hpc.codon_expander import (
    CODONS_BY_AA,
    HUMAN_CODON_FREQUENCIES,
    codon_frequency,
    codons_for_aa,
)


def test_codons_by_aa_covers_twenty_amino_acids_plus_stop():
    # Standard 20 aa plus stop (represented as '*').
    expected_letters = set("ACDEFGHIKLMNPQRSTVWY*")
    assert expected_letters <= set(CODONS_BY_AA)


def test_codons_for_aa_case_insensitive():
    assert codons_for_aa("A") == codons_for_aa("a")


def test_codons_for_aa_returns_expected_codon_counts():
    # Classic genetic-code redundancies.
    assert len(codons_for_aa("M")) == 1  # methionine
    assert len(codons_for_aa("W")) == 1  # tryptophan
    assert len(codons_for_aa("L")) == 6  # leucine
    assert len(codons_for_aa("S")) == 6  # serine
    assert len(codons_for_aa("A")) == 4
    assert len(codons_for_aa("*")) == 3  # stop


def test_codons_are_cytosine_first_ordered():
    # For any aa with multiple codons, codons starting with 'C' should sort
    # before codons starting with other bases in the list.
    for aa in "LVPR":
        codons = codons_for_aa(aa)
        c_first = [c for c in codons if c.startswith("C")]
        others = [c for c in codons if not c.startswith("C")]
        assert codons == c_first + others, (
            f"{aa} codons not cytosine-first: {codons}"
        )


def test_codon_frequency_between_zero_and_one():
    # Every codon has a frequency in the human table.
    for aa, codons in CODONS_BY_AA.items():
        for codon in codons:
            f = codon_frequency(aa, codon)
            assert 0.0 <= f <= 1.0, f"frequency {f} out of range for {codon}"


def test_codon_frequency_sums_to_one_per_aa():
    # The h_sapiens_9606 table rounds to 2 decimals, so per-aa sums are 1.0
    # within rounding (typically ±0.02).
    for aa, codons in CODONS_BY_AA.items():
        total = sum(codon_frequency(aa, c) for c in codons)
        assert total == pytest.approx(1.0, abs=0.05), (
            f"codon frequencies for {aa} sum to {total}, not close to 1.0"
        )


def test_human_codon_frequencies_is_mapping_of_codons():
    # Uppercase 3-letter DNA codons only.
    assert all(len(c) == 3 and set(c) <= set("ACGT") for c in HUMAN_CODON_FREQUENCIES)
    # And 64 codons total.
    assert len(HUMAN_CODON_FREQUENCIES) == 64


def test_unknown_amino_acid_raises():
    with pytest.raises(KeyError):
        codons_for_aa("Z")
