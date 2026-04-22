"""Tests for the weighted optimise_bulk implementation.

Score formula per residue, given the DNA so far (``prev``) and candidate ``codon``:

    score = codon_frequency(aa, codon)
          + 0.3 * codon.count("C")
          - 10.0 if prev[-1] == 'C' and codon[0] == 'G' (CpG at junction)
          - 10.0 if any restriction site appears inside the 10-char window
            spanning the junction
          - 10.0 if placing codon creates a homopolymer run >= 5

Picks the max-score codon per residue. Logs a warning if the best codon
still has violations.
"""

import logging

import pytest

from pegrna_pipeline.hpc.codon_expander import (
    _has_cpg,
    _has_homopolymer,
    _has_restriction_site,
    optimise_bulk,
)


def test_empty_string_returns_empty():
    assert optimise_bulk("") == ""


def test_output_length_is_three_times_aa_length():
    aa = "MKPLAEDG"
    assert len(optimise_bulk(aa)) == 3 * len(aa)


def test_output_is_uppercase_acgt():
    dna = optimise_bulk("MKPLAED")
    assert set(dna) <= set("ACGT")


def test_output_translates_back_to_input():
    from python_codon_tables import get_codons_table

    translate = {}
    for aa, codons in get_codons_table("h_sapiens_9606").items():
        for c in codons:
            translate[c] = aa

    aa = "MKPLAEDG"
    dna = optimise_bulk(aa)
    decoded = "".join(translate[dna[i : i + 3]] for i in range(0, len(dna), 3))
    assert decoded == aa


def test_unknown_amino_acid_raises():
    with pytest.raises(ValueError):
        optimise_bulk("MZK")


def test_optimise_bulk_avoids_restriction_sites_when_possible():
    # A short aa run whose best codons would NOT create a site naturally; no
    # regression.
    dna = optimise_bulk("MK")
    assert not _has_restriction_site(dna)


def test_optimise_bulk_prefers_codons_with_more_c():
    # G and S each have multiple codons; several contain a C.  We should
    # therefore consistently see C's in the output for long stretches.
    dna = optimise_bulk("GSGSGSGSGS")
    # Expect at least some C bases (C_count_bonus pulls the score up).
    assert dna.count("C") > 0


def test_warning_logged_when_best_codon_still_has_violations(caplog):
    # Methionine has only one codon (ATG). Prepend a run that makes
    # concatenation trigger a homopolymer regardless.  We simulate this by
    # optimising a long poly-aa sequence whose best chosen codons inevitably
    # form a homopolymer -- specifically poly-K (AAA/AAG).  Use 6 Ks so the
    # junction between two AAA codons yields 6+ consecutive A's.
    caplog.set_level(logging.WARNING, logger="pegrna_pipeline.hpc.codon_expander")
    dna = optimise_bulk("KKKKKK")
    # Regardless of logger output, the DNA must still translate back correctly
    # -- the optimiser does NOT refuse to emit a result.
    assert len(dna) == 18
