"""Tests for :func:`assemble_tj_pegrna` and :data:`STANDARD_SGRNA_SCAFFOLD`.

A TJ-PE pegRNA is laid out 5' -> 3':

    spacer  +  scaffold  +  insert_dna  +  pbs1  +  rc_pbs2

The scaffold is a constant so tests can pin its identity.
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.tjpe_workflow import (
    STANDARD_SGRNA_SCAFFOLD,
    assemble_tj_pegrna,
)


def test_scaffold_constant_matches_plan_spec():
    expected = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
    assert STANDARD_SGRNA_SCAFFOLD == expected
    assert len(STANDARD_SGRNA_SCAFFOLD) == 76
    assert set(STANDARD_SGRNA_SCAFFOLD) <= set("ACGT")


def test_assembly_concatenates_in_order():
    out = assemble_tj_pegrna(
        spacer="A" * 20,
        insert_dna="C" * 30,
        pbs1="G" * 13,
        rc_pbs2="T" * 17,
    )
    assert out == ("A" * 20) + STANDARD_SGRNA_SCAFFOLD + ("C" * 30) + ("G" * 13) + ("T" * 17)
    # Sanity: length = spacer + scaffold + insert + pbs1 + rc_pbs2.
    assert len(out) == 20 + len(STANDARD_SGRNA_SCAFFOLD) + 30 + 13 + 17


def test_custom_scaffold_supported():
    custom = "TTTTGGGGAAAACCCC"
    out = assemble_tj_pegrna(
        spacer="AAAA",
        insert_dna="CC",
        pbs1="GG",
        rc_pbs2="TT",
        scaffold=custom,
    )
    assert custom in out
    assert out.startswith("AAAA")
    assert out.endswith("TT")


def test_assembly_rejects_non_acgt_in_any_component():
    with pytest.raises(ValueError):
        assemble_tj_pegrna(spacer="AACG", insert_dna="ACGN", pbs1="ACG", rc_pbs2="ACG")
    with pytest.raises(ValueError):
        assemble_tj_pegrna(spacer="AAcG", insert_dna="ACGT", pbs1="ACG", rc_pbs2="ACG")
    with pytest.raises(ValueError):
        assemble_tj_pegrna(spacer="ACGT", insert_dna="ACGT", pbs1="", rc_pbs2="ACG")


def test_assembly_uppercase_only():
    with pytest.raises(ValueError):
        assemble_tj_pegrna(spacer="acgt", insert_dna="ACGT", pbs1="ACGT", rc_pbs2="ACGT")
