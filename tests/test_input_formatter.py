"""Tests for :func:`pegrna_pipeline.hpc.input_formatter.format_bracket_insertion`.

The bracket-insertion string fed to PRIDICT 2.0 is
``{upstream_flank}(+{insertion}){downstream_flank}`` — the parentheses +
leading ``+`` are PRIDICT's marker for an *insertion* edit between two
existing bases.
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.data_types import GeneRecord
from pegrna_pipeline.hpc.input_formatter import (
    format_bracket_insertion,
    format_tjpe_insertion_editseq,
)


def _make_gene(upstream: str = "A" * 100, downstream: str = "G" * 100) -> GeneRecord:
    return GeneRecord(
        gene_name="X",
        refseq_accession="NM_X",
        mrna_sequence=upstream + "TAA" + downstream,
        cds_start=0,
        cds_end=len(upstream) + 3,
        stop_codon_start=len(upstream),
        stop_codon="TAA",
        upstream_flank=upstream,
        downstream_flank=downstream,
        flanking_bp=len(upstream),
    )


def test_format_wraps_insertion_in_parens_with_plus():
    gene = _make_gene(upstream="ACGT", downstream="TTTT")
    assert format_bracket_insertion(gene, "GGG") == "ACGT(+GGG)TTTT"


def test_format_uses_gene_flanks_verbatim():
    gene = _make_gene()
    out = format_bracket_insertion(gene, "CCC")
    assert out.startswith(gene.upstream_flank + "(+")
    assert out.endswith(")" + gene.downstream_flank)


def test_format_empty_insertion_raises():
    gene = _make_gene(upstream="A", downstream="G")
    with pytest.raises(ValueError):
        format_bracket_insertion(gene, "")


def test_format_rejects_non_acgt():
    gene = _make_gene(upstream="A", downstream="G")
    with pytest.raises(ValueError):
        format_bracket_insertion(gene, "ACGU")  # U invalid
    with pytest.raises(ValueError):
        format_bracket_insertion(gene, "acgt")  # lowercase invalid
    with pytest.raises(ValueError):
        format_bracket_insertion(gene, "AC GT")  # spaces invalid
    with pytest.raises(ValueError):
        format_bracket_insertion(gene, "ACNT")  # N invalid


def test_format_accepts_single_base():
    gene = _make_gene(upstream="A", downstream="G")
    assert format_bracket_insertion(gene, "C") == "A(+C)G"


def test_format_accepts_long_insertion():
    gene = _make_gene()
    big = "ACGT" * 50  # 200 bp
    out = format_bracket_insertion(gene, big)
    assert f"(+{big})" in out


# --------------------------------------------------------------------------
# format_tjpe_insertion_editseq: tj-PE insertion editseq from a genome window.
#
# Bracket convention for an INSERTION is ``FLANK(+INSERT)FLANK`` (reference
# empty), CONFIRMED against a live PRIDICT2 run on the DSI box (2026-07-10):
# ``(+SEQ)`` yields a valid pegRNA, ``(/SEQ)`` errors. This matches the legacy
# format_bracket_insertion marker.
# --------------------------------------------------------------------------


def test_tjpe_insertion_places_bracket_at_site_s():
    genome = "A" * 10 + "C" * 10  # len 20
    out = format_tjpe_insertion_editseq(genome, 10, "GGG", min_flank=5)
    assert out == "A" * 10 + "(+GGG)" + "C" * 10


def test_tjpe_insertion_uses_plus_bracket():
    genome = "A" * 10 + "C" * 10
    out = format_tjpe_insertion_editseq(genome, 10, "GGG", min_flank=5)
    assert "(+GGG)" in out
    assert "(/GGG)" not in out


def test_tjpe_insertion_site_s_at_start_of_right_flank():
    # site_s slices genome into [:site_s] and [site_s:] verbatim.
    genome = "ACGT" * 30  # len 120
    out = format_tjpe_insertion_editseq(genome, 60, "TT", min_flank=10)
    assert out == genome[:60] + "(+TT)" + genome[60:]


def test_tjpe_insertion_uppercases_genome_and_insert():
    genome = ("a" * 10) + ("c" * 10)
    out = format_tjpe_insertion_editseq(genome, 10, "ggg", min_flank=5)
    assert out == "A" * 10 + "(+GGG)" + "C" * 10


def test_tjpe_insertion_empty_insert_raises():
    genome = "A" * 200
    with pytest.raises(ValueError):
        format_tjpe_insertion_editseq(genome, 100, "")


def test_tjpe_insertion_rejects_non_acgt_insert():
    genome = "A" * 200
    with pytest.raises(ValueError):
        format_tjpe_insertion_editseq(genome, 100, "ACGU")  # U invalid
    with pytest.raises(ValueError):
        format_tjpe_insertion_editseq(genome, 100, "ACNT")  # N invalid
    with pytest.raises(ValueError):
        format_tjpe_insertion_editseq(genome, 100, "AC GT")  # space invalid
    with pytest.raises(ValueError):
        format_tjpe_insertion_editseq(genome, 100, "AC1T")  # digit invalid


def test_tjpe_insertion_short_left_flank_raises():
    genome = "A" * 100
    # site_s=3 → left flank 3 < min_flank 5.
    with pytest.raises(ValueError):
        format_tjpe_insertion_editseq(genome, 3, "GGG", min_flank=5)


def test_tjpe_insertion_short_right_flank_raises():
    genome = "A" * 100
    # site_s=98 → right flank 2 < min_flank 5.
    with pytest.raises(ValueError):
        format_tjpe_insertion_editseq(genome, 98, "GGG", min_flank=5)


def test_tjpe_insertion_default_min_flank_is_100():
    # A 150 bp window centred cannot satisfy 100 bp each side by default.
    genome = "A" * 150
    with pytest.raises(ValueError):
        format_tjpe_insertion_editseq(genome, 75, "GGG")
    # 300 bp with site_s=150 gives exactly 150 bp each side → OK.
    ok = format_tjpe_insertion_editseq("A" * 300, 150, "GGG")
    assert ok == "A" * 150 + "(+GGG)" + "A" * 150


def test_tjpe_insertion_exact_min_flank_boundary_ok():
    genome = "A" * 10  # len 10
    out = format_tjpe_insertion_editseq(genome, 5, "T", min_flank=5)
    assert out == "A" * 5 + "(+T)" + "A" * 5


def test_tjpe_insertion_site_s_out_of_range_raises():
    genome = "A" * 200
    with pytest.raises(ValueError):
        format_tjpe_insertion_editseq(genome, -1, "GGG", min_flank=0)
    with pytest.raises(ValueError):
        format_tjpe_insertion_editseq(genome, len(genome) + 1, "GGG", min_flank=0)
