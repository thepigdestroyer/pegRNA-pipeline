"""Tests for :func:`pegrna_pipeline.hpc.input_formatter.format_bracket_insertion`.

The bracket-insertion string fed to PRIDICT 2.0 is
``{upstream_flank}(+{insertion}){downstream_flank}`` — the parentheses +
leading ``+`` are PRIDICT's marker for an *insertion* edit between two
existing bases.
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.data_types import GeneRecord
from pegrna_pipeline.hpc.input_formatter import format_bracket_insertion


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
