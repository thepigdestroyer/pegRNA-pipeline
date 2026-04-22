"""Tests for the shared dataclasses in :mod:`pegrna_pipeline.hpc.data_types`.

GeneRecord is exercised implicitly through sequence_fetcher / input_formatter;
this file covers PegRNADesign, NgRNADesign, and VariantScore (the composite
that also carries the ``best_pegRNA`` / ``top_score`` helpers).
"""

from __future__ import annotations

import dataclasses

import pytest

from pegrna_pipeline.hpc.data_types import (
    GeneRecord,
    NgRNADesign,
    PegRNADesign,
    VariantScore,
)


def _peg(name: str = "v1", score: float = 0.5) -> PegRNADesign:
    return PegRNADesign(
        sequence_name=name,
        spacer="A" * 20,
        pbs="C" * 13,
        rt="G" * 30,
        pbs_length=13,
        rt_length=30,
        target_strand="Fw",
        editing_position=50,
        pe2_hek293t=score,
        pe2_k562=score * 0.8,
        pegrna_full="A" * 140,
    )


def _ng(distance: int = 50, score: float = 0.7) -> NgRNADesign:
    return NgRNADesign(
        spacer="T" * 20,
        distance_to_edit=distance,
        target_strand="Rv",
        deepcas9_score=score,
        pe3b=False,
        pam_disrupt=False,
    )


# -------------------------------------------------------------- frozenness


def test_pegrna_design_is_frozen():
    p = _peg()
    with pytest.raises(dataclasses.FrozenInstanceError):
        p.pe2_hek293t = 1.0  # type: ignore[misc]


def test_ngrna_design_is_frozen():
    n = _ng()
    with pytest.raises(dataclasses.FrozenInstanceError):
        n.deepcas9_score = 1.0  # type: ignore[misc]


def test_gene_record_is_frozen():
    g = GeneRecord(
        gene_name="X",
        refseq_accession="NM_X",
        mrna_sequence="ATGTAA",
        cds_start=0,
        cds_end=6,
        stop_codon_start=3,
        stop_codon="TAA",
        upstream_flank="",
        downstream_flank="",
        flanking_bp=0,
    )
    with pytest.raises(dataclasses.FrozenInstanceError):
        g.gene_name = "Y"  # type: ignore[misc]


def test_variant_score_is_not_frozen():
    """VariantScore is intentionally mutable so we can append designs as they arrive."""
    vs = VariantScore(
        sequence_name="v1",
        codon_variant_dna="ACG" * 10,
        pegRNA_designs=[],
        ngRNA_candidates=[],
    )
    vs.pegRNA_designs = [_peg()]  # must not raise
    assert len(vs.pegRNA_designs) == 1


# ------------------------------------------------------------- best_pegRNA


def test_best_pegrna_picks_max_by_pe2_hek293t():
    vs = VariantScore(
        sequence_name="v1",
        codon_variant_dna="ACG",
        pegRNA_designs=[_peg("a", 0.2), _peg("b", 0.7), _peg("c", 0.5)],
        ngRNA_candidates=[],
    )
    assert vs.best_pegRNA is not None
    assert vs.best_pegRNA.sequence_name == "b"


def test_best_pegrna_is_none_when_empty():
    vs = VariantScore(
        sequence_name="v1",
        codon_variant_dna="ACG",
        pegRNA_designs=[],
        ngRNA_candidates=[],
    )
    assert vs.best_pegRNA is None


def test_top_score_reflects_best_design():
    vs = VariantScore(
        sequence_name="v1",
        codon_variant_dna="ACG",
        pegRNA_designs=[_peg("a", 0.2), _peg("b", 0.9)],
        ngRNA_candidates=[],
    )
    assert vs.top_score == pytest.approx(0.9)


def test_top_score_zero_when_no_designs():
    vs = VariantScore(
        sequence_name="v1",
        codon_variant_dna="ACG",
        pegRNA_designs=[],
        ngRNA_candidates=[],
    )
    # Guarantee: top_score returns a real number (0.0) when there are no designs,
    # so callers can sort without guarding against None.
    assert vs.top_score == 0.0


# ------------------------------------------------------------ target_strand literal


def test_pegrna_target_strand_literal():
    # Not enforced at runtime by dataclass, but the field is typed as
    # ``Literal["Fw", "Rv"]``. Sanity-check that both values round-trip.
    for strand in ("Fw", "Rv"):
        p = dataclasses.replace(_peg(), target_strand=strand)
        assert p.target_strand == strand
