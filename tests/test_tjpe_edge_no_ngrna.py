"""T12: the winning variant has pegRNAs but every ngRNA fails filters.

Expected: :class:`TJPEResult` with ``pegrna_pe2_hek293t`` populated from
the best pegRNA but ``ngrna_*`` fields empty and no assembled pegRNA
string. The caller can then fall back to the next-best variant.
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.data_types import (
    GeneRecord,
    NgRNADesign,
    PegRNADesign,
    VariantScore,
)
from pegrna_pipeline.hpc.tjpe_workflow import run_tjpe_workflow


def _gene() -> GeneRecord:
    upstream = "ACGT" * 250
    downstream = "TGCA" * 250
    mrna = upstream + "TAA" + downstream
    return GeneRecord(
        gene_name="FAKE",
        refseq_accession="NM_FAKE.1",
        mrna_sequence=mrna,
        cds_start=0,
        cds_end=1003,
        stop_codon_start=1000,
        stop_codon="TAA",
        upstream_flank=upstream[-100:],
        downstream_flank=downstream[:100],
        flanking_bp=100,
    )


def _peg() -> PegRNADesign:
    return PegRNADesign(
        sequence_name="v1_pegRNA1",
        spacer="GAATTCGTGTTCCTAAGCGT",
        pbs="CGCTTAGGAACAC",
        rt="GCACGCTGATTGCACGCTGATTGCACGCTG",
        pbs_length=13,
        rt_length=30,
        target_strand="Fw",
        editing_position=50,
        pe2_hek293t=0.7,
        pe2_k562=0.5,
        pegrna_full="A" * 140,
    )


def _bad_ng(distance: int, strand: str = "Rv") -> NgRNADesign:
    # Spacer has TTTT -> fails sanity filter.
    return NgRNADesign(
        spacer="ATTTTACCTGGTACGGATCA",
        distance_to_edit=distance,
        target_strand=strand,
        deepcas9_score=0.9,
        pe3b=False,
        pam_disrupt=False,
    )


def test_no_surviving_ngrna_returns_empty_assembly():
    variant = VariantScore(
        sequence_name="v1",
        codon_variant_dna="ACG" * 10,
        pegRNA_designs=[_peg()],
        # Every candidate is filtered by the sanity checks (TTTT in spacer).
        ngRNA_candidates=[_bad_ng(50), _bad_ng(70)],
    )
    result = run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna="ATGCATGCATGCATGC",
        variant_scores=[variant],
    )
    assert result.winning_variant is variant
    # Default ranking is K562; pe2_k562=0.5 in this fixture.
    assert result.pegrna_pe2_score == pytest.approx(0.5)
    assert result.pegrna_pe2_cell_line == "k562"
    # pegRNA fields populated from the winning design…
    assert result.spacer == "GAATTCGTGTTCCTAAGCGT"
    assert result.pbs1 == "CGCTTAGGAACAC"
    # …but no ngRNA picked, so no assembly and no rc_pbs2.
    assert result.tjpe_pegRNA_sequence == ""
    assert result.ngrna_spacer == ""
    assert result.ngrna_distance_to_edit == 0
    assert result.rc_pbs2_length == 0


def test_ngrna_outside_distance_window_is_dropped():
    variant = VariantScore(
        sequence_name="v1",
        codon_variant_dna="ACG",
        pegRNA_designs=[_peg()],
        ngRNA_candidates=[
            NgRNADesign(
                spacer="AGTGCTACCTGGTACGGATC",
                distance_to_edit=10,  # < 40
                target_strand="Rv",
                deepcas9_score=0.9,
                pe3b=False,
                pam_disrupt=False,
            )
        ],
    )
    result = run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna="ATGCATGCATGCATGC",
        variant_scores=[variant],
    )
    assert result.ngrna_spacer == ""
    assert result.tjpe_pegRNA_sequence == ""
