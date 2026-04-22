"""T13: a winning ngRNA with negative ``distance_to_edit``.

A negative distance means the nick is 5' of the edit. The workflow must:

* Accept it (|d| in [40,100] is what matters for filter).
* Compute ``ngrna_nick_in_mrna = edit_position + distance`` (so the nick
  falls UPSTREAM in mRNA coordinates — a smaller index).
* Extract PBS2 correctly from that upstream region — not crash.
"""

from __future__ import annotations

import pytest
from Bio.Seq import Seq

from pegrna_pipeline.hpc.data_types import (
    GeneRecord,
    NgRNADesign,
    PegRNADesign,
    VariantScore,
)
from pegrna_pipeline.hpc.tjpe_workflow import run_tjpe_workflow


def _gene() -> GeneRecord:
    upstream = "ACGT" * 250  # 1000 bp
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
        sequence_name="v1_pegRNA_rv",
        spacer="CGTCGACGCTTAGGAACACG",
        pbs="TTCGTGTTCCTAAGCG",
        rt="CAGCGTGCAATCAGCGTGCAATCAGCGTGC",
        pbs_length=12,
        rt_length=30,
        target_strand="Rv",
        editing_position=50,
        pe2_hek293t=0.6,
        pe2_k562=0.5,
        pegrna_full="A" * 140,
    )


def test_negative_distance_nick_extracts_correct_rc_pbs2():
    ng = NgRNADesign(
        spacer="AGTGCTACCTGGTACGGATC",
        distance_to_edit=-50,  # nick at mRNA 950
        target_strand="Fw",  # opposite of peg (Rv)
        deepcas9_score=0.8,
        pe3b=False,
        pam_disrupt=False,
    )
    variant = VariantScore(
        sequence_name="v1",
        codon_variant_dna="ACG",
        pegRNA_designs=[_peg()],
        ngRNA_candidates=[ng],
    )
    result = run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna="ATGCATGCATGCATGC",
        variant_scores=[variant],
    )
    # Nick arithmetic must preserve sign.
    assert result.ngrna_distance_to_edit == -50
    assert result.ngrna_nick_in_mrna == 1000 - 50  # == 950

    # rc_pbs2 == revcomp( mrna[950 - 17 : 950] )
    expected_window = _gene().mrna_sequence[950 - 17 : 950]
    expected_rc = str(Seq(expected_window).reverse_complement())
    assert result.rc_pbs2 == expected_rc
    assert result.rc_pbs2_length == 17


def test_negative_distance_large_magnitude_still_valid():
    """|-95| = 95, within [40, 100] — must be accepted."""
    ng = NgRNADesign(
        spacer="AGTGCTACCTGGTACGGATC",
        distance_to_edit=-95,
        target_strand="Fw",
        deepcas9_score=0.8,
        pe3b=False,
        pam_disrupt=False,
    )
    variant = VariantScore(
        sequence_name="v1",
        codon_variant_dna="ACG",
        pegRNA_designs=[_peg()],
        ngRNA_candidates=[ng],
    )
    result = run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna="ATGCATGCATGCATGC",
        variant_scores=[variant],
    )
    assert result.ngrna_distance_to_edit == -95
    assert result.ngrna_nick_in_mrna == 905
    assert result.rc_pbs2  # non-empty
