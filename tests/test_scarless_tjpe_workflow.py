"""Tests for scarless tjPE workflow wiring (assembly, distance, ranking, gates)."""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.data_types import GeneRecord, NgRNADesign
from pegrna_pipeline.hpc.tjpe_workflow import (
    ScarlessNoValidDesign,
    assemble_tj_pegrna,
    combine_two_end_scores,
    require_b_at_s,
    run_tjpe_workflow,
    select_tjpe_ngrna,
)


def _gene(seq="ACGT" * 40, stop=100):
    return GeneRecord(
        gene_name="TESTGENE",
        refseq_accession="NM_000000.0",
        mrna_sequence=seq,
        cds_start=0,
        cds_end=stop + 3,
        stop_codon_start=stop,
        stop_codon="TAA",
        upstream_flank="",
        downstream_flank="",
        flanking_bp=0,
    )


# ------------------------------------------------------------------ Task 5
def test_assemble_inserts_restored_segment_between_insert_and_pbs1():
    out = assemble_tj_pegrna(
        spacer="ACGTACGTACGTACGTACGT",
        insert_dna="ATGAAA",
        restored_segment="TACG",
        pbs1="CCCCGGGG",
        rc_pbs2="TTTAAAA",
    )
    assert "ATGAAA" + "TACG" + "CCCCGGGG" + "TTTAAAA" in out


def test_assemble_restored_empty_is_backwards_compatible():
    out = assemble_tj_pegrna(
        spacer="ACGTACGTACGTACGTACGT",
        insert_dna="ATGAAA",
        restored_segment="",
        pbs1="CCCCGGGG",
        rc_pbs2="TTTAAAA",
    )
    assert "ATGAAA" + "CCCCGGGG" + "TTTAAAA" in out


# ------------------------------------------------------------------ Task 6
def _ng(dist, strand="Rv", spacer="ACGTACGTACGTACGTACGT", score=50.0):
    return NgRNADesign(
        spacer=spacer,
        distance_to_edit=dist,
        target_strand=strand,
        deepcas9_score=score,
        pe3b=False,
        pam_disrupt=False,
    )


def test_select_ngrna_40bp_window():
    cands = [_ng(40), _ng(70), _ng(95)]
    chosen = select_tjpe_ngrna(cands, min_distance=36, max_distance=44, pegrna_strand="Fw")
    assert chosen is not None and abs(chosen.distance_to_edit) == 40


def test_select_ngrna_excludes_outside_window():
    cands = [_ng(70), _ng(95)]
    chosen = select_tjpe_ngrna(cands, min_distance=36, max_distance=44, pegrna_strand="Fw")
    assert chosen is None


def test_reverse_strand_requires_b_at_s():
    assert require_b_at_s(distance_to_edit=0, pegrna_strand="Rv") is True
    assert require_b_at_s(distance_to_edit=5, pegrna_strand="Rv") is False
    assert require_b_at_s(distance_to_edit=40, pegrna_strand="Fw") is True


def test_require_b_at_s_bad_strand():
    with pytest.raises(ValueError, match="pegrna_strand"):
        require_b_at_s(distance_to_edit=0, pegrna_strand="X")


# ------------------------------------------------------------------ Task 7
def test_combine_uses_minimum_bottleneck():
    assert combine_two_end_scores(insert_score=80.0, pbs2_score=30.0) == 30.0
    assert combine_two_end_scores(insert_score=10.0, pbs2_score=90.0) == 10.0


# ------------------------------------------------------------------ Task 8
def test_scarless_no_variants_fails_loudly():
    # No variants → scarless must raise, never return a silent empty design (R4).
    with pytest.raises(ScarlessNoValidDesign, match="(?i)no valid scarless design"):
        run_tjpe_workflow(
            gene=_gene(),
            insert_full_dna="ATGAAA",
            variant_scores=[],
            scarless=True,
        )


def test_nonscarless_no_variants_returns_empty_result():
    # Back-compat: non-scarless path still returns a zeroed result, no raise.
    res = run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna="ATGAAA",
        variant_scores=[],
        scarless=False,
    )
    assert res.tjpe_pegRNA_sequence == ""
    assert res.scarless is False


def test_scarless_success_path_zero_deletion_and_warnings():
    import random

    from pegrna_pipeline.hpc.data_types import PegRNADesign, VariantScore

    random.seed(1)
    seq = "".join(random.choice("ACGT") for _ in range(200))
    gene = _gene(seq=seq, stop=100)
    peg = PegRNADesign(
        "v1",
        "ACGTACGTACGTACGTACGT",
        "CCCAGTAGCTAGCTAGC",
        "RT",
        13,
        20,
        "Fw",
        0,
        60.0,
        55.0,
        "full",
    )
    ng = NgRNADesign(
        spacer="GCAGTCGATCGATCGATCGA",
        distance_to_edit=40,
        target_strand="Rv",
        deepcas9_score=70.0,
        pe3b=False,
        pam_disrupt=False,
    )
    vs = VariantScore("v1", "ATGAAACCCGGG", [peg], [ng])
    cargo = "ATGAAACCCGGGTTTAAA"
    res = run_tjpe_workflow(
        gene=gene,
        insert_full_dna=cargo,
        variant_scores=[vs],
        ngrna_distance_min=36,
        ngrna_distance_max=44,
        scarless=True,
    )
    assert res.scarless is True
    assert res.deletion_size == 0
    assert len(res.restored_segment) == 40  # the A→B span
    # INV-0: predicted sense allele == reference + cargo at S, nothing deleted.
    assert len(res.predicted_edited_allele) == len(seq) + len(cargo)
    assert res.predicted_edited_allele == seq[:100] + cargo + seq[100:]
    # Advisory warnings surfaced (microhomology + sub-70bp distance), never blocking.
    assert res.microhomology_flag is True
    assert any("microhomology" in w for w in res.warnings)
    assert any("70 bp" in w for w in res.warnings)
    # The restored segment is embedded in the assembled pegRNA before PBS1.
    assert res.restored_segment in res.tjpe_pegRNA_sequence
