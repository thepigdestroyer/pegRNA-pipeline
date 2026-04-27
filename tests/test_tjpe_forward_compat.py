"""T13: forward-compat test exercising the ``variant_scores`` entry point.

Plan 3 will build :class:`VariantScore` objects in-memory (not from CSV
fixtures). This test proves that pathway — no fixtures directory, just
a hand-rolled list — produces the same shape of result as the
fixture-driven path.
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.data_types import (
    GeneRecord,
    NgRNADesign,
    PegRNADesign,
    VariantScore,
)
from pegrna_pipeline.hpc.tjpe_workflow import (
    STANDARD_SGRNA_SCAFFOLD,
    TJPEResult,
    run_tjpe_workflow,
)


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


def _peg(name: str, score: float) -> PegRNADesign:
    return PegRNADesign(
        sequence_name=name,
        spacer="GAATTCGTGTTCCTAAGCGT",
        pbs="CGCTTAGGAACAC",
        rt="GCACGCTGATTGCACGCTGATTGCACGCTG",
        pbs_length=13,
        rt_length=30,
        target_strand="Fw",
        editing_position=50,
        pe2_hek293t=score,
        pe2_k562=score * 0.8,
        pegrna_full="A" * 140,
    )


def _ng(distance: int = 60, score: float = 0.8) -> NgRNADesign:
    return NgRNADesign(
        spacer="AGTGCTACCTGGTACGGATC",
        distance_to_edit=distance,
        target_strand="Rv",  # opposite of Fw
        deepcas9_score=score,
        pe3b=False,
        pam_disrupt=False,
    )


def test_programmatic_variant_scores_produces_valid_result():
    v1 = VariantScore(
        sequence_name="prog_v1",
        codon_variant_dna="GCAGCAGCA",
        pegRNA_designs=[_peg("prog_v1_pegRNA1", 0.3)],
        ngRNA_candidates=[_ng()],
    )
    v2 = VariantScore(
        sequence_name="prog_v2",
        codon_variant_dna="GCCGCCGCC",
        pegRNA_designs=[_peg("prog_v2_pegRNA1", 0.8)],
        ngRNA_candidates=[_ng()],
    )

    result = run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna="ATGCATGCATGCATGC",
        variant_scores=[v1, v2],
    )

    assert isinstance(result, TJPEResult)
    # v2 scores higher, so it wins.
    assert result.winning_variant is v2
    # Default ranking is K562; pe2_k562 = score * 0.8 → 0.8 * 0.8 = 0.64
    assert result.pegrna_pe2_score == pytest.approx(0.64)
    assert result.pegrna_pe2_cell_line == "k562"
    assert STANDARD_SGRNA_SCAFFOLD in result.tjpe_pegRNA_sequence
    # Evaluated list preserves both variants.
    assert len(result.all_evaluated_variants) == 2


def test_programmatic_and_fixture_paths_are_mutually_exclusive(tmp_path):
    with pytest.raises(ValueError):
        run_tjpe_workflow(
            gene=_gene(),
            insert_full_dna="ATGC",
            variant_fixtures_dir=tmp_path,
            variant_scores=[],
        )
