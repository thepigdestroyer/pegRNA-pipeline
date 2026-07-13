"""Core unit tests for :func:`run_tjpe_workflow`.

Edge-case suites live in ``test_tjpe_edge_*``; here we cover the happy
path and the mutually-exclusive inputs contract.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pegrna_pipeline.hpc.data_types import GeneRecord
from pegrna_pipeline.hpc.tjpe_workflow import (
    STANDARD_SGRNA_SCAFFOLD,
    TJPEResult,
    run_tjpe_workflow,
)

FIXTURES = Path(__file__).parent / "fixtures" / "pridict"


# A synthetic gene with enough flanking to accommodate ngRNA distances up to 100
# bp past the stop codon. Stop codon at position 1000; 200 bp upstream, 200 bp
# downstream. mRNA is deterministic so PBS2 extraction is reproducible.
def _gene() -> GeneRecord:
    upstream = "ACGT" * 250  # 1000 bp upstream
    downstream = "TGCA" * 250  # 1000 bp downstream
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


INSERT_DNA = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"  # 40 bp fake tag


def test_run_tjpe_workflow_end_to_end_with_fixtures():
    result = run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna=INSERT_DNA,
        variant_fixtures_dir=FIXTURES,
    )
    assert isinstance(result, TJPEResult)
    assert result.full_insert_dna == INSERT_DNA
    # Scaffold must appear in the assembled pegRNA.
    assert STANDARD_SGRNA_SCAFFOLD in result.tjpe_pegRNA_sequence
    # Default ranking is K562; pegRNA1 still wins (K562=0.511 vs 0.388, 0.301 in fixture).
    assert result.pegrna_pe2_score == pytest.approx(0.511)
    assert result.pegrna_pe2_cell_line == "k562"
    # The winning pegRNA is pegRNA1 (Fw) — so ngRNA must be Rv.
    assert result.ngrna_spacer
    # Arithmetic per v1.1: edit_position_in_mrna == stop_codon_start.
    assert result.edit_position_in_mrna == 1000
    # ngrna_nick_in_mrna = edit_position_in_mrna + distance_to_edit.
    assert result.ngrna_nick_in_mrna == 1000 + result.ngrna_distance_to_edit
    # rc_pbs2 length matches pbs2_length (default 17).
    assert result.rc_pbs2_length == 17
    assert len(result.rc_pbs2) == 17
    # Full pegRNA length is the sum of components in the authoritative v2 layout:
    # spacer + scaffold + rcPBS2 + RTT(full_insert) + PBS1 + PolyU(UUU).
    expected_len = (
        len(result.spacer)
        + len(STANDARD_SGRNA_SCAFFOLD)
        + len(result.rc_pbs2)
        + len(result.full_insert_dna)
        + len(result.pbs1)
        + 3  # PolyU tail (2-3 U); assemble_tj_pegrna_v2 default is "UUU"
    )
    assert len(result.tjpe_pegRNA_sequence) == expected_len


def test_run_tjpe_requires_exactly_one_input_source():
    """Exactly one of variant_fixtures_dir / variant_scores must be set."""
    with pytest.raises(ValueError):
        run_tjpe_workflow(gene=_gene(), insert_full_dna=INSERT_DNA)
    with pytest.raises(ValueError):
        run_tjpe_workflow(
            gene=_gene(),
            insert_full_dna=INSERT_DNA,
            variant_fixtures_dir=FIXTURES,
            variant_scores=[],
        )


def test_run_tjpe_records_all_variants():
    result = run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna=INSERT_DNA,
        variant_fixtures_dir=FIXTURES,
    )
    # fixture has one variant prefix
    assert len(result.all_evaluated_variants) == 1
    assert result.all_evaluated_variants[0].sequence_name.startswith("ctcf_tjpe_v1")


def test_run_tjpe_winning_variant_has_best_score():
    result = run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna=INSERT_DNA,
        variant_fixtures_dir=FIXTURES,
    )
    assert result.winning_variant is not None
    # Default is K562 ranking → pegRNA1's K562 score = 0.511.
    assert result.winning_variant.top_score == pytest.approx(0.511)
