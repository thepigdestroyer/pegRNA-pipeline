"""T14: full TJ-PE workflow on the bundled PRIDICT fixture pair.

This is *not* a simple re-run of the T10 unit test — it additionally
verifies internal consistency between the individual result fields
(pegRNA components reassemble back to the full sequence, nick position
computed correctly, RC-PBS2 verifiable from the mRNA).
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Seq import Seq

from pegrna_pipeline.hpc.data_types import GeneRecord
from pegrna_pipeline.hpc.tjpe_workflow import (
    STANDARD_SGRNA_SCAFFOLD,
    TJPEResult,
    run_tjpe_workflow,
)

FIXTURES = Path(__file__).parent.parent / "fixtures" / "pridict"


def _gene() -> GeneRecord:
    upstream = "ACGT" * 250
    downstream = "TGCA" * 250
    mrna = upstream + "TAA" + downstream
    return GeneRecord(
        gene_name="CTCF",
        refseq_accession="NM_006565.4",
        mrna_sequence=mrna,
        cds_start=0,
        cds_end=1003,
        stop_codon_start=1000,
        stop_codon="TAA",
        upstream_flank=upstream[-100:],
        downstream_flank=downstream[:100],
        flanking_bp=100,
    )


INSERT = "ATGGTGAGCAAGGGCGAGGAGCTG"  # 24 bp fake mini-tag


@pytest.fixture(scope="module")
def result() -> TJPEResult:
    return run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna=INSERT,
        variant_fixtures_dir=FIXTURES,
    )


def test_result_is_internally_consistent(result: TJPEResult):
    # pegRNA sequence is spacer + scaffold + insert + pbs1 + rc_pbs2.
    reassembled = (
        result.spacer
        + STANDARD_SGRNA_SCAFFOLD
        + result.full_insert_dna
        + result.pbs1
        + result.rc_pbs2
    )
    assert reassembled == result.tjpe_pegRNA_sequence


def test_rc_pbs2_matches_mrna_window(result: TJPEResult):
    """rc_pbs2 == revcomp(mrna[nick - len(rc_pbs2) : nick])."""
    mrna = _gene().mrna_sequence
    nick = result.ngrna_nick_in_mrna
    window = mrna[nick - result.rc_pbs2_length : nick]
    expected_rc = str(Seq(window).reverse_complement())
    assert result.rc_pbs2 == expected_rc


def test_nick_arithmetic(result: TJPEResult):
    assert result.edit_position_in_mrna == 1000  # stop_codon_start
    assert result.ngrna_nick_in_mrna == result.edit_position_in_mrna + result.ngrna_distance_to_edit


def test_winner_has_best_score(result: TJPEResult):
    evaluated = result.all_evaluated_variants
    top_scores = [v.top_score for v in evaluated]
    assert result.pegrna_pe2_score == max(top_scores)


def test_ngrna_strand_is_opposite_of_pegrna(result: TJPEResult):
    # Confirm the {Fw -> Rv} mapping actually held.
    winner = result.winning_variant
    assert winner is not None
    best_peg = winner.best_pegRNA
    assert best_peg is not None
    # We don't know exactly which ngRNA won, but we can reconstruct the strand
    # from the spacer sequence by scanning candidates.
    for ng in winner.ngRNA_candidates:
        if ng.spacer == result.ngrna_spacer:
            assert ng.target_strand != best_peg.target_strand
            break
    else:
        pytest.fail("winning ngRNA spacer not found in candidate list")
