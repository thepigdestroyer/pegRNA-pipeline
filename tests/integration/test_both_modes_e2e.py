"""End-to-end: run both TJ-PE and PA workflows on the same gene.

This test exists to catch any hidden coupling between the two workflows
and to document the contract that both produce self-consistent results
on the shared fixture data.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Seq import Seq

from pegrna_pipeline.hpc.data_types import GeneRecord
from pegrna_pipeline.hpc.pa_workflow import PAResult, run_pa_workflow
from pegrna_pipeline.hpc.tjpe_workflow import (
    STANDARD_SGRNA_SCAFFOLD,
    TJPEResult,
    run_tjpe_workflow,
)

FIXTURES = Path(__file__).parent.parent / "fixtures" / "pridict"
PA_UP = FIXTURES / "pa_upstream"
PA_DOWN = FIXTURES / "pa_downstream"


UPSTREAM_FLAP = "GCAGCAGCAGCAGCAGCAGCAGCAGCAGCA"
DOWNSTREAM_FLAP = "GCCGCCGCCGCCGCCGCCGCCGCCGCCGCC"
BULK = "AAACCCGGGTTTAAACCCGGGTTT"
TJPE_INSERT = UPSTREAM_FLAP + BULK + DOWNSTREAM_FLAP  # 84 bp


def _gene() -> GeneRecord:
    upstream = "ACGT" * 250
    downstream = "TGCA" * 250
    return GeneRecord(
        gene_name="CTCF",
        refseq_accession="NM_006565.4",
        mrna_sequence=upstream + "TAA" + downstream,
        cds_start=0,
        cds_end=1003,
        stop_codon_start=1000,
        stop_codon="TAA",
        upstream_flank=upstream[-100:],
        downstream_flank=downstream[:100],
        flanking_bp=100,
    )


@pytest.fixture(scope="module")
def tjpe_result() -> TJPEResult:
    return run_tjpe_workflow(
        gene=_gene(),
        insert_full_dna=TJPE_INSERT,
        variant_fixtures_dir=FIXTURES,
    )


@pytest.fixture(scope="module")
def pa_result() -> PAResult:
    return run_pa_workflow(
        gene=_gene(),
        upstream_fixtures_dir=PA_UP,
        downstream_fixtures_dir=PA_DOWN,
        upstream_flap_dna=UPSTREAM_FLAP,
        downstream_flap_dna=DOWNSTREAM_FLAP,
        bulk_dna=BULK,
    )


# -------------------------------------------------------------- TJ-PE checks


def test_tjpe_result_internally_consistent(tjpe_result: TJPEResult):
    assert STANDARD_SGRNA_SCAFFOLD in tjpe_result.tjpe_pegRNA_sequence
    # Re-derive rc_pbs2 independently.
    mrna = _gene().mrna_sequence
    nick = tjpe_result.ngrna_nick_in_mrna
    expected = str(Seq(mrna[nick - tjpe_result.rc_pbs2_length : nick]).reverse_complement())
    assert tjpe_result.rc_pbs2 == expected


# -------------------------------------------------------------- PA checks


def test_pa_donor_contains_flaps(pa_result: PAResult):
    assert pa_result.linear_donor_dna.startswith(UPSTREAM_FLAP)
    assert pa_result.linear_donor_dna.endswith(DOWNSTREAM_FLAP)


def test_pa_pair_opposite_strands(pa_result: PAResult):
    assert pa_result.winning_pair is not None
    assert (
        pa_result.winning_pair.upstream.target_strand
        != pa_result.winning_pair.downstream.target_strand
    )


# ---------------------------------------------------------- cross-workflow


def test_both_workflows_independent(tjpe_result: TJPEResult, pa_result: PAResult):
    """TJ-PE and PA operate on the same gene but produce different artefacts."""
    # Different artefact types.
    assert isinstance(tjpe_result, TJPEResult)
    assert isinstance(pa_result, PAResult)
    # TJ-PE assembles a single pegRNA; PA assembles a donor DNA.
    assert STANDARD_SGRNA_SCAFFOLD in tjpe_result.tjpe_pegRNA_sequence
    assert STANDARD_SGRNA_SCAFFOLD not in pa_result.linear_donor_dna  # donor has no scaffold
    # Both reference the same gene's stop codon arithmetic.
    assert tjpe_result.edit_position_in_mrna == _gene().stop_codon_start


def test_both_workflows_see_same_insert_payload(tjpe_result: TJPEResult, pa_result: PAResult):
    """TJ-PE's insert and PA's linear_donor should both carry the same 84 bp payload."""
    assert tjpe_result.full_insert_dna == TJPE_INSERT
    assert pa_result.linear_donor_dna == UPSTREAM_FLAP + BULK + DOWNSTREAM_FLAP
    assert tjpe_result.full_insert_dna == pa_result.linear_donor_dna
