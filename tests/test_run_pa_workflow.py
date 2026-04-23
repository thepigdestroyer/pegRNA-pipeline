"""Core tests for :func:`run_pa_workflow`.

Happy-path runs use the fixture pair. Edge cases (no compatible pairs,
mutually-exclusive inputs) live in dedicated files.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pegrna_pipeline.hpc.data_types import GeneRecord
from pegrna_pipeline.hpc.pa_workflow import PAResult, run_pa_workflow

FIXTURES = Path(__file__).parent / "fixtures" / "pridict"
UPSTREAM_DIR = FIXTURES / "pa_upstream"
DOWNSTREAM_DIR = FIXTURES / "pa_downstream"


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


UPSTREAM_FLAP = "GCAGCAGCAGCAGCAGCAGCAGCAGCAGCA"  # 30 bp = 10 aa
DOWNSTREAM_FLAP = "GCCGCCGCCGCCGCCGCCGCCGCCGCCGCC"  # 30 bp = 10 aa
BULK = "AAACCCGGGTTTAAACCCGGGTTT"  # 24 bp


def test_run_pa_happy_path():
    result = run_pa_workflow(
        gene=_gene(),
        upstream_fixtures_dir=UPSTREAM_DIR,
        downstream_fixtures_dir=DOWNSTREAM_DIR,
        upstream_flap_dna=UPSTREAM_FLAP,
        downstream_flap_dna=DOWNSTREAM_FLAP,
        bulk_dna=BULK,
    )
    assert isinstance(result, PAResult)
    # Donor is assembled from flaps + bulk.
    assert result.linear_donor_dna == UPSTREAM_FLAP + BULK + DOWNSTREAM_FLAP
    assert result.linear_donor_dna.startswith(UPSTREAM_FLAP)
    assert result.linear_donor_dna.endswith(DOWNSTREAM_FLAP)
    # Assembly must be internally consistent.
    assert result.upstream_flap_dna == UPSTREAM_FLAP
    assert result.downstream_flap_dna == DOWNSTREAM_FLAP
    assert result.bulk_dna == BULK

    # Winning pair survived.
    assert result.winning_pair is not None
    assert (
        result.winning_pair.upstream.target_strand != result.winning_pair.downstream.target_strand
    )

    # Top pairs returned + counts consistent.
    assert result.pairs_after_compat_filter == len(result.top_pairs)
    assert result.total_pairs_evaluated >= result.pairs_after_compat_filter


def test_run_pa_uses_donor_dna_if_provided():
    """If the caller supplies a pre-built donor_dna, skip assembly."""
    pre_built = "A" * 50
    result = run_pa_workflow(
        gene=_gene(),
        upstream_fixtures_dir=UPSTREAM_DIR,
        downstream_fixtures_dir=DOWNSTREAM_DIR,
        donor_dna=pre_built,
    )
    assert result.linear_donor_dna == pre_built


def test_run_pa_requires_exactly_one_upstream_source():
    """Caller must supply exactly one of upstream_fixtures_dir /
    upstream_variant_scores, never both or neither."""
    with pytest.raises(ValueError, match="upstream"):
        run_pa_workflow(
            gene=_gene(),
            donor_dna="A" * 50,
            # neither upstream source provided — downstream OK
            downstream_fixtures_dir=DOWNSTREAM_DIR,
        )


def test_run_pa_requires_exactly_one_downstream_source():
    with pytest.raises(ValueError, match="downstream"):
        run_pa_workflow(
            gene=_gene(),
            donor_dna="A" * 50,
            upstream_fixtures_dir=UPSTREAM_DIR,
            # neither downstream source
        )


def test_run_pa_rejects_both_fixture_and_scores_per_flap():
    with pytest.raises(ValueError):
        run_pa_workflow(
            gene=_gene(),
            donor_dna="A" * 50,
            upstream_fixtures_dir=UPSTREAM_DIR,
            upstream_variant_scores=[],
            downstream_fixtures_dir=DOWNSTREAM_DIR,
        )


def test_run_pa_counts_are_coherent():
    result = run_pa_workflow(
        gene=_gene(),
        upstream_fixtures_dir=UPSTREAM_DIR,
        downstream_fixtures_dir=DOWNSTREAM_DIR,
        upstream_flap_dna=UPSTREAM_FLAP,
        downstream_flap_dna=DOWNSTREAM_FLAP,
        bulk_dna=BULK,
    )
    # Fixture: 2 upstream x 2 downstream = 4 combos pre-filter.
    assert result.total_pairs_evaluated == 4
    # Some may be filtered out (same strand, distance) — but count >= 1.
    assert result.pairs_after_compat_filter >= 1
