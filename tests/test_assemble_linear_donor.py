"""Tests for :func:`assemble_linear_donor`.

The linear donor template for Prime Assembly is just the concatenation
of (upstream_flap_dna, bulk_dna, downstream_flap_dna). All three must be
non-empty uppercase ACGT.
"""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.pa_workflow import assemble_linear_donor


def test_donor_is_concatenation():
    out = assemble_linear_donor(
        upstream_flap_dna="ACGT",
        bulk_dna="TTTT",
        downstream_flap_dna="GGGG",
    )
    assert out == "ACGTTTTTGGGG"


def test_donor_starts_with_upstream_flap():
    out = assemble_linear_donor(
        upstream_flap_dna="ACGTACGT", bulk_dna="C" * 10, downstream_flap_dna="G" * 8
    )
    assert out.startswith("ACGTACGT")


def test_donor_ends_with_downstream_flap():
    out = assemble_linear_donor(
        upstream_flap_dna="A" * 8, bulk_dna="C" * 10, downstream_flap_dna="GTACGTAC"
    )
    assert out.endswith("GTACGTAC")


def test_empty_flap_raises():
    with pytest.raises(ValueError):
        assemble_linear_donor(
            upstream_flap_dna="",
            bulk_dna="ACGT",
            downstream_flap_dna="ACGT",
        )


def test_empty_bulk_raises():
    with pytest.raises(ValueError):
        assemble_linear_donor(
            upstream_flap_dna="ACGT",
            bulk_dna="",
            downstream_flap_dna="ACGT",
        )


def test_non_acgt_raises():
    with pytest.raises(ValueError):
        assemble_linear_donor(
            upstream_flap_dna="ACGU",
            bulk_dna="ACGT",
            downstream_flap_dna="ACGT",
        )
    with pytest.raises(ValueError):
        assemble_linear_donor(
            upstream_flap_dna="acgt",
            bulk_dna="ACGT",
            downstream_flap_dna="ACGT",
        )
