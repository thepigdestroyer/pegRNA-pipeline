"""Tests for the :class:`PAResult` dataclass field set.

Plan 2 v1.1 spec — these field names are canonical and tests below
should fail loudly if any get renamed or dropped.
"""

from __future__ import annotations

import dataclasses

from pegrna_pipeline.hpc.pa_workflow import PAResult


def test_pa_result_has_all_required_fields():
    fields = {f.name for f in dataclasses.fields(PAResult)}
    expected = {
        "winning_pair",
        "top_pairs",
        "linear_donor_dna",
        "upstream_pegRNA_sequence",
        "downstream_pegRNA_sequence",
        "upstream_flap_dna",
        "downstream_flap_dna",
        "bulk_dna",
        "full_integrated_product",
        "upstream_pegRNA_off_target",
        "downstream_pegRNA_off_target",
        "total_pairs_evaluated",
        "pairs_after_compat_filter",
    }
    assert expected <= fields, f"missing fields: {expected - fields}"
