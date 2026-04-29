"""Tests for the :class:`TJPEResult` dataclass (v1.1 field names).

Per Plan 2 v1.1 the result carries:
    * ``edit_position_in_mrna`` (previously misnamed ``first_nick_position``)
    * ``ngrna_nick_in_mrna`` (previously misnamed ``second_nick_position``)
    * ``ngrna_distance_to_edit`` (previously misnamed ``ngrna_distance``)

These tests guard the exact field names so the renaming sticks.
"""

from __future__ import annotations

import dataclasses

from pegrna_pipeline.hpc.tjpe_workflow import TJPEResult


def test_tjpe_result_has_v11_field_names():
    fields = {f.name for f in dataclasses.fields(TJPEResult)}
    # Renamed fields (v1.1):
    assert "edit_position_in_mrna" in fields
    assert "ngrna_nick_in_mrna" in fields
    assert "ngrna_distance_to_edit" in fields
    # Old field names must NOT exist (catches copy-paste regressions).
    assert "first_nick_position" not in fields
    assert "second_nick_position" not in fields
    assert "ngrna_distance" not in fields


def test_tjpe_result_carries_the_full_spec():
    fields = {f.name for f in dataclasses.fields(TJPEResult)}
    expected = {
        "winning_variant",
        "full_insert_dna",
        "tjpe_pegRNA_sequence",
        "spacer",
        "pbs1",
        "rc_pbs2",
        "rc_pbs2_length",
        "ngrna_spacer",
        "ngrna_distance_to_edit",
        "ngrna_deepcas9_score",
        "pegrna_pe2_score",
        "pegrna_pe2_cell_line",
        "edit_position_in_mrna",
        "ngrna_nick_in_mrna",
        "all_evaluated_variants",
    }
    assert expected <= fields, f"missing fields: {expected - fields}"
