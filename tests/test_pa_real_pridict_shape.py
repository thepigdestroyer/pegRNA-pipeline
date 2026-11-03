"""The legacy Prime Assembly pair filter cannot pair real PRIDICT output (issue #29, ADR 0027).

These tests pin a verified FACT, not a desired behaviour. Running real PRIDICT2 on the CTCF window
for two flaps inserted at the same site (``scripts/analysis/probe_legacy_pa.py``, on the DSI box)
returned 765 pegRNAs per flap, and the same ``(strand, Editing_Position)`` set for both:

    Fw: 0, 6, 11     Rv: 2, 17

``Editing_Position`` is the edit's distance FROM THE NICK, not a coordinate. The legacy filter
treats it as a coordinate on a shared axis and needs ``|up - down|`` of 30-80, so nothing can pass.
The hand-written fixtures (50/58 against 100/110) never could have come from PRIDICT, and hid it.

When Prime Assembly gets a real designer these tests should be deleted, not "fixed": they describe
the prototype's limit. Until then they stop anyone mistaking an empty result for a bad locus.
"""

import itertools

import pytest

from pegrna_pipeline.hpc.data_types import PegRNADesign
from pegrna_pipeline.hpc.pair_scorer import MIN_NICK_SEPARATION_BP, filter_compatible_pairs

# Observed from real PRIDICT2 on the CTCF window (same for both flap runs).
OBSERVED = {"Fw": (0, 6, 11), "Rv": (2, 17)}
LEGACY_MIN_DISTANCE = 30


def _design(strand: str, position: int, pbs: int, score: float) -> PegRNADesign:
    return PegRNADesign(
        sequence_name=f"d_{strand}_{position}_{pbs}",
        spacer=("ACGT" * 5) if strand == "Fw" else ("TGCA" * 5),
        pbs="A" * pbs,
        rt="C" * 44,
        pbs_length=pbs,
        rt_length=44,
        target_strand=strand,
        editing_position=position,
        pe2_hek293t=score,
        pe2_k562=score,
        pegrna_full="A" * 150,
    )


def _real_shaped_run(score: float) -> list:
    return [
        _design(strand, pos, pbs, score)
        for strand, positions in OBSERVED.items()
        for pos in positions
        for pbs in range(8, 16)
    ]


def test_real_pridict_positions_are_distances_from_the_nick_not_coordinates():
    assert max(p for ps in OBSERVED.values() for p in ps) < LEGACY_MIN_DISTANCE
    assert LEGACY_MIN_DISTANCE > MIN_NICK_SEPARATION_BP


def test_the_legacy_filter_finds_no_pair_in_real_shaped_output():
    up, dn = _real_shaped_run(9.0), _real_shaped_run(8.0)
    assert len(up) == len(dn) == 40
    assert filter_compatible_pairs(upstream=up, downstream=dn, stop_codon_position=147) == []


def test_no_choice_of_designs_could_ever_pass_the_legacy_window():
    """Not an unlucky draw: the largest possible difference between two observed values is
    below the smallest allowed distance, so every one of the cross-product pairs fails."""
    positions = [p for ps in OBSERVED.values() for p in ps]
    widest = max(abs(a - b) for a, b in itertools.product(positions, positions))
    assert widest < LEGACY_MIN_DISTANCE


def test_nick_spacing_for_one_insertion_site_is_bounded_by_what_pridict_proposes():
    """PRIDICT only proposes pegRNAs whose nick is within ~17 bp of the edit, so two of them for
    the same site have nicks at most ~28 bp apart. A Prime Assembly pair needs nicks flanking
    the region to be replaced; that spacing cannot come from PRIDICT's insert-at-S designs."""
    widest_nick_to_nick = max(OBSERVED["Fw"]) + max(OBSERVED["Rv"])
    assert widest_nick_to_nick < LEGACY_MIN_DISTANCE


@pytest.mark.parametrize("up_pos, dn_pos", [(50, 100), (58, 110)])
def test_the_hand_written_fixture_values_do_pass_which_is_why_the_fixtures_hid_it(up_pos, dn_pos):
    up = [_design("Fw", up_pos, 13, 9.0)]
    dn = [_design("Rv", dn_pos, 13, 8.0)]
    assert len(filter_compatible_pairs(upstream=up, downstream=dn, stop_codon_position=0)) == 1
