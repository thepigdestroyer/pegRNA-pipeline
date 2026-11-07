"""Task 10: full-pegRNA 2D structure check (RiNALMo primary, scaffold-excluded)."""

import pytest

from pegrna_pipeline.hpc.full_structure import (
    FullStructureFake,
    functional_region_pairing,
    region_accessibility,
)

# authoritative region order: spacer, scaffold, pbs2, rtt, pbs1 (decision 0)
_REGIONS = {
    "spacer": (0, 20),
    "scaffold": (20, 80),
    "pbs2": (80, 95),
    "rtt": (95, 120),
    "pbs1": (120, 133),
}


def test_functional_pairing_excludes_scaffold():
    fake = FullStructureFake(pairs={(2, 90): 0.9, (25, 30): 0.9})  # (25,30) in scaffold
    score = functional_region_pairing("A" * 140, regions=_REGIONS, scorer=fake, threshold=0.5)
    assert score.cross_region_pairs == 1  # only spacer<->pbs2


def test_clean_pegrna_passes():
    fake = FullStructureFake(pairs={(25, 30): 0.9})  # scaffold self-pair only
    score = functional_region_pairing("A" * 140, regions=_REGIONS, scorer=fake, threshold=0.5)
    assert score.cross_region_pairs == 0


def test_below_threshold_pairs_ignored():
    fake = FullStructureFake(pairs={(2, 90): 0.3})  # below 0.5
    score = functional_region_pairing("A" * 140, regions=_REGIONS, scorer=fake, threshold=0.5)
    assert score.cross_region_pairs == 0


def test_intra_region_pairs_not_counted():
    fake = FullStructureFake(pairs={(96, 110): 0.9})  # both inside rtt
    score = functional_region_pairing("A" * 140, regions=_REGIONS, scorer=fake, threshold=0.5)
    assert score.cross_region_pairs == 0


def test_region_accessibility_from_prob_matrix():
    fake = FullStructureFake(pairs={(0, 50): 0.8})
    acc = region_accessibility("A" * 60, region=(0, 5), scorer=fake)
    assert 0.0 <= acc <= 1.0
    # base 0 is 0.8 paired => ~0.2 unpaired; bases 1-4 fully accessible => mean high
    assert acc > 0.8


def test_real_rinalmo_adapter_errors_without_gpu():
    from pegrna_pipeline.hpc.full_structure import RiNALMoScorer

    with pytest.raises((ImportError, RuntimeError)):
        RiNALMoScorer().pair_probs("ACGT" * 10)


# ---- expected_pairs: position pairs a design makes inevitable are not counted, and nothing else is -----------


def _count(pairs, **kw):
    fake = FullStructureFake(pairs=pairs)
    return functional_region_pairing("A" * 140, regions=_REGIONS, scorer=fake, threshold=0.5, **kw)


def test_an_expected_pair_is_not_counted():
    assert _count({(2, 100): 0.9}).cross_region_pairs == 1  # spacer <-> rtt
    assert _count({(2, 100): 0.9}, expected_pairs=[(2, 100)]).cross_region_pairs == 0


def test_an_expected_pair_is_matched_whichever_way_round_it_is_given():
    assert _count({(2, 100): 0.9}, expected_pairs=[(100, 2)]).cross_region_pairs == 0
    assert _count({(100, 2): 0.9}, expected_pairs=[(2, 100)]).cross_region_pairs == 0


def test_only_the_exact_pair_is_exempt_not_its_neighbours_or_its_regions():
    pairs = {(2, 100): 0.9, (3, 100): 0.9, (2, 101): 0.9, (4, 99): 0.9}
    got = _count(pairs, expected_pairs=[(2, 100)])
    assert (
        got.cross_region_pairs == 3
    ), "one pair exempted; the three that differ by a position are still counted"


def test_an_expected_pair_does_not_hide_a_second_helix_in_the_same_regions():
    got = _count({(2, 100): 0.9, (10, 110): 0.9}, expected_pairs=[(2, 100)])
    assert got.cross_region_pairs == 1 and got.helices[0].a_lo == 10


def test_without_expected_pairs_nothing_changes():
    pairs = {(2, 90): 0.9, (25, 30): 0.9, (2, 100): 0.9, (96, 110): 0.9, (2, 125): 0.9}
    assert _count(pairs).cross_region_pairs == _count(pairs, expected_pairs=None).cross_region_pairs
    assert (
        _count(pairs).cross_region_pairs == _count(pairs, expected_pairs=[]).cross_region_pairs == 2
    )
