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
