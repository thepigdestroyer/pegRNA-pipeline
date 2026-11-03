"""The pegRNA-folding gate's threshold belongs to the folder that produced it.

The generic default (0.5) is RiNALMo's pair-calling boundary. ViennaRNA's
partition-function probabilities are distributed differently, and under it a 0.5
cutoff on SINGLE cross-region base pairs rejected ~90% of RANDOM tjPE-sized
pegRNAs, and 94-100% of the real candidates for the CTCF example: a gate that
rejects nearly everything is not screening out bad designs. At 0.99 it removes
only the worst tail. The ViennaRNA scorer therefore carries its own calibrated
default, and the workflow uses it unless told otherwise.
"""

import random

import pytest

from pegrna_pipeline.hpc.full_structure import FullStructureFake, functional_region_pairing
from pegrna_pipeline.hpc.gate_thresholds import THRESHOLDS
from pegrna_pipeline.hpc.recoded_tjpe_workflow import resolve_structure_threshold
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD
from pegrna_pipeline.hpc.vienna_structure import ViennaRNAStructureScorer


def test_explicit_threshold_always_wins():
    assert resolve_structure_threshold(0.7, ViennaRNAStructureScorer()) == 0.7
    assert resolve_structure_threshold(0.0, ViennaRNAStructureScorer()) == 0.0  # falsy but explicit


def test_vienna_scorer_supplies_its_own_default():
    assert resolve_structure_threshold(None, ViennaRNAStructureScorer()) == 0.99


def test_vienna_default_is_the_recorded_calibrated_value():
    assert ViennaRNAStructureScorer.recommended_threshold == THRESHOLDS.vienna_pair_prob
    assert (
        "vienna_pair_prob" in THRESHOLDS.citations
    ), "a threshold without its evidence is arbitrary"


def test_scorers_without_a_recommendation_keep_the_generic_default():
    assert resolve_structure_threshold(None, FullStructureFake(pairs={})) == 0.5


# --- the calibration evidence itself ---------------------------------------------------

SPACER = "TCAGCATGATGGACCGGTGA"


def _random_design(rng, n_rcpbs2=12, n_rtt=70, n_pbs1=13):
    def rnd(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    rcpbs2, rtt, pbs1 = rnd(n_rcpbs2), rnd(n_rtt), rnd(n_pbs1)
    scaf = STANDARD_SGRNA_SCAFFOLD
    i1 = len(SPACER)
    i2 = i1 + len(scaf)
    i3 = i2 + len(rcpbs2)
    i4 = i3 + len(rtt)
    seq = SPACER + scaf + rcpbs2 + rtt + pbs1
    regions = {
        "spacer": (0, i1),
        "scaffold": (i1, i2),
        "rcpbs2": (i2, i3),
        "rtt": (i3, i4),
        "pbs1": (i4, len(seq)),
    }
    return seq, regions


def test_the_old_threshold_rejects_most_random_designs_and_the_new_one_does_not():
    pytest.importorskip("RNA")  # ViennaRNA; skipped where it is not installed (CI)
    vienna = ViennaRNAStructureScorer()
    rng = random.Random(7)
    designs = [_random_design(rng) for _ in range(80)]

    def pass_rate(threshold):
        ok = sum(
            functional_region_pairing(
                s, regions=r, scorer=vienna, threshold=threshold
            ).cross_region_pairs
            == 0
            for s, r in designs
        )
        return ok / len(designs)

    old, new = pass_rate(0.5), pass_rate(0.99)
    # Measured on 300 designs: 9% and 80%. Generous bounds for an 80-design sample.
    assert old < 0.30, f"0.5 should reject most random tjPE-sized designs (passed {old:.0%})"
    assert new > 0.65, f"0.99 should reject only the worst tail (passed {new:.0%})"
