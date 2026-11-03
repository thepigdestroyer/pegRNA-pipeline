"""Extension-to-scaffold pairing is measured and reported, but never rejects (issue #27, ADR 0026).

The folding gate leaves every pair that touches the scaffold out, by design. The docking gate
catches extensive sequestration but only 4 of 12 24-nt complements. So nothing flagged a
3' extension that ties up the scaffold. ViennaRNA separates the cases cleanly at the folding
gate's own cutoff: designs built to sequester pair 18-29 scaffold nucleotides with the extension,
synthetic neutral designs at most 11, and none of 130 real candidates of the CTCF pool above 9.
"""

import sys
from pathlib import Path

import pytest
import test_pridict_ngrna_wiring as W

from pegrna_pipeline.hpc.full_structure import (
    FullStructureFake,
    ScaffoldPairing,
    functional_region_pairing,
    scaffold_extension_pairing,
)
from pegrna_pipeline.hpc.gate_thresholds import THRESHOLDS
from pegrna_pipeline.hpc.recoded_tjpe_workflow import run_recoded_tjpe
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "analysis"))
import docking_reference as dr  # noqa: E402

ADVISORY = THRESHOLDS.scaffold_pairing_advisory_nt
REGIONS = {"spacer": (0, 20), "scaffold": (20, 96), "rtt": (96, 150)}


def _count(pairs, threshold=0.5):
    return scaffold_extension_pairing(
        "A" * 150, regions=REGIONS, scorer=FullStructureFake(pairs), threshold=threshold
    )


# ---- what is counted ------------------------------------------------------------------


def test_scaffold_to_extension_pairs_are_counted_by_distinct_scaffold_nucleotide():
    r = _count({(30, 100): 1.0, (31, 101): 1.0, (31, 120): 0.9})
    assert r == ScaffoldPairing(paired_nt=2, pairs=3)  # nt 31 has two partners but counts once


def test_only_pairs_between_the_scaffold_and_the_extension_count():
    r = _count(
        {
            (40, 60): 1.0,  # scaffold-scaffold: the scaffold's own fold
            (5, 110): 1.0,  # spacer-extension: the folding gate's business
            (5, 50): 1.0,  # spacer-scaffold
            (100, 140): 1.0,  # extension-extension
            (30, 100): 1.0,  # the one that counts
        }
    )
    assert r == ScaffoldPairing(paired_nt=1, pairs=1)


def test_the_threshold_is_respected_and_inclusive():
    pairs = {(30, 100): 0.6, (31, 101): 0.4}
    assert _count(pairs, threshold=0.5).paired_nt == 1
    assert _count(pairs, threshold=0.4).paired_nt == 2
    assert _count(pairs, threshold=0.7).paired_nt == 0


def test_a_design_without_a_scaffold_region_is_zero():
    r = scaffold_extension_pairing(
        "A" * 50, regions={"spacer": (0, 20)}, scorer=FullStructureFake({(1, 30): 1.0})
    )
    assert r.paired_nt == 0


def test_a_supplied_probability_map_is_reused_not_recomputed():
    class Boom:
        def pair_probs(self, seq):
            raise AssertionError("must not fold again when pairs are supplied")

    pairs = {(30, 100): 1.0}
    assert (
        scaffold_extension_pairing("A", regions=REGIONS, scorer=Boom(), pairs=pairs).paired_nt == 1
    )
    assert (
        functional_region_pairing(
            "A", regions=REGIONS, scorer=Boom(), pairs=pairs
        ).cross_region_pairs
        == 0
    )


# ---- it is invisible to the folding gate, which is the point ----------------------------


def test_the_folding_gate_does_not_see_what_the_advisory_sees():
    pairs = {(20 + k, 100 + k): 1.0 for k in range(25)}
    fake = FullStructureFake(pairs)
    gate = functional_region_pairing("A" * 150, regions=REGIONS, scorer=fake, threshold=0.5)
    assert gate.cross_region_pairs == 0
    assert _count(pairs).paired_nt == 25


# ---- in the walk: advisory only ---------------------------------------------------------


def _walk_with_scaffold_pairs(n):
    """The fake folder reports ``n`` scaffold<->RTT pairs on every candidate."""
    pairs = {(20 + k, 120 + k): 1.0 for k in range(n)}
    return run_recoded_tjpe(**W._wf_inputs(full_structure_scorer=FullStructureFake(pairs)))


def _advisories(result):
    return [w for w in result.warnings if "scaffold nucleotides" in w]


def test_heavy_scaffold_pairing_warns_and_does_not_reject():
    res = _walk_with_scaffold_pairs(ADVISORY + 10)
    assert res.pegRNA_sequence, "an advisory must never remove a design"
    (w,) = _advisories(res)
    assert f"{ADVISORY + 10} scaffold nucleotides" in w and "before ordering" in w


def test_the_advisory_starts_exactly_at_the_documented_threshold():
    assert _advisories(_walk_with_scaffold_pairs(ADVISORY - 1)) == []
    assert len(_advisories(_walk_with_scaffold_pairs(ADVISORY))) == 1


def test_no_pairing_no_warning():
    res = _walk_with_scaffold_pairs(0)
    assert _advisories(res) == []
    assert "extension pairs with 0 scaffold nt" in dict(res.gate_results)["structure"]


def test_the_structure_row_states_the_scaffold_count():
    res = _walk_with_scaffold_pairs(7)
    assert "extension pairs with 7 scaffold nt" in dict(res.gate_results)["structure"]


# ---- the calibration behind the threshold, pinned to the committed data ----------------


def _vienna():
    pytest.importorskip("RNA")
    from pegrna_pipeline.hpc.vienna_structure import ViennaRNAStructureScorer

    return ViennaRNAStructureScorer()


def _scaffold_nt(vienna, rna):
    seq = rna.replace("U", "T")
    regions = {"spacer": (0, 20), "scaffold": (20, 96), "rtt": (96, len(seq))}
    return scaffold_extension_pairing(
        seq, regions=regions, scorer=vienna, threshold=THRESHOLDS.vienna_pair_prob
    ).paired_nt


def test_the_threshold_separates_the_calibration_classes():
    vienna = _vienna()
    recs = dr.load()
    neutral = [
        _scaffold_nt(vienna, r["rna"]) for r in recs if r["class"] == "neutral" and r.get("rna")
    ]
    sweep = [
        _scaffold_nt(vienna, r["rna"])
        for r in recs
        if r["class"] == "sequestering" and r.get("rna")
    ]
    assert len(neutral) >= 30 and len(sweep) >= 10
    assert max(neutral) < ADVISORY <= min(sweep), (max(neutral), ADVISORY, min(sweep))


def test_the_pipelines_own_designs_are_not_flagged():
    vienna = _vienna()
    designed = [r for r in dr.load() if r["class"] == "designed" and r.get("rna")]
    assert designed
    for r in designed:
        assert _scaffold_nt(vienna, r["rna"]) < ADVISORY


def test_a_cargo_that_contains_the_scaffold_is_flagged():
    """The realistic way to hit this: an insert that carries a guide-RNA scaffold, whose reverse
    complement then sits in the extension and pairs with the pegRNA's own scaffold."""
    vienna = _vienna()
    spacer = "TCAGCATGATGGACCGGTGA"
    rc = STANDARD_SGRNA_SCAFFOLD[::-1].translate(str.maketrans("ACGT", "TGCA"))
    pegrna = spacer + STANDARD_SGRNA_SCAFFOLD + "CGTACGTACGTACG" + rc[:60] + "GCATCGATCGGCATCGA"
    assert _scaffold_nt(vienna, pegrna) >= ADVISORY


def test_none_of_the_real_candidate_pool_is_flagged():
    """The advisory is only useful if it stays quiet on ordinary designs. The pegRNAs the pipeline
    itself produces for CTCF (several cargos, all four 5' nicks) never come near the threshold:
    130 distinct designs gave a maximum of 9. (At p >= 0.5 the same designs would be flagged,
    which is why the advisory uses the folding gate's own 0.99.)"""
    vienna = _vienna()
    import json

    import sample_pool_candidates as sp

    cfg = json.loads(sp.LOCUS_CONFIG.read_text())
    cargos = {name: parts for name, parts in sp.CARGOS}
    seen = []
    for name in ("FLAG", "V5-P2A"):
        seen += sp.collect_pool(cfg, sp.SPACERS[0], sp.build_cargo(cargos[name]))
    assert len(seen) >= 12
    worst = max(_scaffold_nt(vienna, peg) for peg in seen)
    assert worst < ADVISORY, f"a real design reaches {worst} scaffold nt paired"
