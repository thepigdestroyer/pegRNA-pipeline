"""A folding-gate failure says WHAT paired, in words, and names the folder that said so.

The CTCF cassette run was rejected on every candidate with five copies of "RiNALMo cross-region pairing". The folder
was ViennaRNA, and nothing said which parts of the pegRNA were involved. The cause was one 4-bp stem: the first
four bases of the rcPBS2 (GGGG) paired at p >= 0.996 with the first four of the RTT (CCCC), which is the reverse
complement of the cassette's last four bases (a T7 promoter remnant ending in GGGG). Trimming two bases of that tail
cleared it. The user could not have found that from the old message.

The gate's COUNT and its verdict are unchanged; only what it reports is new.
"""

import pytest
import test_pridict_ngrna_wiring as W

from pegrna_pipeline.hpc.full_structure import (
    FullStructureFake,
    Helix,
    StructureScore,
    functional_region_pairing,
)
from pegrna_pipeline.hpc.recoded_tjpe_workflow import (
    NoValidPegRNA,
    _folder_name,
    _is_insert_side_failure,
    run_recoded_tjpe,
)

# spacer 0-20, scaffold 20-96, rcPBS2 96-108, RTT 108-150, PBS1 150-163
REGIONS = {
    "spacer": (0, 20),
    "scaffold": (20, 96),
    "rcpbs2": (96, 108),
    "rtt": (108, 150),
    "pbs1": (150, 163),
}


def _peg() -> str:
    seq = ["A"] * 163
    seq[96:100] = "GGGG"  # rcPBS2 5' end
    seq[108:112] = "CCCC"  # RTT 5' end
    return "".join(seq)


# the four stacked pairs rcPBS2[0..3] x RTT[3..0]: (96,111) (97,110) (98,109) (99,108)
STEM = {(96, 111): 0.996, (97, 110): 0.998, (98, 109): 0.999, (99, 108): 0.997}


def _score(pairs, threshold=0.99, seq=None):
    return functional_region_pairing(
        seq if seq is not None else _peg(),
        regions=REGIONS,
        scorer=FullStructureFake(pairs),
        threshold=threshold,
    )


def test_a_stem_between_two_regions_is_reported_as_one_helix_with_its_sequences():
    score = _score(STEM)
    assert score.cross_region_pairs == 4
    (h,) = score.helices
    assert (h.region_a, h.a_lo, h.a_hi) == ("rcpbs2", 0, 3)
    assert (h.region_b, h.b_lo, h.b_hi) == ("rtt", 0, 3)
    assert h.bp == 4 and h.min_prob == pytest.approx(0.996)
    assert (h.seq_a, h.seq_b) == ("GGGG", "CCCC")


def test_the_description_reads_as_a_sentence_a_person_can_check():
    assert _score(STEM).describe() == "rcPBS2[0..3] GGGG x RTT[0..3] CCCC (4 bp, p >= 0.996)"


def test_two_separate_stems_are_two_helices_in_position_order():
    pairs = {**STEM, (20 - 5, 150): 0.999}  # spacer[15] x PBS1[0]: the intrinsic pair, not counted
    pairs.update(
        {(112, 130): 0.995, (113, 129): 0.995, (114, 128): 0.995}
    )  # RTT x RTT: not counted
    pairs.update({(116, 155): 0.999})  # RTT[8] x PBS1[5]
    score = _score(pairs)
    assert score.cross_region_pairs == 5
    assert [h.describe().split(" ")[0] for h in score.helices] == ["rcPBS2[0..3]", "RTT[8..8]"]
    assert [h.bp for h in score.helices] == [4, 1]


def test_pairs_the_gate_does_not_count_are_not_described():
    spacer_pbs1 = {(10, 155): 0.999}  # the expected spacer x PBS1 pair
    scaffold = {(30, 120): 0.999}  # anything touching the scaffold
    inside = {(112, 130): 0.999}  # within one region
    low = {(96, 111): 0.5}  # under the cutoff
    score = _score({**spacer_pbs1, **scaffold, **inside, **low})
    assert score == StructureScore(cross_region_pairs=0, helices=())
    assert score.describe() == ""


def test_the_count_is_exactly_what_it_was_before_helices_existed():
    for pairs in ({}, STEM, {**STEM, (116, 155): 0.999}, {(96, 111): 0.989}):
        assert _score(pairs).cross_region_pairs == sum(
            1 for (i, j), p in pairs.items() if p >= 0.99 and i != j
        )


def test_a_long_list_is_cut_off_and_says_how_many_more():
    pairs = {(96 + k, 150 + k): 0.999 for k in range(0, 12, 2)}  # six unstacked single pairs
    score = _score(pairs)
    assert len(score.helices) == 6
    assert score.describe(limit=3).endswith("+3 more")


def test_a_helix_is_a_plain_value():
    h = Helix("rtt", 1, 2, "pbs1", 3, 4, 2, 0.99, "AC", "GT")
    assert h.describe() == "RTT[1..2] AC x PBS1[3..4] GT (2 bp, p >= 0.990)"


# ---- the walk ------------------------------------------------------------------------------------


def _l_rcpbs2() -> int:
    return len(run_recoded_tjpe(**W._wf_inputs()).rcpbs2)


def _failing_inputs():
    """The same stem, at the real coordinates of this fixture's pegRNA (rcPBS2 starts at 96)."""
    r_start = 96 + _l_rcpbs2()  # the RTT begins where the rcPBS2 ends
    pairs = {(96 + k, r_start + 3 - k): 0.998 for k in range(4)}
    return W._wf_inputs(full_structure_scorer=FullStructureFake(pairs=pairs))


def test_the_failure_names_the_folder_and_the_helix():
    with pytest.raises(NoValidPegRNA) as exc:
        run_recoded_tjpe(**_failing_inputs())
    first = next(d for d in exc.value.diagnostics if d.startswith("cand1:"))
    assert first.startswith("cand1: cross-region pairing (FullStructureFake): rcPBS2[0..3] ")
    assert "x RTT[0..3]" in first and "(4 bp, p >= 0.998)" in first


def test_the_failure_is_still_an_insert_side_failure_so_the_insert_lever_still_fires():
    with pytest.raises(NoValidPegRNA) as exc:
        run_recoded_tjpe(**_failing_inputs())
    assert _is_insert_side_failure(exc.value.diagnostics)


def test_the_progress_event_carries_the_helix_where_the_web_readout_shows_it():
    """The readout prints ``extra`` under the count when there is a cutoff, and ``info`` when there is not."""
    events = []
    with pytest.raises(NoValidPegRNA):
        run_recoded_tjpe(
            **_failing_inputs(),
            on_gate=lambda cand, gate, info, passed, metrics=None: events.append(
                (gate, passed, info, metrics or {})
            ),
        )
    gate, passed, info, metrics = next(e for e in events if e[0] == "structure")
    assert passed is False and "4 cross-region pairs" in info and "rcPBS2[0..3]" in info
    assert "rcPBS2[0..3]" in metrics["extra"] and metrics["value"] == 4


def test_a_passing_structure_says_nothing_extra():
    events = []
    run_recoded_tjpe(
        **W._wf_inputs(),
        on_gate=lambda cand, gate, info, passed, metrics=None: events.append((gate, passed, info)),
    )
    gate, passed, info = next(e for e in events if e[0] == "structure")
    assert passed is True and info.endswith("cross-region pairs (p>=0.5)") and ":" not in info


@pytest.mark.parametrize(
    "scorer_name, expected",
    [
        ("ViennaRNAStructureScorer", "ViennaRNA"),
        ("RiNALMoScorer", "RiNALMo"),
        ("RNAformerScorer", "RNAformer"),
        ("FullStructureFake", "FullStructureFake"),
    ],
)
def test_the_folder_is_named_for_what_it_is(scorer_name, expected):
    assert _folder_name(type(scorer_name, (), {})()) == expected
