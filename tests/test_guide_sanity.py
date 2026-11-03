"""One definition of "a usable guide spacer", shared by every path that picks a guide.

Before this module the PRIDICT path rejected a spacer with a homopolymer run of 5 or more
(``select_tjpe_ngrna``) while the genome-enumerated path applied no check at all. The same real
guide, ``GCTGCCACAGATGCCCCCAA`` (it contains CCCCC), was therefore rejected on one path and
accepted on the other, and a guide with a TTTT run, which Pol III reads as a stop signal, could
get through the enumerated path untouched.
"""

import pytest

from pegrna_pipeline.hpc.data_types import NgRNADesign
from pegrna_pipeline.hpc.guide_sanity import assess_spacer, spacer_passes_sanity
from pegrna_pipeline.hpc.nick_search import NickCandidate, filter_nick_sanity
from pegrna_pipeline.hpc.tjpe_workflow import select_tjpe_ngrna

NICK_101 = "GCTGCCACAGATGCCCCCAA"  # the CTCF example's clean nick; contains CCCCC
NICK_232 = "TTCTCCCTTCTTTCTTTTTT"  # the CTCF example's poly-T guide


def _gc(n_gc: int, length: int = 20) -> str:
    """A 20-mer with exactly ``n_gc`` G/C and no run longer than 2."""
    pattern = ["GA", "CT"]
    out = []
    for i in range(length):
        out.append(pattern[i % 2][0] if i < n_gc else pattern[i % 2][1])
    s = "".join(out)
    assert (s.count("G") + s.count("C")) == n_gc
    return s


# ---- hard rules ----------------------------------------------------------------------


def test_tttt_is_a_hard_failure_and_says_why():
    v = assess_spacer("ACGTTTTACGCAGGCATCGA")
    assert not v.ok
    assert any("TTTT" in f and "polymerase III" in f for f in v.failures)


def test_a_run_of_three_t_is_fine():
    assert spacer_passes_sanity("ACGTTTACGCAGGCATCGAC")


@pytest.mark.parametrize("n_gc, ok", [(5, False), (6, True), (15, True), (16, False)])
def test_gc_content_bounds_are_30_to_75_percent_inclusive(n_gc, ok):
    assert spacer_passes_sanity(_gc(n_gc)) is ok


def test_an_empty_spacer_fails():
    assert not assess_spacer("").ok


def test_case_does_not_matter():
    assert assess_spacer("acgtttt acgcaggcatcga".replace(" ", "")).ok is False
    assert assess_spacer(NICK_101.lower()).ok


# ---- advisories ----------------------------------------------------------------------


def test_a_homopolymer_run_is_an_advisory_not_a_rejection():
    v = assess_spacer(NICK_101)
    assert v.ok, "CCCCC must not disqualify a guide: only TTTT has a mechanism"
    assert len(v.advisories) == 1
    # G C T G C C A C A G A T G [C C C C C] A A: 1-based positions 14-18
    assert NICK_101[13:18] == "CCCCC"
    assert "CCCCC" in v.advisories[0] and "positions 14-18" in v.advisories[0]


@pytest.mark.parametrize("base", "ACG")
def test_runs_of_any_base_are_reported_not_rejected(base):
    spacer = "GCATGCA" + base * 6 + "TGCATG"
    spacer = spacer[:20]
    v = assess_spacer(spacer)
    assert v.ok and v.advisories


def test_a_run_of_four_is_not_even_an_advisory():
    assert assess_spacer("GCATGCACCCCTGCATGCAC").advisories == ()


def test_each_run_is_reported_separately():
    v = assess_spacer("AAAAAGCATGCGGGGGCATG")
    assert len(v.advisories) == 2


def test_advisories_never_affect_the_boolean_verdict():
    assert spacer_passes_sanity("AAAAAGCATGCGGGGGCATG")


# ---- the two paths agree -------------------------------------------------------------

SPACERS = [
    NICK_101,
    NICK_232,
    "AAAAAACCTGGTACGGATCA",  # A-run, otherwise fine
    "GGGGGACCTGGTACGGATCA",  # G-run
    "ATTTTACCTGGTACGGATCA",  # TTTT
    "GCGCGCGCGCGCGCGCGCGC",  # 100% GC
    "AAAATAATAATAAAATAAAT",  # ~0% GC
    "AGTGCTACCTGGTACGGATC",  # ordinary
]


def _via_pridict_path(spacer: str) -> bool:
    ng = NgRNADesign(
        spacer=spacer,
        distance_to_edit=50,
        target_strand="Fw",
        deepcas9_score=0.5,
        pe3b=False,
        pam_disrupt=False,
    )
    return select_tjpe_ngrna([ng]) is not None


def _via_genome_path(spacer: str) -> bool:
    cand = NickCandidate(
        protospacer=spacer, pam="TGG", nick_position=200, distance_to_edit=50, target_strand="Fw"
    )
    kept, _ = filter_nick_sanity([cand])
    return bool(kept)


@pytest.mark.parametrize("spacer", SPACERS)
def test_both_candidate_paths_judge_every_guide_identically(spacer):
    expected = spacer_passes_sanity(spacer)
    assert _via_pridict_path(spacer) is expected
    assert _via_genome_path(spacer) is expected


def test_the_guide_the_paths_used_to_disagree_about_is_now_accepted_by_both():
    assert _via_pridict_path(NICK_101) and _via_genome_path(NICK_101)


def test_a_tttt_guide_no_longer_slips_through_the_genome_path():
    assert not _via_genome_path(NICK_232)


def test_filter_reports_why_each_guide_was_dropped():
    cands = [
        NickCandidate(NICK_101, "TGG", 200, 50, "Fw"),
        NickCandidate(NICK_232, "TGG", 210, 60, "Fw"),
    ]
    kept, dropped = filter_nick_sanity(cands)
    assert [c.protospacer for c in kept] == [NICK_101]
    assert len(dropped) == 1 and dropped[0][0].protospacer == NICK_232
    assert any("TTTT" in reason for reason in dropped[0][1])
