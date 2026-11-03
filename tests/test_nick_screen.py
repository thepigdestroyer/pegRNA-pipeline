"""Task 6: off-target / on-target screen for nick candidates (+other-allele)."""

from pegrna_pipeline.hpc.nick_search import NickCandidate, screen_nick_candidates


def _fake_scorer(protospacers):
    # {protospacer: (ontarget, offtarget_count, has_offtarget)}
    return {p: (0.7, 0, False) if "ACGT" in p else (0.2, 9, True) for p in protospacers}


def test_screen_drops_offtarget_and_low_ontarget():
    cands = [
        NickCandidate("ACGT" * 5, "AGG", 10, 45, "Fw"),  # clean
        NickCandidate("TTTT" * 5, "AGG", 12, 47, "Fw"),  # offtarget-heavy
    ]
    kept = screen_nick_candidates(cands, scorer=_fake_scorer, min_ontarget=0.5, max_offtargets=0)
    assert len(kept) == 1
    assert kept[0].protospacer == "ACGT" * 5


def test_screen_drops_candidate_matching_other_allele():
    # other allele identical => would cut both alleles => not allele-specific => drop
    clean = NickCandidate("ACGT" * 5, "AGG", 10, 45, "Fw")
    kept = screen_nick_candidates(
        [clean],
        scorer=_fake_scorer,
        min_ontarget=0.5,
        max_offtargets=0,
        other_allele_seqs={"ACGT" * 5: "ACGT" * 5},
    )
    assert kept == []


def test_screen_keeps_candidate_when_other_allele_differs():
    clean = NickCandidate("ACGT" * 5, "AGG", 10, 45, "Fw")
    kept = screen_nick_candidates(
        [clean],
        scorer=_fake_scorer,
        min_ontarget=0.5,
        max_offtargets=0,
        other_allele_seqs={"ACGT" * 5: "ACGA" * 5},  # SNP differs => discriminating
    )
    assert len(kept) == 1


def test_screen_sorts_by_ontarget_desc():
    def scorer(ps):
        return {p: (0.6 if p.startswith("ACGTA") else 0.9, 0, False) for p in ps}

    a = NickCandidate("ACGTA" + "C" * 15, "AGG", 10, 45, "Fw")  # 0.6
    b = NickCandidate("CCCCC" + "G" * 15, "AGG", 11, 46, "Fw")  # 0.9
    kept = screen_nick_candidates([a, b], scorer=scorer)
    # Survivors are copies (they now carry the off-target count), so compare by content.
    assert [c.protospacer for c in kept] == [b.protospacer, a.protospacer]
    assert [c.offtargets for c in kept] == [0, 0]
