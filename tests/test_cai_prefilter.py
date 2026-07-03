"""Stage B: CAI enumerate/rank + top-N PRIDICT pre-filter (scalable k-best).

`rank_insert_variants_by_cai` must return the top synonymous variants of an insert
in CAI-descending order WITHOUT enumerating the full product (a 40 bp end is ~10^5-10^7
combos; a 500 bp middle is ~3^140). Correctness is pinned against brute force on a small
case; scalability is pinned by a long-sequence case that would hang under full enumeration.
"""

from __future__ import annotations

import time

import pytest
from Bio.Seq import Seq

from pegrna_pipeline.hpc.codon_expander import (
    CodonVariant,
    _surviving_variants,
    human_cai,
    rank_insert_variants_by_cai,
    select_top_n_for_pridict,
    split_insert,
)


def _translate(dna: str) -> str:
    return str(Seq(dna).translate())


# --------------------------------------------------------------- B1: rank by CAI
def test_rank_empty_returns_empty():
    assert rank_insert_variants_by_cai("") == []


def test_rank_descending_cai_and_translates():
    ranked = rank_insert_variants_by_cai("LR", max_variants=1000)
    assert ranked
    cais = [human_cai(v.dna) for v in ranked]
    assert cais == sorted(cais, reverse=True)  # non-increasing CAI
    for v in ranked:
        assert isinstance(v, CodonVariant)
        assert _translate(v.dna) == "LR"


def test_rank_matches_bruteforce_set_over_all_passers():
    # With max_variants covering everything, k-best must return exactly the set of
    # all filter-passing synonymous variants (completeness), ordered by CAI desc.
    aa = "LR"
    ranked = rank_insert_variants_by_cai(aa, max_variants=10_000)
    brute = {dna for dna, _c, _gc in _surviving_variants(aa)}
    assert {v.dna for v in ranked} == brute
    assert human_cai(ranked[0].dna) == max(human_cai(d) for d in brute)  # #1 is global max


def test_rank_applies_quality_filters():
    ranked = rank_insert_variants_by_cai("LLLL", max_variants=500)
    assert ranked
    for v in ranked:
        assert "CG" not in v.dna  # no CpG
        assert not any(  # no homopolymer run >= 5
            len(set(v.dna[i : i + 5])) == 1 for i in range(len(v.dna) - 4)
        )


def test_rank_caps_at_max_variants():
    assert len(rank_insert_variants_by_cai("LLLL", max_variants=3)) <= 3


def test_rank_deterministic():
    a = rank_insert_variants_by_cai("LRLR", max_variants=20)
    b = rank_insert_variants_by_cai("LRLR", max_variants=20)
    assert [v.dna for v in a] == [v.dna for v in b]


def test_rank_scales_without_full_enumeration():
    # 3^30 synonymous combos; must return top-50 fast via k-best, not enumerate all.
    t0 = time.time()
    ranked = rank_insert_variants_by_cai("L" * 30, max_variants=50)
    elapsed = time.time() - t0
    assert len(ranked) == 50
    cais = [human_cai(v.dna) for v in ranked]
    assert cais == sorted(cais, reverse=True)
    assert elapsed < 10, f"took {elapsed:.1f}s -- likely enumerating the full product"


def test_rank_assigns_sequential_ranks():
    ranked = rank_insert_variants_by_cai("LR", max_variants=5)
    assert [v.rank for v in ranked] == list(range(len(ranked)))


def test_rank_truncation_returns_fewer_without_hanging():
    # A tiny exploration cap forces truncation: returns fewer than requested, bounded,
    # no hang/crash (the documented restrictive-filter path).
    ranked = rank_insert_variants_by_cai("LLLL", max_variants=1000, max_examined=5)
    assert 0 <= len(ranked) < 1000


def test_rank_unknown_aa_raises_regardless_of_max_variants():
    # Invalid residue must fail loud whether or not max_variants short-circuits.
    with pytest.raises(KeyError):
        rank_insert_variants_by_cai("X", max_variants=50)
    with pytest.raises(KeyError):
        rank_insert_variants_by_cai("X", max_variants=0)


# ----------------------------------------------------- B2: select_top_n_for_pridict
def _fake_ranked(n: int):
    return [CodonVariant(dna="ACG", c_count=1, gc_content=0.66, rank=i) for i in range(n)]


def test_select_top_n_default_is_2000():
    assert len(select_top_n_for_pridict(_fake_ranked(6000))) == 2000


def test_select_top_n_clamps_below_1000():
    assert len(select_top_n_for_pridict(_fake_ranked(6000), n=10)) == 1000


def test_select_top_n_clamps_above_5000():
    assert len(select_top_n_for_pridict(_fake_ranked(9000), n=99999)) == 5000


def test_select_top_n_cannot_exceed_available():
    assert len(select_top_n_for_pridict(_fake_ranked(1200), n=2000)) == 1200


def test_select_top_n_preserves_order_prefix():
    ranked = _fake_ranked(6000)
    assert select_top_n_for_pridict(ranked, n=1500) == ranked[:1500]


# ------------------------------------------------------------------ B3: split_insert
def test_split_reconstructs_long_insert():
    insert = "ACGT" * 30  # 120 nt
    first, middle, last = split_insert(insert, end_len=40)
    assert first + middle + last == insert
    assert (len(first), len(middle), len(last)) == (40, 40, 40)


def test_split_reconstructs_medium_insert():
    insert = "ACGTACGTAC" * 6  # 60 nt, between end_len and 2*end_len
    first, middle, last = split_insert(insert, end_len=40)
    assert first + middle + last == insert
    assert len(first) == 40 and middle == "" and len(last) == 20


def test_split_reconstructs_short_insert():
    insert = "ACGTACGT"  # <= end_len
    first, middle, last = split_insert(insert, end_len=40)
    assert (first, middle, last) == (insert, "", "")
    assert first + middle + last == insert
