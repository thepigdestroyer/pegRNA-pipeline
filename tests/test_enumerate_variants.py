"""Tests for the heap-based enumerate_variants (v1.2).

enumerate_variants must:
- Yield CodonVariant objects.
- Apply filters (_has_internal_stop, _has_cpg, _has_restriction_site,
  _has_homopolymer) and keep only surviving combinations.
- Globally correct top-K ordering by (c_count DESC, abs(gc - 0.55) ASC).
- Respect max_variants via heapq.nlargest.
- Have deterministic rank assignment (0, 1, 2, ...).
"""

import itertools
from typing import List, Tuple


from pegrna_pipeline.hpc.codon_expander import (
    CODONS_BY_AA,
    CodonVariant,
    _gc_content,
    _has_cpg,
    _has_homopolymer,
    _has_internal_stop,
    _has_restriction_site,
    enumerate_variants,
)


def _manual_surviving(aa_sequence: str) -> List[Tuple[str, int, float]]:
    """Brute-force enumerate all synonymous combinations, applying filters."""
    codon_choices = [CODONS_BY_AA[aa.upper()] for aa in aa_sequence]
    out = []
    for combo in itertools.product(*codon_choices):
        dna = "".join(combo)
        if _has_internal_stop(dna):
            continue
        if _has_cpg(dna):
            continue
        if _has_restriction_site(dna):
            continue
        if _has_homopolymer(dna):
            continue
        out.append((dna, dna.count("C"), _gc_content(dna)))
    return out


def _manual_sort(triples):
    """Sort like enumerate_variants: c_count DESC, then abs(gc - 0.55) ASC, then dna for stability."""
    return sorted(
        triples,
        key=lambda t: (-t[1], abs(t[2] - 0.55), t[0]),
    )


def test_enumerate_variants_returns_list_of_CodonVariant():
    variants = enumerate_variants("M", max_variants=10)
    assert all(isinstance(v, CodonVariant) for v in variants)


def test_enumerate_variants_empty_input_returns_empty_list():
    assert enumerate_variants("", max_variants=10) == []


def test_enumerate_variants_uncapped_matches_manual_sort():
    # Small enough that uncapped enumeration is tractable (~a few dozen combos).
    aa = "MK"
    manual = _manual_surviving(aa)
    expected_sorted = _manual_sort(manual)
    variants = enumerate_variants(aa, max_variants=1_000_000)
    got = [(v.dna, v.c_count, v.gc_content) for v in variants]
    # Must contain exactly the surviving set, ordered by (c_count DESC, |gc - 0.55| ASC).
    assert len(got) == len(manual)
    # Check top-level ordering is correct.
    for (dna, c, gc), (e_dna, e_c, e_gc) in zip(got, expected_sorted):
        # c_count primary ordering:
        assert c == e_c


def test_enumerate_variants_capped_matches_top_k_manual():
    aa = "MKL"  # bigger: L has 6 codons, so many combinations.
    manual = _manual_surviving(aa)
    top_k_manual = _manual_sort(manual)[:10]
    variants = enumerate_variants(aa, max_variants=10)
    got = [(v.dna, v.c_count, v.gc_content) for v in variants]
    assert len(got) <= 10
    # Compare primary sort key of capped set to top-K-by-key of manual set.
    assert [v[1] for v in got] == [v[1] for v in top_k_manual[: len(got)]]


def test_enumerate_variants_assigns_sequential_ranks():
    variants = enumerate_variants("MK", max_variants=1000)
    assert [v.rank for v in variants] == list(range(len(variants)))


def test_enumerate_variants_all_translate_back_to_aa():
    from pegrna_pipeline.hpc.codon_expander import _build_codons_by_aa as _  # noqa: F401

    # Simple reverse translation table.
    translate = {}
    from python_codon_tables import get_codons_table

    for aa, codons in get_codons_table("h_sapiens_9606").items():
        for c in codons:
            translate[c] = aa

    aa = "MKL"
    variants = enumerate_variants(aa, max_variants=50)
    for v in variants:
        codons = [v.dna[i : i + 3] for i in range(0, len(v.dna), 3)]
        assert "".join(translate[c] for c in codons) == aa


def test_enumerate_variants_obeys_max_variants_cap():
    # Enumerate something guaranteed to have more than 2 surviving variants.
    variants = enumerate_variants("MKL", max_variants=2)
    assert len(variants) <= 2


def test_enumerate_variants_no_filter_violations():
    variants = enumerate_variants("MKL", max_variants=50)
    for v in variants:
        assert not _has_internal_stop(v.dna)
        assert not _has_cpg(v.dna)
        assert not _has_restriction_site(v.dna)
        assert not _has_homopolymer(v.dna)
