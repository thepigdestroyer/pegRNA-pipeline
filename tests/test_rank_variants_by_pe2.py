"""Tests for :func:`rank_variants_by_pe2`.

Sorts a list of :class:`VariantScore` objects by ``top_score`` descending,
stably (so equal scores preserve input order).
"""

from __future__ import annotations

from pegrna_pipeline.hpc.data_types import PegRNADesign, VariantScore
from pegrna_pipeline.hpc.tjpe_workflow import rank_variants_by_pe2


def _peg(score: float, name: str = "peg") -> PegRNADesign:
    return PegRNADesign(
        sequence_name=name,
        spacer="A" * 20,
        pbs="C" * 13,
        rt="G" * 30,
        pbs_length=13,
        rt_length=30,
        target_strand="Fw",
        editing_position=50,
        pe2_hek293t=score,
        pe2_k562=score,
        pegrna_full="A" * 140,
    )


def _variant(name: str, top: float, extras: tuple = ()) -> VariantScore:
    designs = [_peg(top, f"{name}_top")]
    designs.extend(_peg(s, f"{name}_extra_{i}") for i, s in enumerate(extras))
    return VariantScore(
        sequence_name=name, codon_variant_dna="ACG", pegRNA_designs=designs
    )


def test_rank_variants_sorts_descending_by_top_score():
    a = _variant("a", 0.1)
    b = _variant("b", 0.8)
    c = _variant("c", 0.5)
    out = rank_variants_by_pe2([a, b, c])
    assert [v.sequence_name for v in out] == ["b", "c", "a"]


def test_rank_variants_is_stable_on_ties():
    a = _variant("a", 0.7)
    b = _variant("b", 0.7)
    c = _variant("c", 0.7)
    # Same top score; insertion order must be preserved.
    assert rank_variants_by_pe2([a, b, c]) == [a, b, c]
    assert rank_variants_by_pe2([c, a, b]) == [c, a, b]


def test_rank_variants_handles_empty_designs():
    empty = VariantScore(sequence_name="e", codon_variant_dna="ACG")
    a = _variant("a", 0.3)
    out = rank_variants_by_pe2([empty, a])
    # Empty variant has top_score == 0.0 so should end up last.
    assert out[0].sequence_name == "a"
    assert out[1].sequence_name == "e"


def test_rank_variants_empty_input():
    assert rank_variants_by_pe2([]) == []


def test_rank_variants_returns_a_new_list():
    """Callers should be able to mutate the returned list without affecting input."""
    a = _variant("a", 0.1)
    b = _variant("b", 0.8)
    original = [a, b]
    out = rank_variants_by_pe2(original)
    assert out is not original
