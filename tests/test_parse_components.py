"""Tests for parse_components: token stream -> list of RawComponent."""

import pytest

from pegrna_pipeline.agent.construct_parser import (
    RawComponent,
    RawComponentKind,
    parse_components,
    tokenize,
)


def _parse(src):
    return parse_components(tokenize(src))


def test_single_name():
    components = _parse("FLAG")
    assert components == [RawComponent(RawComponentKind.MOTIF, "FLAG")]


def test_dashed_chain_preserves_order():
    components = _parse("GGGGS-FLAG-P2A")
    assert [c.value for c in components] == ["GGGGS", "FLAG", "P2A"]


def test_repeat_with_x_expands_inline():
    components = _parse("GGGGSx3")
    assert components == [
        RawComponent(RawComponentKind.MOTIF, "GGGGS"),
        RawComponent(RawComponentKind.MOTIF, "GGGGS"),
        RawComponent(RawComponentKind.MOTIF, "GGGGS"),
    ]


def test_repeat_without_x_in_parenthesised_group():
    components = _parse("(GGGGS)3")
    assert [c.value for c in components] == ["GGGGS", "GGGGS", "GGGGS"]


def test_repeat_with_x_in_parenthesised_group():
    components = _parse("(GGGGS)x3")
    assert [c.value for c in components] == ["GGGGS", "GGGGS", "GGGGS"]


def test_multi_motif_group_repeat():
    components = _parse("(GGGGS-FLAG)2")
    assert [c.value for c in components] == [
        "GGGGS",
        "FLAG",
        "GGGGS",
        "FLAG",
    ]


def test_raw_aa_component_preserved():
    components = _parse("FLAG-aa:DYKDDDDK-GGGGS")
    assert components == [
        RawComponent(RawComponentKind.MOTIF, "FLAG"),
        RawComponent(RawComponentKind.RAW_AA, "DYKDDDDK"),
        RawComponent(RawComponentKind.MOTIF, "GGGGS"),
    ]


def test_raw_dna_component_preserved():
    components = _parse("dna:ATGCGAT-FLAG")
    assert components == [
        RawComponent(RawComponentKind.RAW_DNA, "ATGCGAT"),
        RawComponent(RawComponentKind.MOTIF, "FLAG"),
    ]


def test_zero_repeat_rejected():
    with pytest.raises(ValueError):
        _parse("GGGGSx0")


def test_negative_implicit_empty_repeat():
    # e.g. (FLAG) with no following digit should just yield ["FLAG"].
    components = _parse("(FLAG)")
    assert [c.value for c in components] == ["FLAG"]


def test_dangling_x_followed_by_eof_is_name_not_operator():
    # "FLAGx" has no digit after 'x', so 'x' is absorbed into the NAME.
    components = _parse("FLAGx")
    assert components == [RawComponent(RawComponentKind.MOTIF, "FLAGx")]


def test_repeat_of_a_repeat():
    # (FLAGx2)2 -> FLAG FLAG FLAG FLAG
    components = _parse("(FLAGx2)2")
    assert [c.value for c in components] == ["FLAG"] * 4


def test_trailing_dash_rejected():
    with pytest.raises(ValueError):
        _parse("FLAG-")
