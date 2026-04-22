"""Tokenizer tests for the construct grammar.

Grammar summary:
- NAME:    [A-Za-z][A-Za-z0-9_]*   (NO hyphens; underscores fine)
- DASH:    -                       (component separator)
- LPAREN:  (
- RPAREN:  )
- X:       x                       (repeat operator, only between a NAME/RPAREN and a digit)
- NUMBER:  digits                  (repeat count; follows X or paren group)
- COLON:   :
- RAW:     after "aa:" or "dna:", greedy up to next DASH/EOL
"""

import pytest

from pegrna_pipeline.agent.construct_parser import (
    Token,
    TokenKind,
    tokenize,
)


def _kinds(tokens):
    return [t.kind for t in tokens]


def test_single_name():
    toks = tokenize("FLAG")
    assert _kinds(toks) == [TokenKind.NAME]
    assert toks[0].value == "FLAG"


def test_dashed_chain_of_names():
    toks = tokenize("GGGGS-FLAG-P2A")
    assert _kinds(toks) == [
        TokenKind.NAME, TokenKind.DASH,
        TokenKind.NAME, TokenKind.DASH,
        TokenKind.NAME,
    ]
    assert [t.value for t in toks] == ["GGGGS", "-", "FLAG", "-", "P2A"]


def test_repeat_with_explicit_x_and_digit():
    toks = tokenize("GGGGSx3")
    assert _kinds(toks) == [TokenKind.NAME, TokenKind.X, TokenKind.NUMBER]
    assert toks[-1].value == "3"


def test_name_ending_in_x_still_a_name():
    # Single-letter M followed by x1 -> Mx1 means M repeated once.
    toks = tokenize("Mx1")
    assert _kinds(toks) == [TokenKind.NAME, TokenKind.X, TokenKind.NUMBER]
    assert toks[0].value == "M"


def test_name_with_internal_x_not_followed_by_digit_stays_in_name():
    # Mxx -> single NAME "Mxx" because trailing 'x' has no digit after.
    toks = tokenize("Mxx")
    assert _kinds(toks) == [TokenKind.NAME]
    assert toks[0].value == "Mxx"


def test_x_inside_name_between_letters_is_name_char():
    # Actually 'E' X 'E' X 1 -> (NAME "ExE") (X) (NUMBER "1")? No:
    # the grammar says "x" is a repeat only when preceded by NAME/RPAREN AND
    # followed by a digit. So in "ExEx1": first 'x' is followed by 'E' (not a
    # digit), so it's part of the name; second 'x' is followed by '1', so it
    # splits off as X. Result: NAME("ExE"), X, NUMBER("1").
    toks = tokenize("ExEx1")
    assert _kinds(toks) == [TokenKind.NAME, TokenKind.X, TokenKind.NUMBER]
    assert toks[0].value == "ExE"
    assert toks[-1].value == "1"


def test_parenthesised_group_repeat():
    toks = tokenize("(GGGGS)3")
    assert _kinds(toks) == [
        TokenKind.LPAREN, TokenKind.NAME, TokenKind.RPAREN, TokenKind.NUMBER
    ]


def test_parenthesised_group_with_explicit_x():
    toks = tokenize("(GGGGS)x3")
    assert _kinds(toks) == [
        TokenKind.LPAREN, TokenKind.NAME, TokenKind.RPAREN,
        TokenKind.X, TokenKind.NUMBER,
    ]


def test_raw_aa_payload_greedy_to_dash():
    toks = tokenize("aa:DYKDDDDK-FLAG")
    assert _kinds(toks) == [
        TokenKind.RAW_AA, TokenKind.DASH, TokenKind.NAME
    ]
    assert toks[0].value == "DYKDDDDK"


def test_raw_aa_payload_to_end_of_input():
    toks = tokenize("aa:DYKDDDDK")
    assert _kinds(toks) == [TokenKind.RAW_AA]
    assert toks[0].value == "DYKDDDDK"


def test_raw_dna_payload():
    toks = tokenize("dna:ATGCGAT")
    assert _kinds(toks) == [TokenKind.RAW_DNA]
    assert toks[0].value == "ATGCGAT"


def test_name_with_underscore_ok():
    toks = tokenize("NLS_SV40-bNLS")
    assert _kinds(toks) == [TokenKind.NAME, TokenKind.DASH, TokenKind.NAME]
    assert toks[0].value == "NLS_SV40"


def test_multiple_features_combine():
    toks = tokenize("(GGGGS)3-FLAGx2-aa:DYK-ZeocinR")
    kinds = _kinds(toks)
    # LPAREN NAME RPAREN NUMBER DASH NAME X NUMBER DASH RAW_AA DASH NAME
    assert kinds == [
        TokenKind.LPAREN, TokenKind.NAME, TokenKind.RPAREN, TokenKind.NUMBER,
        TokenKind.DASH, TokenKind.NAME, TokenKind.X, TokenKind.NUMBER,
        TokenKind.DASH, TokenKind.RAW_AA, TokenKind.DASH, TokenKind.NAME,
    ]


def test_tokens_are_dataclasses_with_position():
    toks = tokenize("FLAG")
    assert isinstance(toks[0], Token)
    assert hasattr(toks[0], "position")
    assert toks[0].position == 0


def test_unknown_char_raises():
    with pytest.raises(ValueError):
        tokenize("FLAG@HA")
