"""Construct grammar: tokenizer, raw-component parser, and resolver.

Grammar (informal):

    construct ::= component ("-" component)*
    component ::= atom ("x" NUMBER | NUMBER)?          # repeat suffix
                | "(" construct ")" ("x" NUMBER | NUMBER)   # parenthesised repeat
    atom      ::= NAME                                 # motif name
                | "aa:" <greedy protein string>
                | "dna:" <greedy DNA string>

Notes:
- NAMEs are `[A-Za-z][A-Za-z0-9_]*`. No hyphens (reserved as separator).
- 'x' is only the X token when:
    (a) the previous character emitted a NAME or RPAREN, AND
    (b) the next character is a digit.
  Otherwise 'x' is just a name/identifier character.
- Raw `aa:` / `dna:` payloads run greedily until the next '-' or end-of-input.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto


class TokenKind(Enum):
    NAME = auto()
    DASH = auto()
    LPAREN = auto()
    RPAREN = auto()
    X = auto()
    NUMBER = auto()
    COLON = auto()  # not emitted externally — absorbed into RAW_* payloads
    RAW_AA = auto()
    RAW_DNA = auto()


@dataclass(frozen=True)
class Token:
    kind: TokenKind
    value: str
    position: int


def _is_name_start(ch: str) -> bool:
    return ch.isalpha()


def _is_name_cont(ch: str) -> bool:
    return ch.isalnum() or ch == "_"


def tokenize(source: str) -> list[Token]:
    """Tokenise a construct string into a flat list of Tokens.

    Raises ValueError on unrecognised characters.
    """
    tokens: list[Token] = []
    i = 0
    n = len(source)

    while i < n:
        ch = source[i]

        # Whitespace is ignored.
        if ch.isspace():
            i += 1
            continue

        # Structural characters.
        if ch == "-":
            tokens.append(Token(TokenKind.DASH, "-", i))
            i += 1
            continue
        if ch == "(":
            tokens.append(Token(TokenKind.LPAREN, "(", i))
            i += 1
            continue
        if ch == ")":
            tokens.append(Token(TokenKind.RPAREN, ")", i))
            i += 1
            continue

        # Numeric repeat count.
        if ch.isdigit():
            start = i
            while i < n and source[i].isdigit():
                i += 1
            tokens.append(Token(TokenKind.NUMBER, source[start:i], start))
            continue

        # After an RPAREN, an 'x' followed by a digit is an X token, not the
        # start of a new NAME.
        if (
            ch == "x"
            and tokens
            and tokens[-1].kind is TokenKind.RPAREN
            and i + 1 < n
            and source[i + 1].isdigit()
        ):
            tokens.append(Token(TokenKind.X, "x", i))
            i += 1
            continue

        # NAME (possibly absorbed by later 'x' / payload logic).
        if _is_name_start(ch):
            start = i
            # Walk a maximal run of name characters.
            j = i + 1
            while j < n and _is_name_cont(source[j]):
                j += 1
            candidate = source[start:j]

            # Raw payload prefix?
            if (
                candidate in ("aa", "dna")
                and j < n
                and source[j] == ":"
            ):
                kind = TokenKind.RAW_AA if candidate == "aa" else TokenKind.RAW_DNA
                # Payload: greedy up to next '-' or end.
                payload_start = j + 1
                k = payload_start
                while k < n and source[k] != "-":
                    k += 1
                tokens.append(Token(kind, source[payload_start:k], start))
                i = k
                continue

            # Otherwise we need to walk the name character by character so we
            # can peel off a trailing 'x' when it's followed by a digit AND
            # there's at least one preceding character emitted as part of the
            # name. I.e. scan greedily, but if we see "x<digit>" at a position
            # where we've already consumed at least one name character, stop
            # the name there and let X + NUMBER take over.
            k = start
            while k < n:
                c = source[k]
                if not _is_name_cont(c):
                    break
                # Potential X boundary: lowercase 'x' followed by a digit, and
                # we've already consumed >= 1 name character.
                if c == "x" and k > start and k + 1 < n and source[k + 1].isdigit():
                    break
                k += 1
            tokens.append(Token(TokenKind.NAME, source[start:k], start))
            i = k

            # Immediately after a NAME, check for an X token.
            if i < n and source[i] == "x" and i + 1 < n and source[i + 1].isdigit():
                tokens.append(Token(TokenKind.X, "x", i))
                i += 1
            continue

        raise ValueError(
            f"unexpected character {ch!r} at position {i} in construct {source!r}"
        )

    return tokens
