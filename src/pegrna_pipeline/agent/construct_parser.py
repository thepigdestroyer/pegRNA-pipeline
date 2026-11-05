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
from typing import Iterable, Sequence

from pegrna_pipeline.hpc.cargo_spec import library_part_is_coding


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
            if candidate in ("aa", "dna") and j < n and source[j] == ":":
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

        raise ValueError(f"unexpected character {ch!r} at position {i} in construct {source!r}")

    return tokens


# ---------------------------------------------------------------- RawComponent


class RawComponentKind(Enum):
    MOTIF = auto()
    RAW_AA = auto()
    RAW_DNA = auto()


@dataclass(frozen=True)
class RawComponent:
    kind: RawComponentKind
    value: str


def parse_components(tokens: Sequence[Token]) -> list[RawComponent]:
    """Convert a token stream into a flat list of RawComponent entries.

    Expands x<N> and (group)<N> repeats inline. Raises ValueError on any
    syntactic violation (dangling -, unmatched (, repeat count of zero, etc.).
    """
    components, consumed = _parse_sequence(tokens, 0, terminator=None)
    if consumed != len(tokens):
        raise ValueError(f"unexpected trailing tokens at position {consumed}")
    return components


def _parse_sequence(
    tokens: Sequence[Token], start: int, *, terminator: TokenKind | None
) -> tuple[list[RawComponent], int]:
    """Parse component (DASH component)* stopping at `terminator` or EOF."""
    out: list[RawComponent] = []
    i = start

    # Empty sequences are not valid constructs.
    if i >= len(tokens) or tokens[i].kind is terminator:
        raise ValueError("empty construct/group")

    while True:
        chunk, i = _parse_component(tokens, i)
        out.extend(chunk)

        if i >= len(tokens):
            if terminator is not None:
                raise ValueError("unterminated parenthesised group")
            return out, i

        if tokens[i].kind is terminator:
            return out, i

        if tokens[i].kind is TokenKind.DASH:
            i += 1
            if i >= len(tokens) or tokens[i].kind is terminator:
                raise ValueError(f"trailing '-' at position {tokens[i - 1].position}")
            continue

        raise ValueError(
            f"expected '-' or end of construct at token position {tokens[i].position}, "
            f"got {tokens[i].kind.name}"
        )


def _parse_component(tokens: Sequence[Token], i: int) -> tuple[list[RawComponent], int]:
    """Parse a single component (plus any repeat suffix). Returns the expanded
    list of RawComponents and the new cursor.
    """
    tok = tokens[i]

    # Parenthesised group.
    if tok.kind is TokenKind.LPAREN:
        inner, j = _parse_sequence(tokens, i + 1, terminator=TokenKind.RPAREN)
        # Consume RPAREN.
        j += 1
        # Optional repeat suffix: either X NUMBER or bare NUMBER.
        count, j = _consume_repeat(tokens, j, allow_bare_number=True)
        return inner * count, j

    # Raw aa:/dna: payload.
    if tok.kind is TokenKind.RAW_AA:
        return [RawComponent(RawComponentKind.RAW_AA, tok.value)], i + 1
    if tok.kind is TokenKind.RAW_DNA:
        return [RawComponent(RawComponentKind.RAW_DNA, tok.value)], i + 1

    # Motif NAME with optional Xn suffix (bare numbers aren't allowed after a
    # NAME; they must come through explicit 'x').
    if tok.kind is TokenKind.NAME:
        count, j = _consume_repeat(tokens, i + 1, allow_bare_number=False)
        return [RawComponent(RawComponentKind.MOTIF, tok.value)] * count, j

    raise ValueError(f"unexpected token {tok.kind.name} at position {tok.position}")


def _consume_repeat(tokens: Sequence[Token], i: int, *, allow_bare_number: bool) -> tuple[int, int]:
    """If the next token(s) form a repeat suffix, consume them and return the
    count plus the new cursor. Otherwise return (1, i) unchanged.
    """
    if i >= len(tokens):
        return 1, i

    if tokens[i].kind is TokenKind.X:
        if i + 1 >= len(tokens) or tokens[i + 1].kind is not TokenKind.NUMBER:
            raise ValueError(f"'x' at position {tokens[i].position} is not followed by a number")
        count = int(tokens[i + 1].value)
        if count <= 0:
            raise ValueError(f"repeat count must be > 0, got {count}")
        return count, i + 2

    if allow_bare_number and tokens[i].kind is TokenKind.NUMBER:
        count = int(tokens[i].value)
        if count <= 0:
            raise ValueError(f"repeat count must be > 0, got {count}")
        return count, i + 1

    return 1, i


def _ignored() -> Iterable[None]:  # pragma: no cover - typing anchor
    yield None


# ---------------------------------------------------------- Resolved construct


@dataclass(frozen=True)
class ResolvedComponent:
    """A RawComponent after motif lookup / codon optimisation."""

    name: str  # motif name, or a synthetic label for raw payloads
    dna: str  # uppercase ACGT
    aa: str  # single-letter amino acid sequence ("" for raw_dna / Kozak)
    source: str  # "motif" | "raw_aa" | "raw_dna"
    #: Never recoded (ADR 0029): not protein (a promoter, a 5'UTR, a stop codon), or raw DNA typed in as is.
    locked: bool = False


@dataclass(frozen=True)
class Construct:
    components: tuple[ResolvedComponent, ...]

    @property
    def dna(self) -> str:
        return "".join(c.dna for c in self.components)

    @property
    def aa(self) -> str:
        return "".join(c.aa for c in self.components)

    def cargo_parts(self) -> list[dict]:
        """The ``cargo_parts`` a design config takes: each part's DNA and whether it is locked."""
        return [{"name": c.name, "dna": c.dna, "locked": c.locked} for c in self.components]


def _codon_optimise_simple(aa_sequence: str) -> str:
    """Codon-optimise a raw aa: payload via :func:`optimise_bulk`.

    Lazy import so the agent layer doesn't pull in hpc internals at import
    time (and so tests can monkeypatch the expander module cleanly).
    """
    from pegrna_pipeline.hpc.codon_expander import optimise_bulk

    return optimise_bulk(aa_sequence)


def resolve_construct(
    components: Iterable[RawComponent],
    *,
    library: "object",
) -> Construct:
    """Resolve raw components to a Construct.

    The ``library`` argument is a :class:`MotifLibrary` — typed loosely to
    avoid an import cycle.
    """
    resolved: list[ResolvedComponent] = []
    for raw in components:
        if raw.kind is RawComponentKind.MOTIF:
            motif = library.lookup(raw.value)
            resolved.append(
                ResolvedComponent(
                    name=motif.name,
                    dna=motif.dna.upper(),
                    aa=motif.aa,
                    source="motif",
                    locked=not library_part_is_coding(aa=motif.aa, category=motif.category),
                )
            )
        elif raw.kind is RawComponentKind.RAW_AA:
            aa = raw.value.upper()
            dna = _codon_optimise_simple(aa)
            resolved.append(
                ResolvedComponent(
                    name=f"aa:{aa}",
                    dna=dna,
                    aa=aa,
                    source="raw_aa",
                )
            )
        elif raw.kind is RawComponentKind.RAW_DNA:
            dna = raw.value.upper()
            if set(dna) - set("ACGT"):
                raise ValueError(f"raw DNA payload contains non-ACGT: {raw.value!r}")
            resolved.append(
                ResolvedComponent(
                    name=f"dna:{dna}",
                    dna=dna,
                    aa="",
                    source="raw_dna",
                    locked=True,  # DNA typed in as is is kept as is
                )
            )
        else:  # pragma: no cover - defensive
            raise ValueError(f"unhandled RawComponentKind {raw.kind!r}")
    return Construct(tuple(resolved))
