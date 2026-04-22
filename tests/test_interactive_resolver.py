"""Tests for the interactive resolver for unknown motif names.

The resolver presents numbered choices (fuzzy suggestions first, then `aa:`,
`dna:`, and `abort`) and resolves the user's pick into a RawComponent or
raises an abort exception.
"""

import pytest

from pegrna_pipeline.agent.construct_parser import RawComponent, RawComponentKind
from pegrna_pipeline.agent.interactive_resolver import (
    AbortedByUser,
    build_interactive_resolver,
)
from pegrna_pipeline.agent.motif_library import MotifLibrary


@pytest.fixture
def library():
    return MotifLibrary.load_bundled()


def _scripted_input(responses):
    """Return a prompt callable that replays responses in order."""
    it = iter(responses)

    def _prompt(msg: str) -> str:
        try:
            return next(it)
        except StopIteration:  # pragma: no cover - shouldn't happen in these tests
            raise AssertionError(f"interactive resolver asked for more input than scripted: {msg}")

    return _prompt


def test_resolver_offers_fuzzy_suggestion_and_picks_it(library, capsys):
    resolver = build_interactive_resolver(library, prompt=_scripted_input(["1"]))
    # "FLAC" is 1 edit from "FLAG".
    result = resolver("FLAC")
    assert isinstance(result, RawComponent)
    assert result.kind is RawComponentKind.MOTIF
    assert result.value == "FLAG"
    out = capsys.readouterr().out
    assert "FLAC" in out
    assert "FLAG" in out


def test_resolver_picks_raw_aa(library):
    # Two suggestions at most — then aa, dna, abort. pick "aa" by its index.
    # Start with a query that has no suggestions.
    resolver = build_interactive_resolver(
        library,
        prompt=_scripted_input(["aa", "DYKDDDDK"]),
    )
    result = resolver("ZZZZZ")
    assert result.kind is RawComponentKind.RAW_AA
    assert result.value == "DYKDDDDK"


def test_resolver_picks_raw_dna(library):
    resolver = build_interactive_resolver(
        library,
        prompt=_scripted_input(["dna", "ATGCGA"]),
    )
    result = resolver("ZZZZZ")
    assert result.kind is RawComponentKind.RAW_DNA
    assert result.value == "ATGCGA"


def test_resolver_abort_raises(library):
    resolver = build_interactive_resolver(library, prompt=_scripted_input(["abort"]))
    with pytest.raises(AbortedByUser):
        resolver("ZZZZZ")


def test_resolver_retries_on_invalid_input(library):
    # First answer is invalid, second is "abort".
    resolver = build_interactive_resolver(
        library,
        prompt=_scripted_input(["not-a-choice", "abort"]),
    )
    with pytest.raises(AbortedByUser):
        resolver("ZZZZZ")


def test_resolver_numbered_choice_includes_aa_and_dna(library, capsys):
    resolver = build_interactive_resolver(library, prompt=_scripted_input(["abort"]))
    with pytest.raises(AbortedByUser):
        resolver("ZZZZZ")
    out = capsys.readouterr().out
    assert "aa:" in out.lower() or "aa " in out.lower()
    assert "dna:" in out.lower() or "dna " in out.lower()
    assert "abort" in out.lower()
