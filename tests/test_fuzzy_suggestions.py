"""Tests for Levenshtein-based motif name suggestions."""

from pegrna_pipeline.agent.motif_library import MotifLibrary, fuzzy_suggestions


def test_exact_match_unlikely_suggested_because_not_mistyped():
    # If the user typed exactly a valid motif, the query wouldn't have failed.
    lib = MotifLibrary.load_bundled()
    suggestions = fuzzy_suggestions("FLAG", lib)
    assert "FLAG" in suggestions


def test_one_character_typo_suggested():
    lib = MotifLibrary.load_bundled()
    suggestions = fuzzy_suggestions("FLAC", lib)
    assert "FLAG" in suggestions


def test_two_character_typo_suggested():
    lib = MotifLibrary.load_bundled()
    suggestions = fuzzy_suggestions("FLCK", lib)
    assert "FLAG" in suggestions


def test_very_dissimilar_name_returns_no_suggestions():
    lib = MotifLibrary.load_bundled()
    suggestions = fuzzy_suggestions("ZZZZZZZZZZZZ", lib)
    assert suggestions == []


def test_suggestions_are_case_insensitive():
    lib = MotifLibrary.load_bundled()
    # lowercase "flac" is still close to "FLAG".
    suggestions = fuzzy_suggestions("flac", lib)
    assert "FLAG" in suggestions


def test_suggestions_return_unique_names():
    lib = MotifLibrary.load_bundled()
    suggestions = fuzzy_suggestions("GGGGS", lib)
    assert len(suggestions) == len(set(suggestions))
