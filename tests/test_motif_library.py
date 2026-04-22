"""Tests for MotifLibrary class."""

import pytest

from pegrna_pipeline.agent.motif_library import (
    Motif,
    MotifLibrary,
    MotifNotFoundError,
)


def test_load_bundled_returns_library_with_entries():
    lib = MotifLibrary.load_bundled()
    assert isinstance(lib, MotifLibrary)
    assert len(lib.list_names()) >= 30


def test_lookup_returns_motif():
    lib = MotifLibrary.load_bundled()
    flag = lib.lookup("FLAG")
    assert isinstance(flag, Motif)
    assert flag.name == "FLAG"
    assert flag.dna == "GACTACAAGGACGACGATGACAAG"
    assert flag.aa == "DYKDDDDK"


def test_lookup_is_case_insensitive():
    lib = MotifLibrary.load_bundled()
    assert lib.lookup("flag").name == "FLAG"
    assert lib.lookup("Flag").name == "FLAG"
    assert lib.lookup("FLAG").name == "FLAG"
    # An alias differing only in casing (ZEOCINr) preserves its own casing.
    assert lib.lookup("zeocinr").name in {"ZeocinR", "ZEOCINr"}


def test_lookup_unknown_raises_motif_not_found():
    lib = MotifLibrary.load_bundled()
    with pytest.raises(MotifNotFoundError) as exc:
        lib.lookup("DoesNotExist")
    assert "DoesNotExist" in str(exc.value)


def test_list_names_returns_all_motif_names_sorted_case_insensitive():
    lib = MotifLibrary.load_bundled()
    names = lib.list_names()
    assert "FLAG" in names
    assert "GGGGS" in names
    # Sorted case-insensitively.
    assert names == sorted(names, key=str.lower)
    # Library intentionally includes the ZEOCINr casing alias alongside ZeocinR.
    assert "ZeocinR" in names
    assert "ZEOCINr" in names
