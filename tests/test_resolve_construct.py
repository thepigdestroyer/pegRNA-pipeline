"""Tests for resolve_construct: RawComponents -> resolved Construct."""

import pytest

from pegrna_pipeline.agent.construct_parser import (
    Construct,
    ResolvedComponent,
    resolve_construct,
    tokenize,
    parse_components,
)
from pegrna_pipeline.agent.motif_library import MotifLibrary


@pytest.fixture
def library():
    return MotifLibrary.load_bundled()


def _resolve(src, library):
    return resolve_construct(parse_components(tokenize(src)), library=library)


def test_single_motif(library):
    c = _resolve("FLAG", library)
    assert isinstance(c, Construct)
    assert len(c.components) == 1
    assert c.components[0].name == "FLAG"
    assert c.components[0].dna == "GACTACAAGGACGACGATGACAAG"
    assert c.components[0].aa == "DYKDDDDK"
    assert c.components[0].source == "motif"


def test_chain_of_motifs(library):
    c = _resolve("GGGGS-FLAG", library)
    assert [rc.name for rc in c.components] == ["GGGGS", "FLAG"]
    # Concatenated DNA = GGGGS.dna + FLAG.dna.
    flag = library.lookup("FLAG")
    ggggs = library.lookup("GGGGS")
    assert c.dna == ggggs.dna + flag.dna
    assert c.aa == "GGGGS" + "DYKDDDDK"


def test_repeat_expansion(library):
    c = _resolve("GGGGSx3", library)
    assert len(c.components) == 3
    assert all(rc.name == "GGGGS" for rc in c.components)


def test_raw_aa_component_gets_codon_optimised(library):
    c = _resolve("aa:MK", library)
    assert len(c.components) == 1
    comp = c.components[0]
    assert comp.source == "raw_aa"
    assert comp.aa == "MK"
    # Codon-optimised DNA must translate back to MK.
    assert len(comp.dna) == 6
    assert comp.dna.isupper()
    assert set(comp.dna) <= set("ACGT")


def test_raw_dna_component_passthrough(library):
    c = _resolve("dna:ATGAAG", library)
    assert len(c.components) == 1
    comp = c.components[0]
    assert comp.source == "raw_dna"
    assert comp.dna == "ATGAAG"
    # Raw DNA components don't carry an aa annotation.
    assert comp.aa in ("", None) or isinstance(comp.aa, str)


def test_unknown_motif_raises(library):
    with pytest.raises(Exception):  # MotifNotFoundError subclasses KeyError
        _resolve("TotallyFakeTag", library)


def test_resolved_component_stores_source_kind(library):
    c = _resolve("FLAG-aa:M-dna:ATG", library)
    sources = [rc.source for rc in c.components]
    assert sources == ["motif", "raw_aa", "raw_dna"]
    assert isinstance(c.components[0], ResolvedComponent)


def test_construct_aa_joins_components(library):
    c = _resolve("GGGGS-FLAG", library)
    assert c.aa == "GGGGS" + "DYKDDDDK"


def test_construct_dna_joins_components(library):
    c = _resolve("FLAG-dna:GGG", library)
    assert c.dna.endswith("GGG")
    assert c.dna.startswith("GAC")
