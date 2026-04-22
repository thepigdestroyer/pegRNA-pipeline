"""End-to-end integration tests for construct resolution.

Builds non-trivial constructs from the bundled library and verifies that the
resolved aa / dna exactly match the hand-concatenation of motif DNAs / aas.
"""

import pytest

from pegrna_pipeline.agent.construct_parser import (
    parse_components,
    resolve_construct,
    tokenize,
)
from pegrna_pipeline.agent.motif_library import MotifLibrary


@pytest.fixture
def lib():
    return MotifLibrary.load_bundled()


def _resolve(src, lib):
    return resolve_construct(parse_components(tokenize(src)), library=lib)


def test_real_world_c_terminal_flag_linker_puro_construct(lib):
    # A real C-terminal tagging layout: 5-GGGGS linker, FLAG tag, P2A, puromycin.
    construct = _resolve("GGGGSx2-FLAG-P2A-PuroR", lib)
    parts = ["GGGGS", "GGGGS", "FLAG", "P2A", "PuroR"]
    assert [c.name for c in construct.components] == parts
    # AA concatenation matches.
    assert construct.aa == "".join(lib.lookup(p).aa for p in parts)
    # DNA concatenation matches.
    assert construct.dna == "".join(lib.lookup(p).dna for p in parts)


def test_mixed_motifs_and_raw_aa_payload_round_trip(lib):
    construct = _resolve("FLAG-aa:GSGSG-HA", lib)
    # Component kinds are motif / raw_aa / motif.
    sources = [c.source for c in construct.components]
    assert sources == ["motif", "raw_aa", "motif"]
    # AA sequence interleaves FLAG, GSGSG, HA.
    assert construct.aa == "DYKDDDDK" + "GSGSG" + "YPYDVPDYA"
    # Raw aa payload DNA encodes GSGSG (15 nt).
    raw_dna = construct.components[1].dna
    assert len(raw_dna) == 15
    assert set(raw_dna) <= set("ACGT")


def test_nested_repeat_group(lib):
    construct = _resolve("(GGGGS-FLAG)2", lib)
    names = [c.name for c in construct.components]
    assert names == ["GGGGS", "FLAG", "GGGGS", "FLAG"]


def test_kozak_regulatory_included_in_dna_but_not_aa(lib):
    construct = _resolve("Kozak-FLAG", lib)
    # Kozak adds DNA but no aa.
    assert construct.aa == "DYKDDDDK"
    assert construct.dna.startswith(lib.lookup("Kozak").dna)


def test_stop_codon_included(lib):
    construct = _resolve("FLAG-STOP_TAA", lib)
    assert construct.dna.endswith("TAA")


def test_user_override_is_merged_via_load_full(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    import yaml as _yaml
    (tmp_path / ".pegrna").mkdir()
    (tmp_path / ".pegrna" / "motifs.yaml").write_text(
        _yaml.safe_dump({"motifs": [{
            "name": "LabTag", "category": "custom",
            "dna": "ATGAAG", "aa": "MK", "description": ""
        }]})
    )
    lib = MotifLibrary.load_full()
    construct = resolve_construct(
        parse_components(tokenize("FLAG-LabTag")),
        library=lib,
    )
    assert construct.aa.endswith("MK")
    assert construct.dna.endswith("ATGAAG")
