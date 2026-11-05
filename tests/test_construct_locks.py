"""A construct string resolves to parts that say whether they may be recoded (ADR 0029)."""

import json

from click.testing import CliRunner

from pegrna_pipeline.agent.cli import cli
from pegrna_pipeline.agent.construct_parser import parse_components, resolve_construct, tokenize
from pegrna_pipeline.agent.motif_library import MotifLibrary
from pegrna_pipeline.hpc.cargo_spec import build_cargo_spec, free_codons


def _resolve(text):
    return resolve_construct(parse_components(tokenize(text)), library=MotifLibrary.load_full())


def test_protein_parts_are_coding_and_a_5utr_a_stop_and_raw_dna_are_locked():
    c = {p.name: p.locked for p in _resolve("FLAG-Kozak-dna:ACGTAC-aa:MKP").components}
    assert c == {"FLAG": False, "Kozak": True, "dna:ACGTAC": True, "aa:MKP": False}


def test_a_stop_motif_is_locked():
    stops = [m.name for m in MotifLibrary.load_full() if m.category == "stop"]
    assert stops
    assert all(_resolve(name).components[0].locked for name in stops)


def test_cargo_parts_round_trip_into_the_lock_model():
    construct = _resolve("FLAG-GGGGS-Kozak-aa:MKP")
    parts = construct.cargo_parts()
    assert [p["name"] for p in parts] == ["FLAG", "GGGGS", "Kozak", "aa:MKP"]
    spec = build_cargo_spec(construct.dna, parts=parts)
    assert [p.locked for p in spec.parts] == [False, False, True, False]
    assert len(free_codons(spec)) == 8 + 5 + 3  # FLAG, GGGGS and MKP; the Kozak is not touched


def test_the_validate_command_marks_locked_parts_and_prints_the_parts_for_a_config():
    out = CliRunner().invoke(cli, ["validate", "FLAG-Kozak"]).output
    assert "coding aa=DYKDDDDK" in out and "LOCKED aa=-" in out
    parts = json.loads(out.strip().splitlines()[-1])
    assert parts == [
        {"name": "FLAG", "dna": "GACTACAAGGACGACGATGACAAG", "locked": False},
        {"name": "Kozak", "dna": "GCCACC", "locked": True},
    ]
