"""Guard against the pegRNA / ngRNA CSVs silently drifting in schema.

Both parsers must work on the fixture pair **simultaneously** — that's how
the real pipeline consumes them (pegRNA scores + nicking guides come out of
one PRIDICT run). This test is a regression trap: if either required
column set changes and someone forgets to update the constants, the
integration chain breaks here.
"""

from __future__ import annotations

from pathlib import Path

from pegrna_pipeline.hpc.pridict_parser import parse_ngrna_csv, parse_pegrna_csv

FIXTURES = Path(__file__).parent / "fixtures" / "pridict"
PEGRNA_CSV = FIXTURES / "ctcf_tjpe_v1_pegRNA_Pridict_full.csv"
NGRNA_CSV = FIXTURES / "ctcf_tjpe_v1_nicking_guides.csv"


def test_both_parsers_work_on_fixture_pair():
    pegrnas = parse_pegrna_csv(PEGRNA_CSV)
    ngrnas = parse_ngrna_csv(NGRNA_CSV)
    assert pegrnas, "pegRNA fixture parsed to empty list"
    assert ngrnas, "ngRNA fixture parsed to empty list"


def test_target_strand_column_names_differ_between_files():
    """Document the annoying hyphen-vs-underscore quirk in PRIDICT output."""
    import csv

    with PEGRNA_CSV.open(newline="") as fh:
        pegrna_headers = next(csv.reader(fh))
    with NGRNA_CSV.open(newline="") as fh:
        ngrna_headers = next(csv.reader(fh))

    # pegRNA uses hyphen; ngRNA uses underscore. This is a property of the
    # files themselves, not of our code — the parser has to cope with both.
    assert "Target-Strand" in pegrna_headers
    assert "Target_Strand" in ngrna_headers
    # And *not* the other way around.
    assert "Target_Strand" not in pegrna_headers
    assert "Target-Strand" not in ngrna_headers
