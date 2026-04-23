"""Tests for :mod:`pegrna_pipeline.hpc.pridict_parser`.

Two parsers live in the module:

* :func:`parse_pegrna_csv` reads ``*_pegRNA_Pridict_full.csv``.
* :func:`parse_ngrna_csv` reads ``*_nicking_guides.csv``.

Both raise :class:`ValueError` when a required column is missing, so a
caller never silently processes partial data.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from pegrna_pipeline.hpc.data_types import NgRNADesign, PegRNADesign
from pegrna_pipeline.hpc.pridict_parser import (
    _NGRNA_REQUIRED,
    _PEGRNA_REQUIRED,
    parse_ngrna_csv,
    parse_pegrna_csv,
)

FIXTURES = Path(__file__).parent / "fixtures" / "pridict"
PEGRNA_CSV = FIXTURES / "ctcf_tjpe_v1_pegRNA_Pridict_full.csv"
NGRNA_CSV = FIXTURES / "ctcf_tjpe_v1_nicking_guides.csv"


# ----------------------------------------------------------- pegRNA parser


def test_parse_pegrna_csv_returns_designs():
    designs = parse_pegrna_csv(PEGRNA_CSV)
    assert len(designs) == 3
    for d in designs:
        assert isinstance(d, PegRNADesign)
        assert d.target_strand in ("Fw", "Rv")
        assert 0.0 <= d.pe2_hek293t <= 1.0
        assert d.editing_position == 50  # see test_pridict_index_convention


def test_parse_pegrna_csv_preserves_field_values():
    designs = parse_pegrna_csv(PEGRNA_CSV)
    by_name = {d.sequence_name: d for d in designs}
    d1 = by_name["ctcf_tjpe_v1_pegRNA1"]
    assert d1.spacer == "GAATTCGTGTTCCTAAGCGT"
    assert d1.pbs == "TCGACGCTTAGGAACAC"
    assert d1.pbs_length == 13
    assert d1.rt_length == 30
    assert d1.target_strand == "Fw"
    assert d1.pe2_hek293t == pytest.approx(0.682)


def test_parse_pegrna_csv_missing_column_raises(tmp_path):
    """Drop one required column and confirm the parser errors."""
    dst = tmp_path / "bad_pegrna.csv"
    with PEGRNA_CSV.open(newline="") as fh:
        reader = csv.reader(fh)
        rows = list(reader)
    header = rows[0]
    # Drop 'Spacer-Sequence' column
    drop_idx = header.index("Spacer-Sequence")
    trimmed = [[c for i, c in enumerate(r) if i != drop_idx] for r in rows]
    with dst.open("w", newline="") as fh:
        csv.writer(fh).writerows(trimmed)
    with pytest.raises(ValueError, match="Spacer-Sequence"):
        parse_pegrna_csv(dst)


def test_pegrna_required_constant_is_the_full_set():
    """Guard against accidental deletions from the required-column list."""
    expected = {
        "sequence_name",
        "Spacer-Sequence",
        "PBSrevcomp",
        "RTrevcomp",
        "PBSlength",
        "RTlength",
        "Target-Strand",
        "Editing_Position",
        "PRIDICT2_0_editing_Score_deep_HEK",
    }
    assert set(_PEGRNA_REQUIRED) == expected


# ----------------------------------------------------------- ngRNA parser


def test_parse_ngrna_csv_returns_designs():
    designs = parse_ngrna_csv(NGRNA_CSV)
    assert len(designs) == 7
    for d in designs:
        assert isinstance(d, NgRNADesign)
        assert d.target_strand in ("Fw", "Rv")
        assert isinstance(d.distance_to_edit, int)


def test_parse_ngrna_csv_preserves_signed_distance():
    """A row with distance -30 must round-trip as an int(-30)."""
    designs = parse_ngrna_csv(NGRNA_CSV)
    distances = [d.distance_to_edit for d in designs]
    assert -30 in distances  # preserved as signed int, not abs()
    assert 55 in distances


def test_parse_ngrna_csv_missing_column_raises(tmp_path):
    dst = tmp_path / "bad_ngrna.csv"
    with NGRNA_CSV.open(newline="") as fh:
        rows = list(csv.reader(fh))
    header = rows[0]
    drop_idx = header.index("DeepCas9score")
    trimmed = [[c for i, c in enumerate(r) if i != drop_idx] for r in rows]
    with dst.open("w", newline="") as fh:
        csv.writer(fh).writerows(trimmed)
    with pytest.raises(ValueError, match="DeepCas9score"):
        parse_ngrna_csv(dst)


def test_ngrna_required_constant_is_the_full_set():
    """Note the *underscore* in Target_Strand — hyphen would be wrong here."""
    expected = {
        "Nicking-Protospacer",
        "Nicking-Position-to-edit",
        "Target_Strand",
        "DeepCas9score",
    }
    assert set(_NGRNA_REQUIRED) == expected
