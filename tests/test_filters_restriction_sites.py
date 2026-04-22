"""Tests for restriction site filter.

The site table lives at ``data/restriction_sites.yaml`` and must include the
twelve commonly-used enzymes below. _has_restriction_site checks for any of
the listed recognition sequences as a substring of the given DNA.
"""

from pathlib import Path

import pytest
import yaml

from pegrna_pipeline.hpc.codon_expander import (
    _has_restriction_site,
    load_restriction_sites,
)

YAML_PATH = (
    Path(__file__).resolve().parent.parent
    / "src"
    / "pegrna_pipeline"
    / "data"
    / "restriction_sites.yaml"
)


EXPECTED_ENZYMES = (
    "BamHI",
    "EcoRI",
    "HindIII",
    "XhoI",
    "NotI",
    "SalI",
    "BglII",
    "NdeI",
    "KpnI",
    "SacI",
    "SpeI",
    "XbaI",
)


def test_restriction_sites_yaml_exists_and_has_required_enzymes():
    assert YAML_PATH.is_file()
    data = yaml.safe_load(YAML_PATH.read_text())
    present = {e["name"] for e in data["enzymes"]}
    for enzyme in EXPECTED_ENZYMES:
        assert enzyme in present, f"missing {enzyme}"


def test_recognition_sequences_are_all_acgt():
    data = yaml.safe_load(YAML_PATH.read_text())
    for e in data["enzymes"]:
        assert set(e["recognition"].upper()) <= set("ACGT"), e


def test_load_restriction_sites_returns_list_of_strings():
    sites = load_restriction_sites()
    assert isinstance(sites, list)
    assert sites
    assert all(isinstance(s, str) and len(s) >= 4 for s in sites)


@pytest.mark.parametrize(
    "enzyme,dna",
    [
        ("BamHI", "AAAGGATCCTTT"),
        ("EcoRI", "ACGGAATTCACG"),
        ("HindIII", "GGGAAGCTTCCC"),
        ("XhoI", "ACCCTCGAGGTT"),
        ("NotI", "AATGCGGCCGCAA"),
        ("SalI", "TTTGTCGACGGG"),
    ],
)
def test_has_restriction_site_detects_known_sites(enzyme, dna):
    assert _has_restriction_site(dna) is True, f"{enzyme} site not detected in {dna}"


def test_clean_dna_has_no_site():
    assert _has_restriction_site("ATGAAAGGGCCCTTT") is False


def test_empty_sequence():
    assert _has_restriction_site("") is False
