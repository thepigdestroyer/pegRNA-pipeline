"""Verify bundled motif_library.yaml integrity.

Every motif must have:
- name (no hyphens, since '-' is the grammar separator)
- category
- dna (ACGT only)
- aa (single-letter protein code, with '*' allowed for stop codons)
- description
And the DNA must translate (without its trailing stop, if any) to the listed aa.
"""

from pathlib import Path

import pytest
import yaml

LIBRARY = (
    Path(__file__).resolve().parent.parent
    / "src"
    / "pegrna_pipeline"
    / "data"
    / "motif_library.yaml"
)

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _translate(dna: str) -> str:
    aa = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i : i + 3]
        aa.append(CODON_TABLE[codon])
    return "".join(aa)


EXPECTED_NAMES = {
    "FLAG", "HA", "Myc", "V5", "ALFA", "T7", "6xHis", "8xHis",
    "StrepII", "TwinStrep",
    "GGGGS", "GS", "EAAAK", "XTEN",
    "P2A", "T2A", "E2A", "F2A",
    "PuroR", "ZeocinR", "ZEOCINr",
    "NLS_SV40", "bNLS", "NES",
    "Kozak",
    "STOP_TAA", "STOP_TAG", "STOP_TGA",
    "EGFP", "mEGFP",
}


def _load():
    assert LIBRARY.is_file(), f"missing {LIBRARY}"
    return yaml.safe_load(LIBRARY.read_text())


def test_library_file_exists_and_is_mapping_of_motifs():
    data = _load()
    assert isinstance(data, dict)
    assert "motifs" in data
    assert isinstance(data["motifs"], list)
    assert len(data["motifs"]) >= 30, f"need >=30 motifs, got {len(data['motifs'])}"


def test_all_expected_motifs_present():
    data = _load()
    names = {m["name"] for m in data["motifs"]}
    missing = EXPECTED_NAMES - names
    assert not missing, f"missing motifs: {missing}"


def test_no_hyphens_in_motif_names():
    data = _load()
    for motif in data["motifs"]:
        assert "-" not in motif["name"], (
            f"motif name '{motif['name']}' contains '-', which is reserved as grammar separator"
        )


def test_every_motif_has_required_fields():
    data = _load()
    for motif in data["motifs"]:
        for field in ("name", "category", "dna", "aa", "description"):
            assert field in motif, f"motif {motif.get('name')} missing '{field}'"


def test_dna_contains_only_acgt():
    data = _load()
    for motif in data["motifs"]:
        assert set(motif["dna"].upper()) <= set("ACGT"), (
            f"motif {motif['name']} has non-ACGT chars in DNA"
        )


@pytest.mark.parametrize("index", range(60))
def test_dna_translates_to_listed_aa(index):
    data = _load()
    if index >= len(data["motifs"]):
        pytest.skip("no motif at this index")
    motif = data["motifs"][index]
    dna = motif["dna"].upper()
    expected_aa = motif["aa"]

    # Skip motifs marked as 5'UTR / regulatory (no translation check)
    if motif.get("category") in {"regulatory", "stop"}:
        # Stops still have 3-char DNA that should translate to '*'.
        if motif.get("category") == "stop":
            assert len(dna) == 3
            assert CODON_TABLE[dna] == "*"
            assert expected_aa == "*"
        return

    # DNA length must be a multiple of 3 for coding motifs.
    assert len(dna) % 3 == 0, f"{motif['name']}: DNA length {len(dna)} not multiple of 3"

    translated = _translate(dna)
    # Strip trailing stop if it's included in the DNA.
    if translated.endswith("*"):
        translated = translated[:-1]
    assert translated == expected_aa, (
        f"{motif['name']}: DNA translates to {translated!r}, expected {expected_aa!r}"
    )
