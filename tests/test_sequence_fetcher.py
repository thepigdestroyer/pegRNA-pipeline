"""Tests for :mod:`pegrna_pipeline.hpc.sequence_fetcher`.

The module wraps Bio.Entrez to discover a RefSeq mRNA for a gene, download
its GenBank record, extract the CDS / stop codon, and pull flanking context
around the stop codon. All network calls go through module-private helpers
so the unit tests can monkey-patch them without touching NCBI.
"""

from __future__ import annotations

import json
import subprocess
import sys

import pytest

from pegrna_pipeline.hpc import sequence_fetcher as sf
from pegrna_pipeline.hpc.data_types import GeneRecord


# ------------------------------------------------------------------ fixtures


# A tiny synthetic mRNA: 30 nt 5'UTR, then CDS of 30 nt ending in TAA, then
# 30 nt 3'UTR. Total length 90. CDS is [30, 60); stop codon starts at 57.
SYNTH_MRNA = (
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # 30 nt 5' UTR
    "ATGAAACCCGGGTTTAAACCCGGGCCCTAA"  # 30 nt CDS: M K P G F K P G P *
    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"  # 30 nt 3' UTR
)
CDS_START = 30
CDS_END = 60
STOP_START = 57  # position of T in TAA


@pytest.fixture
def patched_helpers(monkeypatch):
    """Replace the three network helpers with deterministic stubs."""

    def fake_search(gene_name, *, organism, api_key=None):
        return "NM_TEST_001.1"

    def fake_fetch(accession, *, email, api_key=None):
        return {
            "accession": accession,
            "sequence": SYNTH_MRNA,
            "cds_features": [(CDS_START, CDS_END)],
        }

    def fake_extract_cds(record):
        start, end = record["cds_features"][0]
        stop_codon = record["sequence"][end - 3 : end]
        return start, end, end - 3, stop_codon

    monkeypatch.setattr(sf, "_search_refseq_accession", fake_search)
    monkeypatch.setattr(sf, "_fetch_genbank_record", fake_fetch)
    monkeypatch.setattr(sf, "_extract_cds_coordinates", fake_extract_cds)
    return monkeypatch


# -------------------------------------------------------------- unit tests


def test_fetch_gene_record_returns_gene_record(patched_helpers):
    record = sf.fetch_gene_record("TEST1", email="x@y.z", flanking_bp=10)
    assert isinstance(record, GeneRecord)
    assert record.gene_name == "TEST1"
    assert record.refseq_accession == "NM_TEST_001.1"
    assert record.cds_start == CDS_START
    assert record.cds_end == CDS_END
    assert record.stop_codon_start == STOP_START
    assert record.stop_codon == "TAA"
    assert record.flanking_bp == 10


def test_flanks_anchored_on_stop_codon(patched_helpers):
    """Upstream flank ends immediately before stop; downstream starts after stop."""
    record = sf.fetch_gene_record("TEST1", email="x@y.z", flanking_bp=10)
    assert len(record.upstream_flank) == 10
    assert len(record.downstream_flank) == 10
    assert record.upstream_flank == SYNTH_MRNA[STOP_START - 10 : STOP_START]
    assert record.downstream_flank == SYNTH_MRNA[STOP_START + 3 : STOP_START + 13]


def test_flank_default_is_100_bp(monkeypatch):
    """With default flanking, we need a bigger synthetic mRNA — patch accordingly."""
    big_up = "A" * 200
    big_down = "G" * 200
    seq = big_up + "ATG" + "AAA" * 10 + "TAA" + big_down
    cds_start = 200
    cds_end = 200 + 3 + 30 + 3
    stop_start = cds_end - 3

    def fake_search(*a, **k):
        return "NM_X"

    def fake_fetch(*a, **k):
        return {"accession": "NM_X", "sequence": seq, "cds_features": [(cds_start, cds_end)]}

    def fake_extract(record):
        return cds_start, cds_end, stop_start, seq[stop_start:cds_end]

    monkeypatch.setattr(sf, "_search_refseq_accession", fake_search)
    monkeypatch.setattr(sf, "_fetch_genbank_record", fake_fetch)
    monkeypatch.setattr(sf, "_extract_cds_coordinates", fake_extract)

    record = sf.fetch_gene_record("X", email="x@y.z")
    assert record.flanking_bp == 100
    assert len(record.upstream_flank) == 100
    assert len(record.downstream_flank) == 100


def test_insufficient_flanking_raises(monkeypatch):
    """If mRNA is too short for the requested flank we raise InsufficientFlankingError."""
    short = "A" * 5 + "ATGTAA" + "G" * 5
    cds_start, cds_end = 5, 11
    stop_start = 8

    monkeypatch.setattr(sf, "_search_refseq_accession", lambda *a, **k: "NM_X")
    monkeypatch.setattr(
        sf,
        "_fetch_genbank_record",
        lambda *a, **k: {"accession": "NM_X", "sequence": short, "cds_features": [(cds_start, cds_end)]},
    )
    monkeypatch.setattr(
        sf,
        "_extract_cds_coordinates",
        lambda r: (cds_start, cds_end, stop_start, "TAA"),
    )
    with pytest.raises(sf.InsufficientFlankingError):
        sf.fetch_gene_record("X", email="x@y.z", flanking_bp=100)


def test_gene_not_found_propagates(monkeypatch):
    def fail_search(*a, **k):
        raise sf.RefSeqNotFoundError("nothing returned")

    monkeypatch.setattr(sf, "_search_refseq_accession", fail_search)
    with pytest.raises(sf.RefSeqNotFoundError):
        sf.fetch_gene_record("DOESNOTEXIST", email="x@y.z")


# ----------------------------------------------------- exception hierarchy


def test_exception_classes_inherit_from_sequence_fetch_error():
    assert issubclass(sf.GeneNotFoundError, sf.SequenceFetchError)
    assert issubclass(sf.RefSeqNotFoundError, sf.SequenceFetchError)
    assert issubclass(sf.InsufficientFlankingError, sf.SequenceFetchError)


# --------------------------------------------------------------- __main__

def test_main_block_emits_json(tmp_path, monkeypatch):
    """Running the module as a script with --json should print a GeneRecord as JSON."""
    # Monkey-patching the helpers is easier via env-var stubbing that the
    # subprocess respects, so we invoke a tiny wrapper that imports and
    # patches before calling the module main.
    wrapper = tmp_path / "run_main.py"
    wrapper.write_text(
        "from pegrna_pipeline.hpc import sequence_fetcher as sf\n"
        "from unittest.mock import patch\n"
        "mrna = 'A'*40 + 'ATGAAATAA' + 'G'*40\n"
        "cds_start, cds_end = 40, 49\n"
        "stop_start = 46\n"
        "with patch.object(sf, '_search_refseq_accession', lambda *a, **k: 'NM_TEST'),\\\n"
        "     patch.object(sf, '_fetch_genbank_record',\n"
        "                   lambda *a, **k: {'accession': 'NM_TEST', 'sequence': mrna,\n"
        "                                    'cds_features': [(cds_start, cds_end)]}),\\\n"
        "     patch.object(sf, '_extract_cds_coordinates',\n"
        "                   lambda r: (cds_start, cds_end, stop_start, 'TAA')):\n"
        "    sf._main(['TESTGENE', '--email', 'x@y.z', '--flanking-bp', '20', '--json'])\n"
    )
    result = subprocess.run(
        [sys.executable, str(wrapper)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["gene_name"] == "TESTGENE"
    assert payload["refseq_accession"] == "NM_TEST"
    assert payload["flanking_bp"] == 20
