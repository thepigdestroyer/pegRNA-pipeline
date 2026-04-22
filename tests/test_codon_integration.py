"""End-to-end codon enumeration integration tests.

- enumerate_with_metadata("GGGGSDYK") yields CodonVariants that all translate
  back to the input aa sequence.
- optimise_bulk on the full 124-aa ZeocinR protein round-trips correctly.
"""

import pytest
from python_codon_tables import get_codons_table

from pegrna_pipeline.agent.motif_library import MotifLibrary
from pegrna_pipeline.hpc.codon_expander import (
    enumerate_with_metadata,
    optimise_bulk,
)


@pytest.fixture(scope="module")
def translate():
    table = {}
    for aa, codons in get_codons_table("h_sapiens_9606").items():
        for c in codons:
            table[c] = aa
    return table


def _translate_to_aa(dna: str, table: dict) -> str:
    return "".join(table[dna[i : i + 3]] for i in range(0, len(dna), 3))


def test_ggggsdyk_variants_all_translate_back(translate):
    aa = "GGGGSDYK"
    result = enumerate_with_metadata(aa, max_variants=100)
    assert len(result.variants) > 0, "must have at least one clean variant"
    for v in result.variants:
        assert _translate_to_aa(v.dna, translate) == aa, (
            f"variant {v.dna!r} did not translate back to {aa!r}"
        )


def test_ggggsdyk_theoretical_max_matches_manual_computation():
    aa = "GGGGSDYK"
    # G has 4 codons, S has 6, D has 2, Y has 2, K has 2.  GGGGS = 4*4*4*4*6,
    # then DYK = 2*2*2.
    expected = (4 ** 4) * 6 * 2 * 2 * 2
    result = enumerate_with_metadata(aa, max_variants=5)
    assert result.theoretical_max == expected


def test_zeocinr_optimise_bulk_round_trips(translate):
    lib = MotifLibrary.load_bundled()
    zeo = lib.lookup("ZeocinR")
    aa = zeo.aa
    assert len(aa) == 124
    dna = optimise_bulk(aa)
    assert len(dna) == 3 * len(aa)
    assert _translate_to_aa(dna, translate) == aa


def test_zeocinr_optimise_bulk_produces_clean_acgt(translate):
    lib = MotifLibrary.load_bundled()
    zeo = lib.lookup("ZeocinR")
    dna = optimise_bulk(zeo.aa)
    assert set(dna) <= set("ACGT")
