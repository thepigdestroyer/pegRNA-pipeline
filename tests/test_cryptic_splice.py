"""MaxEntScan cryptic-splice-site gate (ADR 0011).

The gate scans the *variant region* of the edited construct (inserted cassette +
recoded bases + spanning junctions) for a donor/acceptor motif at/above the strong
MaxEnt threshold that the edit created or strengthened vs the reference. The real
MaxEntScan backend is CX3-only (maxentpy); these tests drive the gate logic with a
deterministic fake that maps specific motif windows to scores.
"""

from pegrna_pipeline.hpc.cryptic_splice import (
    MaxEntScanFake,
    passes_cryptic_splice_gate,
)

PREFIX = "CCCCCC"
SUFFIX = "GGGGGG"


def _pair(insert):
    """A reference (prefix+suffix) and an edited allele with `insert` spliced in."""
    return PREFIX + SUFFIX, PREFIX + insert + SUFFIX


def test_no_strong_site_passes():
    ref, edited = _pair("AAAGTAAAAA")
    ok, why = passes_cryptic_splice_gate(ref, edited, scorer=MaxEntScanFake())
    assert ok
    assert "no cryptic" in why


def test_strong_donor_in_variant_region_fails():
    ref, edited = _pair("AAAGTAAAAA")  # donor window "AAAGTAAAA" lands in the insert
    fake = MaxEntScanFake(score5_map={"AAAGTAAAA": 10.0})
    ok, why = passes_cryptic_splice_gate(ref, edited, scorer=fake)
    assert not ok
    assert "cryptic donor" in why
    assert "10.0" in why


def test_strong_acceptor_in_variant_region_fails():
    ref = "C" * 20 + "TTTTTT"
    edited = "C" * 20 + "AGAAA" + "TTTTTT"  # AG at index 20; window C*18 + AGAAA
    window = "C" * 18 + "AGAAA"
    fake = MaxEntScanFake(score3_map={window: 12.0})
    ok, why = passes_cryptic_splice_gate(ref, edited, scorer=fake)
    assert not ok
    assert "cryptic acceptor" in why


def test_below_threshold_site_passes():
    ref, edited = _pair("AAAGTAAAAA")
    fake = MaxEntScanFake(score5_map={"AAAGTAAAA": 7.9})  # just under strong=8
    ok, _why = passes_cryptic_splice_gate(ref, edited, scorer=fake)
    assert ok


def test_native_site_not_stronger_than_reference_passes():
    # The same strong donor motif exists in the shared suffix (so it is present in
    # BOTH ref and edited); the edit does not strengthen it beyond reference -> pass.
    # Insert carries a donor window scored equal to the native one -> not "created".
    ref, edited = _pair("AAAGTAAAAA")
    # Give the insert donor score 8.0 AND make the reference already contain an
    # equally strong donor by mapping a ref-present window to 8.0.
    # ref = "CCCCCCGGGGGG" has no GT, so add one to both via a shared motif:
    ref2 = "CCCGTCCCGGGGGG"  # contains GT at idx 3 -> donor window "CCCGTCCC"...
    edited2 = "CCCGTCCC" + "AAAGTAAAAA" + "GGGGGG"
    shared_window = "CCCGTCCCG"[:9]  # ref donor window around its GT
    fake = MaxEntScanFake(score5_map={"AAAGTAAAA": 8.0, shared_window: 8.0}, default=0.0)
    ok, _why = passes_cryptic_splice_gate(ref2, edited2, scorer=fake)
    # created requires score strictly greater than the best reference site (8.0),
    # so an equal-strength insert site is not flagged.
    assert ok


def test_strong_site_in_shared_flank_not_flagged():
    # A strong donor entirely within the shared prefix is native (present in ref),
    # so ref_best captures it and it is not attributed to the edit.
    ref = "AAGTACGTACAAAAAA"
    edited = "AAGTACGTAC" + "TTTTT" + "AAAAAA"
    # donor window around the prefix GT (index 2): "AAGTACGTA"
    fake = MaxEntScanFake(score5_map={"AAGTACGTA": 11.0})
    ok, _why = passes_cryptic_splice_gate(ref, edited, scorer=fake)
    assert ok
