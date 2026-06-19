"""Task 12: PBS2 joint length/GC/Tm + repetitive gate."""

from pegrna_pipeline.hpc.recoded_tjpe_workflow import passes_pbs2_gate


def test_pbs2_in_band_passes():
    # 13 nt, GC 0.54, Tm ~40.8 C (computed, in 37-42 band)
    ok, why = passes_pbs2_gate("ACGTGACGTGACA", repetitive=False)
    assert ok is True, why


def test_pbs2_rejects_too_short():
    ok, why = passes_pbs2_gate("ACGTGCAC", repetitive=False)  # 8 nt < 11
    assert ok is False and "length" in why


def test_pbs2_rejects_too_long():
    ok, why = passes_pbs2_gate("ACGTGACGTGACAACGTGA", repetitive=False)  # 19 nt > 17
    assert ok is False and "length" in why


def test_pbs2_rejects_repetitive():
    ok, why = passes_pbs2_gate("ACGTGACGTGACA", repetitive=True)
    assert ok is False and "repetitive" in why


def test_pbs2_rejects_out_of_band_tm():
    # very GC-rich 15-mer => Tm well above 42 C
    ok, why = passes_pbs2_gate("GCGCGCGCGCGCGCG", repetitive=False)
    assert ok is False and ("Tm" in why or "GC" in why)


def test_pbs2_rejects_low_gc():
    ok, why = passes_pbs2_gate("ATATATATATATA", repetitive=False)  # GC ~0
    assert ok is False and ("GC" in why or "Tm" in why)
