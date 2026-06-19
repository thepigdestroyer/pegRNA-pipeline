from pegrna_pipeline.hpc.microhomology import longest_common_run, max_identity_run


def test_longest_common_run_finds_shared_substring():
    # 12-nt identical stretch embedded in otherwise different sequences
    a = "AAAACGTACGTACGTAGGGG"
    b = "TTTTCGTACGTACGTATTTT"
    assert longest_common_run(a, b) == len("CGTACGTACGTA")  # 12


def test_max_identity_run_against_self_offset():
    # identical sequences => full-length run
    s = "ACGTACGTACGT"
    assert max_identity_run(s, s) == len(s)


def test_no_shared_run_returns_zero_or_one():
    assert longest_common_run("AAAA", "GGGG") == 0


def test_empty_inputs():
    assert longest_common_run("", "ACGT") == 0
    assert longest_common_run("ACGT", "") == 0


def test_case_insensitive():
    assert longest_common_run("acgtACGT", "ACGTACGT") == 8


def test_single_base_match():
    assert longest_common_run("AAAA", "GGAGG") == 1
