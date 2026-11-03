"""Protospacers are not all the same length, and cas-offinder cannot cope.

cas-offinder takes ONE pattern per run and aborts the entire job if any query
differs from it in length:

    Critical error! The length of target sequences should match with the
    length of pattern sequence.

PRIDICT returns some 19-nt ngRNA spacers alongside the usual 20-nt ones, so a
single mixed batch lost every guide's result, not just the odd one. This only
surfaced once the off-target gate ran a real cas-offinder: with the clean-pass
stand-in, mixed lengths were silently fine.
"""

from pegrna_pipeline.hpc.offtarget_scorer import CasOffinderScorer


class _RecordingScorer(CasOffinderScorer):
    """Captures each cas-offinder invocation instead of running the binary."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.invocations: list[dict] = []

    def _run(self, input_path: str, output_path: str) -> None:
        with open(input_path) as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
        genome, pattern, *queries = lines
        self.invocations.append(
            {
                "pattern": pattern,
                "queries": [q.split()[0] for q in queries],
            }
        )
        # Every guide matches itself exactly once, and nothing else.
        with open(output_path, "w") as fh:
            for q in queries:
                crrna = q.split()[0]
                fh.write(f"{crrna}\tchr1\t1000\t{crrna}\t+\t0\n")


def _scorer():
    return _RecordingScorer(genome="/genomes/hg38.2bit", pam="NGG", max_mismatches=3)


def test_mixed_lengths_are_split_into_separate_runs():
    s = _scorer()
    guides = ["A" * 20, "C" * 19, "G" * 20, "T" * 19]
    s(guides)

    assert len(s.invocations) == 2, "one cas-offinder run per distinct length"
    patterns = sorted(i["pattern"] for i in s.invocations)
    assert patterns == ["N" * 19 + "NGG", "N" * 20 + "NGG"]


def test_each_query_matches_its_own_pattern_length():
    s = _scorer()
    s(["A" * 20, "C" * 19])
    for inv in s.invocations:
        want = len(inv["pattern"])
        for q in inv["queries"]:
            assert len(q) == want, f"{q} does not match pattern length {want}"


def test_every_guide_still_gets_a_result():
    s = _scorer()
    guides = ["A" * 20, "C" * 19, "G" * 20]
    out = s(guides)
    assert set(out) == set(guides)
    # Exactly one exact hit each => on-target found, no off-targets.
    for g in guides:
        on, offtargets, has_off = out[g]
        assert on == 1.0
        assert offtargets == 0
        assert has_off is False


def test_uniform_lengths_still_use_a_single_run():
    """The common case must not become three invocations of one guide each."""
    s = _scorer()
    s(["A" * 20, "C" * 20, "G" * 20])
    assert len(s.invocations) == 1
    assert len(s.invocations[0]["queries"]) == 3
