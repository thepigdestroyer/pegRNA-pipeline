"""cas-offinder is the single most expensive thing a design does.

One whole-genome scan took ~100 s of wall time but 69 CPU-MINUTES, spread over ~40 cores of a
shared 64-core machine. The public demo scores the same handful of guides every time, so:

* results are cached per guide (keyed on the genome file, PAM and mismatch settings), and a
  repeat run does not invoke the binary at all;
* when it does run, it is capped to a thread budget instead of taking every core.
"""

import json

from pegrna_pipeline.hpc.offtarget_scorer import CasOffinderScorer


class _Fake(CasOffinderScorer):
    """Stands in for the binary: every guide matches itself once, plus a scripted number of
    1- and 2-mismatch near-matches. Records each invocation."""

    def __init__(self, *, near=None, **kw):
        super().__init__(**kw)
        self.runs = []
        self.near = near or {}

    def _run(self, input_path, output_path):
        with open(input_path) as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
        _genome, _pattern, *queries = lines
        guides = [q.split()[0] for q in queries]
        self.runs.append(guides)
        with open(output_path, "w") as fh:
            for crrna in guides:
                g = crrna[:-3]
                fh.write(f"{crrna}\tchr1\t1\t{crrna}\t+\t0\n")
                for i in range(self.near.get(g, 0)):
                    fh.write(f"{crrna}\tchr2\t{i}\t{crrna}\t+\t2\n")


G1, G2, G3 = "A" * 20, "C" * 20, "G" * 20


def _scorer(tmp_path, genome="/genomes/hg38.2bit", **kw):
    return _Fake(genome=genome, pam="NGG", max_mismatches=3, cache_dir=str(tmp_path), **kw)


def test_a_repeat_run_does_not_invoke_the_binary(tmp_path):
    s = _scorer(tmp_path)
    first = s([G1, G2])
    second = s([G1, G2])
    assert len(s.runs) == 1, "the second call must be served from the cache"
    assert first == second


def test_only_uncached_guides_are_scanned(tmp_path):
    s = _scorer(tmp_path)
    s([G1])
    s([G1, G2, G3])
    assert s.runs[1] == [G2 + "NNN", G3 + "NNN"], "G1 was cached; only G2 and G3 should be scanned"


def test_cached_results_keep_their_off_target_counts(tmp_path):
    s = _scorer(tmp_path, near={G2: 5})
    fresh = s([G1, G2])
    cached = s([G1, G2])
    assert fresh[G2][1] == 5 and cached[G2][1] == 5
    assert fresh == cached


def test_cache_persists_across_scorer_instances(tmp_path):
    _scorer(tmp_path)([G1])
    again = _scorer(tmp_path)
    again([G1])
    assert again.runs == []


def test_a_different_genome_does_not_share_results(tmp_path):
    a = _scorer(tmp_path, genome="/genomes/hg38.2bit")
    b = _scorer(tmp_path, genome="/genomes/other.2bit")
    a([G1])
    b([G1])
    assert len(b.runs) == 1, "results for one genome must never answer for another"


def test_different_mismatch_settings_do_not_share_results(tmp_path):
    a = _Fake(genome="g", pam="NGG", max_mismatches=3, cache_dir=str(tmp_path))
    b = _Fake(genome="g", pam="NGG", max_mismatches=4, cache_dir=str(tmp_path))
    a([G1])
    b([G1])
    assert len(b.runs) == 1


def test_a_modified_genome_file_invalidates_the_cache(tmp_path):
    g = tmp_path / "genome.2bit"
    g.write_bytes(b"x" * 100)
    s1 = _scorer(tmp_path / "cache", genome=str(g))
    s1([G1])
    g.write_bytes(b"y" * 250)  # size and mtime change
    s2 = _scorer(tmp_path / "cache", genome=str(g))
    s2([G1])
    assert len(s2.runs) == 1


def test_a_corrupt_cache_file_is_ignored_not_fatal(tmp_path):
    s = _scorer(tmp_path)
    s([G1])
    for f in tmp_path.rglob("*.json"):
        f.write_text("{ this is not json")
    again = _scorer(tmp_path)
    assert again([G1])[G1][0] == 1.0
    assert len(again.runs) == 1


def test_cache_is_off_unless_configured(tmp_path, monkeypatch):
    """No env var and no cache_dir -> nothing is written (keeps unit tests hermetic)."""
    monkeypatch.delenv("PEGRNA_CACHE_DIR", raising=False)
    s = _Fake(genome="g", pam="NGG", max_mismatches=3)
    s([G1])
    s([G1])
    assert len(s.runs) == 2


def test_cache_file_is_plain_inspectable_json(tmp_path):
    _scorer(tmp_path)([G1])
    files = list(tmp_path.rglob("*.json"))
    assert files
    data = json.loads(files[0].read_text())
    assert isinstance(data, dict) and data


# ---- thread budget for the real binary -------------------------------------------------


def test_the_binary_runs_with_a_capped_thread_pool(monkeypatch, tmp_path):
    """Without a cap, PoCL (cas-offinder's CPU backend) takes every core."""
    seen = {}

    def fake_run(cmd, **kw):
        seen["env"] = kw.get("env")

        class R:
            returncode = 0
            stderr = stdout = ""

        open(cmd[-1], "w").close()
        return R()

    import shutil
    import subprocess

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(shutil, "which", lambda b: "/usr/bin/cas-offinder")
    monkeypatch.delenv("POCL_MAX_PTHREAD_COUNT", raising=False)
    monkeypatch.delenv("PEGRNA_CASOFFINDER_THREADS", raising=False)
    CasOffinderScorer(genome="g", pam="NGG")._run(str(tmp_path / "in"), str(tmp_path / "out"))
    assert int(seen["env"]["POCL_MAX_PTHREAD_COUNT"]) == 16


def test_the_thread_budget_is_configurable(monkeypatch, tmp_path):
    seen = {}

    def fake_run(cmd, **kw):
        seen["env"] = kw.get("env")

        class R:
            returncode = 0
            stderr = stdout = ""

        open(cmd[-1], "w").close()
        return R()

    import shutil
    import subprocess

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(shutil, "which", lambda b: "/usr/bin/cas-offinder")
    monkeypatch.setenv("PEGRNA_CASOFFINDER_THREADS", "4")
    CasOffinderScorer(genome="g", pam="NGG")._run(str(tmp_path / "in"), str(tmp_path / "out"))
    assert seen["env"]["POCL_MAX_PTHREAD_COUNT"] == "4"


def test_an_explicit_pocl_setting_is_respected(monkeypatch, tmp_path):
    seen = {}

    def fake_run(cmd, **kw):
        seen["env"] = kw.get("env")

        class R:
            returncode = 0
            stderr = stdout = ""

        open(cmd[-1], "w").close()
        return R()

    import shutil
    import subprocess

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(shutil, "which", lambda b: "/usr/bin/cas-offinder")
    monkeypatch.setenv("POCL_MAX_PTHREAD_COUNT", "2")
    CasOffinderScorer(genome="g", pam="NGG")._run(str(tmp_path / "in"), str(tmp_path / "out"))
    assert seen["env"]["POCL_MAX_PTHREAD_COUNT"] == "2"
