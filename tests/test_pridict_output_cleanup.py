"""PRIDICT's per-run output directory must not accumulate on a long-lived host.

Each design leaves predictions/<uuid>/ (~1 MB) that nothing reads once the CSVs are
parsed. On the shared deployment box -- whose system disk is already at 86% -- that
is an unbounded leak, so run_batch removes it. PEGRNA_KEEP_PRIDICT_OUTPUT=1 keeps it.
"""

from pegrna_pipeline.hpc.pridict_local import PridictLocalRunner


class _Runner(PridictLocalRunner):
    """Runs everything except the PRIDICT subprocess; leaves a fake output behind."""

    def _run_subprocess(self, input_basename, reldir):
        out = self.root / reldir
        (out / "marker.txt").write_text("pridict output")

    def _collect(self, outdir, rows):
        assert (outdir / "marker.txt").exists(), "output must still exist while being read"
        return {"s1": {"pegrnas": [], "nicking": []}}


def _runner(tmp_path):
    r = _Runner(pridict_dir=str(tmp_path), pridict_python="/usr/bin/python3")
    r.root = tmp_path
    return r


def test_output_directory_is_removed_after_collection(tmp_path, monkeypatch):
    monkeypatch.delenv("PEGRNA_KEEP_PRIDICT_OUTPUT", raising=False)
    _runner(tmp_path).run_batch([{"sequence_name": "s1", "editseq": "ACGT"}])
    assert list((tmp_path / "predictions").iterdir()) == []


def test_input_csv_is_still_removed(tmp_path, monkeypatch):
    monkeypatch.delenv("PEGRNA_KEEP_PRIDICT_OUTPUT", raising=False)
    _runner(tmp_path).run_batch([{"sequence_name": "s1", "editseq": "ACGT"}])
    assert list((tmp_path / "input").glob("*.csv")) == []


def test_keep_flag_preserves_the_output_for_debugging(tmp_path, monkeypatch):
    monkeypatch.setenv("PEGRNA_KEEP_PRIDICT_OUTPUT", "1")
    _runner(tmp_path).run_batch([{"sequence_name": "s1", "editseq": "ACGT"}])
    kept = list((tmp_path / "predictions").iterdir())
    assert len(kept) == 1 and (kept[0] / "marker.txt").exists()
