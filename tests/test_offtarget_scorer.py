"""Tests for :mod:`pegrna_pipeline.hpc.offtarget_scorer` (Plan 3 T9)."""

from __future__ import annotations

from typing import List

from pegrna_pipeline.hpc.offtarget_scorer import (
    parse_offtarget_tsv,
    submit_offtarget_job,
)


CANNED_TSV = (
    "guide_id\tspacer\tdeepspcas9\tdeephf\tcfd_specificity\t"
    "mm0_count\tmm1_count\tmm2_count\tmm3_count\n"
    "g1\tACGTACGTACGTACGTACGT\t0.85\t0.91\t0.97\t1\t0\t3\t12\n"
    "g2\tTTTTCCCCGGGGAAAACCCC\t0.22\t0.31\t0.50\t1\t2\t5\t20\n"
)


def test_parse_offtarget_tsv_basic():
    reports = parse_offtarget_tsv(CANNED_TSV)
    assert len(reports) == 2
    r = reports[0]
    assert r.guide_id == "g1"
    assert r.spacer == "ACGTACGTACGTACGTACGT"
    assert r.deepspcas9 == 0.85
    assert r.deephf == 0.91
    assert r.cfd_specificity == 0.97
    assert r.mm0_count == 1
    assert r.mm3_count == 12


def test_parse_offtarget_tsv_empty_body():
    tsv = (
        "guide_id\tspacer\tdeepspcas9\tdeephf\tcfd_specificity\t"
        "mm0_count\tmm1_count\tmm2_count\tmm3_count\n"
    )
    assert parse_offtarget_tsv(tsv) == []


def test_parse_offtarget_tsv_handles_na():
    tsv = (
        "guide_id\tspacer\tdeepspcas9\tdeephf\tcfd_specificity\t"
        "mm0_count\tmm1_count\tmm2_count\tmm3_count\n"
        "g1\tACGT\tNA\tNA\tNA\t1\t0\t0\t0\n"
    )
    reports = parse_offtarget_tsv(tsv)
    assert len(reports) == 1
    # NA values are materialised as 0.0 (NaN would also be acceptable;
    # our convention is "treat NA as missing -> 0.0" for now).
    assert reports[0].deepspcas9 == 0.0


class FakeSSH:
    def __init__(self, tsv_output, rscript_rc=0):
        self.commands: List[str] = []
        self.uploads: List[tuple] = []
        self.downloads: List[tuple] = []
        self._tsv = tsv_output
        self._rc = rscript_rc
        self._qsub_ids = iter([f"{2000 + i}.pbs" for i in range(100)])

    def upload(self, local, remote):
        self.uploads.append((local, remote))

    def download(self, remote, local):
        self.downloads.append((remote, local))
        with open(local, "w") as fh:
            fh.write(self._tsv)

    def exec(self, cmd: str):
        self.commands.append(cmd)
        if cmd.startswith("qsub "):
            return (next(self._qsub_ids) + "\n", "", 0)
        if cmd.startswith("qstat"):
            return ("", "", 0)
        if cmd.startswith("grep "):
            return ("1\n", "", 0)  # marker present
        return ("", "", 0)


def test_submit_offtarget_job_returns_reports(tmp_path, monkeypatch):
    ssh = FakeSSH(CANNED_TSV)
    # Force download to a temp path we control
    monkeypatch.chdir(tmp_path)
    guides = [{"guide_id": "g1", "spacer": "ACGT"}, {"guide_id": "g2", "spacer": "TTTT"}]
    reports = submit_offtarget_job(
        ssh=ssh,
        guides=guides,
        remote_run_dir="/work/run",
        rscript_path="/home/mj/pegRNA-pipeline/scripts/score_guides.R",
        conda_env="pegrna-ngrna",
        walltime_hours=1,
    )
    assert len(reports) == 2
    assert reports[0].guide_id == "g1"
    # qsub was called exactly once
    assert sum(1 for c in ssh.commands if c.startswith("qsub ")) == 1
    # guides CSV was uploaded
    csv_uploads = [u for u in ssh.uploads if u[1].endswith(".csv")]
    assert len(csv_uploads) == 1


def test_submit_offtarget_job_uses_configured_nuclease():
    ssh = FakeSSH(CANNED_TSV)
    submit_offtarget_job(
        ssh=ssh,
        guides=[{"guide_id": "g1", "spacer": "ACGT"}],
        remote_run_dir="/work/run",
        rscript_path="/home/mj/score_guides.R",
        conda_env="pegrna-ngrna",
        nuclease="SpCas9",
        walltime_hours=1,
    )
    # The rendered PBS script referenced SpCas9 (uploaded as .pbs via exec heredoc).
    # We check the heredoc catted by exec commands.
    joined = "\n".join(ssh.commands)
    assert "SpCas9" in joined
