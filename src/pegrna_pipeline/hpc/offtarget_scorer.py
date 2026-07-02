"""Off-target scoring submitter.

Submits an ``Rscript scripts/score_guides.R`` PBS job, waits for it,
downloads the output TSV, and parses it into a list of
:class:`OffTargetReport`.
"""

from __future__ import annotations

import csv
import io
import os
import tempfile
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

OFFTARGET_SCORING_DONE = "OFFTARGET_SCORING_DONE"


@dataclass(frozen=True)
class OffTargetReport:
    """One row of :file:`offtarget_scores.tsv`.

    On-target score columns (deepspcas9, deephf, cfd_specificity) and
    mismatch counts are ``Optional[float]`` / ``Optional[int]`` so a
    missing-score / unscanned guide can be distinguished from a real
    zero. Renderers print ``"NA"`` for ``None``.
    """

    guide_id: str
    spacer: str
    deepspcas9: Optional[float]
    deephf: Optional[float]
    cfd_specificity: Optional[float]
    mm0_count: Optional[int]
    mm1_count: Optional[int]
    mm2_count: Optional[int]
    mm3_count: Optional[int]


def _to_float(value: str) -> Optional[float]:
    """Parse a float; ``None`` for empty / NA / NaN / NULL / unparseable."""
    if value is None:
        return None
    v = value.strip()
    if v == "" or v.upper() in {"NA", "NAN", "NULL"}:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _to_int(value: str) -> Optional[int]:
    """Parse an int; ``None`` for empty / NA / NaN / NULL / unparseable."""
    v = (value or "").strip()
    if v == "" or v.upper() in {"NA", "NAN", "NULL"}:
        return None
    try:
        return int(v)
    except ValueError:
        return None


def _fmt_float(value: Optional[float], width: int = 3) -> str:
    """``"%.3f" % value`` or ``"NA"`` for None."""
    if value is None:
        return "NA"
    return f"{value:.{width}f}"


def _fmt_int(value: Optional[int]) -> str:
    """``str(value)`` or ``"NA"`` for None."""
    return "NA" if value is None else str(value)


class CasOffinderScorer:
    """Real genome-wide off-target specificity via the ``cas-offinder`` binary
    (Bae et al. 2014, Bioinformatics 30:1473).

    Callable matching :func:`nick_search.screen_nick_candidates`'s contract:
    ``protospacers -> {protospacer: (ontarget_proxy, offtarget_count, has_offtarget)}``.
    For each guide it enumerates all genomic sites within ``max_mismatches`` of
    ``spacer + PAM`` and aggregates:

      * ``offtarget_count`` = disqualifying genomic hits beyond the single intended
        on-target: extra EXACT (0-mismatch) duplicates PLUS near-tier hits with
        ``0 < mismatches <= count_mismatch_threshold``. Distant hits (mismatches >
        ``count_mismatch_threshold``, e.g. 3-mismatch sites) are still SEARCHED up to
        ``max_mismatches`` but are NOT counted — they are not functional off-targets
        and, against a whole 3 Gb genome, every 20-mer has dozens of them, which would
        otherwise reject every guide,
      * ``has_offtarget`` = ``offtarget_count > 0``,
      * ``ontarget_proxy`` = 1.0 iff the guide has >=1 exact (0-mismatch) genomic
        hit (it exists uniquely at its locus), else 0.0.

    ``count_mismatch_threshold`` (default 2) is the gate stringency: the mismatch
    depth at or below which a hit counts as an off-target. It is decoupled from
    ``max_mismatches`` (the search depth) so the search can characterise the full
    off-target profile while the gate reflects real cleavage risk.

    On-target ACTIVITY (DeepSpCas9/DeepHF) is deliberately NOT this adapter's job
    — Cas-OFFinder scores SPECIFICITY only; layer ``score_guides.R`` for activity.
    Requires the ``cas-offinder`` binary + a genome FASTA (or 2bit dir) on CX3;
    raises loud if the binary is absent (never silently passes a guide).
    """

    def __init__(
        self,
        *,
        genome: str,
        pam: str = "NGG",
        spacer_len: int = 20,
        max_mismatches: int = 3,
        count_mismatch_threshold: int = 2,
        device: str = "C",  # cas-offinder: C=CPU, G=GPU
        binary: str = "cas-offinder",
    ) -> None:
        self.genome = genome
        self.pam = pam.upper()
        self.spacer_len = spacer_len
        self.max_mismatches = max_mismatches
        self.count_mismatch_threshold = count_mismatch_threshold
        self.device = device
        self.binary = binary

    def _run(self, input_path: str, output_path: str) -> None:  # pragma: no cover - CX3 only
        """Invoke cas-offinder. Overridden in tests to inject canned output."""
        import shutil  # noqa: PLC0415
        import subprocess  # noqa: PLC0415

        if shutil.which(self.binary) is None:
            raise ImportError(
                f"CasOffinderScorer requires the {self.binary!r} binary + a genome; "
                "install cas-offinder and stage the genome on CX3 (never silently skipped)."
            )
        proc = subprocess.run(
            [self.binary, input_path, self.device, output_path],
            check=False,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or "").strip()[-800:]
            raise RuntimeError(f"cas-offinder failed (rc={proc.returncode}):\n{tail}")

    def _write_input(self, path: str, protospacers: Sequence[str]) -> None:
        pam_ns = "N" * len(self.pam)
        pattern = "N" * self.spacer_len + self.pam
        lines = [self.genome, pattern]
        for sp in protospacers:
            lines.append(f"{sp}{pam_ns} {self.max_mismatches}")
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")

    def __call__(self, protospacers: List[str]) -> "dict[str, Tuple[float, int, bool]]":
        results: dict = {sp: {"exact": 0, "near": 0} for sp in protospacers}
        if not protospacers:
            return {}
        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "casin.txt")
            outp = os.path.join(tmp, "casout.txt")
            self._write_input(inp, protospacers)
            self._run(inp, outp)
            for row in self._parse_output(outp):
                guide = row["guide"]
                if guide not in results:
                    continue
                mm = row["mismatches"]
                if mm == 0:
                    results[guide]["exact"] += 1
                elif mm <= self.count_mismatch_threshold:
                    results[guide]["near"] += 1
                # mm > count_mismatch_threshold: searched (up to max_mismatches) but
                # NOT counted -- distant hits are not functional off-targets.
        out: dict = {}
        for sp, agg in results.items():
            exact, near = agg["exact"], agg["near"]
            # Disqualifying off-targets = extra EXACT duplicates + near-tier hits.
            offtargets = max(0, exact - 1) + near
            on = 1.0 if exact >= 1 else 0.0
            out[sp] = (on, offtargets, offtargets > 0)
        return out

    def _parse_output(self, path: str):
        """Yield {'guide', 'mismatches'} rows from a cas-offinder output file.

        cas-offinder columns are whitespace/tab separated:
        ``crRNA  Chrom  Position  MatchedSeq  Strand  Mismatches``. The crRNA is
        the query (guide + PAM Ns); we strip the trailing PAM Ns to recover the
        guide. Header lines / blanks are skipped.
        """
        pam_len = len(self.pam)
        with open(path) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 6:
                    continue
                crrna = parts[0].upper()
                try:
                    mism = int(parts[-1])
                except ValueError:
                    continue  # header row
                guide = crrna[: len(crrna) - pam_len] if pam_len else crrna
                yield {"guide": guide, "mismatches": mism}


def parse_offtarget_tsv(tsv_text: str) -> List[OffTargetReport]:
    """Parse the TSV emitted by ``score_guides.R`` into reports."""
    reader = csv.DictReader(io.StringIO(tsv_text), delimiter="\t")
    out: List[OffTargetReport] = []
    for row in reader:
        out.append(
            OffTargetReport(
                guide_id=row["guide_id"],
                spacer=row["spacer"],
                deepspcas9=_to_float(row.get("deepspcas9", "")),
                deephf=_to_float(row.get("deephf", "")),
                cfd_specificity=_to_float(row.get("cfd_specificity", "")),
                mm0_count=_to_int(row.get("mm0_count", "")),
                mm1_count=_to_int(row.get("mm1_count", "")),
                mm2_count=_to_int(row.get("mm2_count", "")),
                mm3_count=_to_int(row.get("mm3_count", "")),
            )
        )
    return out


def _build_offtarget_pbs_body(
    *,
    conda_env: str,
    rscript_path: str,
    guides_csv: str,
    output_tsv: str,
    done_flag: str,
    nuclease: str,
    scores: Sequence[str],
    genome: str,
    max_mismatches: int,
) -> str:
    score_arg = ",".join(scores)
    return (
        "set -euo pipefail\n"
        "module load miniforge/3\n"
        # Imperial CX3 ships only a miniforge bootstrap; conda lives in
        # ~/miniforge3 (created by setup_hpc.sh).
        'source "$HOME/miniforge3/etc/profile.d/conda.sh"\n'
        f"conda activate {conda_env}\n"
        f"Rscript {rscript_path} "
        f"--guides {guides_csv} "
        f"--nuclease {nuclease} "
        f"--scores {score_arg} "
        f"--genome {genome} "
        f"--max-mismatches {max_mismatches} "
        f"--output {output_tsv}\n"
        # Write sentinel to shared RDS filesystem (not PBS log, which is staged
        # asynchronously and may not exist when the driver checks immediately
        # after poll_until_done returns).
        f"touch {done_flag}\n"
    )


def _write_guides_csv_and_upload(ssh, guides, remote_path: str) -> None:
    """Write an in-memory guides CSV and upload it."""
    with tempfile.NamedTemporaryFile(mode="w", newline="", suffix=".csv", delete=False) as tmp:
        writer = csv.writer(tmp)
        writer.writerow(["guide_id", "spacer"])
        for g in guides:
            writer.writerow([g["guide_id"], g["spacer"]])
        local = tmp.name
    try:
        ssh.upload(local, remote_path)
    finally:
        os.unlink(local)


def submit_offtarget_job(
    *,
    ssh,
    guides: Sequence[dict],
    remote_run_dir: str,
    rscript_path: str,
    conda_env: str = "pegrna-ngrna",
    nuclease: str = "SpCas9",
    scores: Tuple[str, ...] = ("DeepSpCas9", "DeepHF", "CFD"),
    genome: str = "GRCh38",
    max_mismatches: int = 3,
    walltime_hours: int = 6,
) -> List[OffTargetReport]:
    """Submit the off-target scorer, wait for it, and return parsed reports."""
    from pegrna_pipeline.hpc.pbs import (
        PBSJob,
        PBSJobSpec,
        poll_until_done,
        qsub,
        render_pbs_script,
    )

    guides_csv = f"{remote_run_dir}/guides.csv"
    output_tsv = f"{remote_run_dir}/offtarget_scores.tsv"
    # Sentinel written to shared RDS filesystem by the PBS job body.
    # Checking this avoids the PBS log staging race (PBS Pro stages stdout
    # asynchronously; the file may not exist immediately after qstat shows done).
    done_flag = f"{remote_run_dir}/OFFTARGET_DONE"
    log_path = f"{remote_run_dir}/offtarget_pbs.log"
    script_path = f"{remote_run_dir}/offtarget.pbs"

    _write_guides_csv_and_upload(ssh, guides, guides_csv)

    body = _build_offtarget_pbs_body(
        conda_env=conda_env,
        rscript_path=rscript_path,
        guides_csv=guides_csv,
        output_tsv=output_tsv,
        done_flag=done_flag,
        nuclease=nuclease,
        scores=scores,
        genome=genome,
        max_mismatches=max_mismatches,
    )
    spec = PBSJobSpec(
        name="offtarget",
        ncpus=4,
        mem_gb=16,
        walltime_hours=walltime_hours,
        log_path=log_path,
        body=body,
    )
    script_text = render_pbs_script(spec)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pbs", delete=False) as tmp:
        tmp.write(script_text)
        local_pbs = tmp.name
    try:
        ssh.upload(local_pbs, script_path)
    finally:
        os.unlink(local_pbs)

    job = PBSJob(spec=spec, remote_script_path=script_path)
    qsub(ssh, job)
    poll_until_done(ssh, [job], max_hours=walltime_hours)

    # Verify success via the RDS sentinel file (written by the PBS job body
    # immediately after the R script finishes, on the shared filesystem).
    out, _err, _rc = ssh.exec(f"test -f {done_flag} && echo 1 || echo 0")
    if out.strip() != "1":
        raise RuntimeError(f"off-target job did not write sentinel {done_flag}")

    # Download the TSV
    fd, local_tsv = tempfile.mkstemp(suffix=".tsv")
    os.close(fd)
    try:
        ssh.download(output_tsv, local_tsv)
        with open(local_tsv) as fh:
            tsv_text = fh.read()
    finally:
        os.unlink(local_tsv)
    return parse_offtarget_tsv(tsv_text)
