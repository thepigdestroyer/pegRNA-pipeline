"""Off-target scoring submitter.

Submits an ``Rscript scripts/score_guides.R`` PBS job, waits for it,
downloads the output TSV, and parses it into a list of
:class:`OffTargetReport`.
"""

from __future__ import annotations

import csv
import io
import tempfile
from dataclasses import dataclass
from typing import List, Sequence, Tuple


OFFTARGET_SCORING_DONE = "OFFTARGET_SCORING_DONE"


@dataclass(frozen=True)
class OffTargetReport:
    """One row of :file:`offtarget_scores.tsv`."""

    guide_id: str
    spacer: str
    deepspcas9: float
    deephf: float
    cfd_specificity: float
    mm0_count: int
    mm1_count: int
    mm2_count: int
    mm3_count: int


def _to_float(value: str) -> float:
    """Parse a float; treat empty / NA / NaN as 0.0."""
    if value is None:
        return 0.0
    v = value.strip()
    if v == "" or v.upper() in {"NA", "NAN", "NULL"}:
        return 0.0
    try:
        return float(v)
    except ValueError:
        return 0.0


def _to_int(value: str) -> int:
    v = (value or "").strip()
    if v == "" or v.upper() in {"NA", "NAN", "NULL"}:
        return 0
    try:
        return int(v)
    except ValueError:
        return 0


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
    nuclease: str,
    scores: Sequence[str],
    genome: str,
    max_mismatches: int,
) -> str:
    score_arg = ",".join(scores)
    return (
        "set -euo pipefail\n"
        "module load miniforge/3\n"
        f"conda activate {conda_env}\n"
        f"Rscript {rscript_path} "
        f"--guides {guides_csv} "
        f"--nuclease {nuclease} "
        f"--scores {score_arg} "
        f"--genome {genome} "
        f"--max-mismatches {max_mismatches} "
        f"--output {output_tsv}\n"
        f"echo {OFFTARGET_SCORING_DONE}\n"
    )


def _write_guides_csv_and_upload(ssh, guides, remote_path: str) -> None:
    """Write an in-memory guides CSV and upload it."""
    with tempfile.NamedTemporaryFile(
        mode="w", newline="", suffix=".csv", delete=False
    ) as tmp:
        writer = csv.writer(tmp)
        writer.writerow(["guide_id", "spacer"])
        for g in guides:
            writer.writerow([g["guide_id"], g["spacer"]])
        local = tmp.name
    ssh.upload(local, remote_path)


def submit_offtarget_job(
    *,
    ssh,
    guides: Sequence[dict],
    remote_run_dir: str,
    rscript_path: str,
    conda_env: str = "pegrna-ngrna",
    nuclease: str = "SpCas9_H840A",
    scores: Tuple[str, ...] = ("DeepSpCas9", "DeepHF", "CFD"),
    genome: str = "GRCh38",
    max_mismatches: int = 3,
    walltime_hours: int = 1,
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
    log_path = f"{remote_run_dir}/offtarget_pbs.log"
    script_path = f"{remote_run_dir}/offtarget.pbs"

    _write_guides_csv_and_upload(ssh, guides, guides_csv)

    body = _build_offtarget_pbs_body(
        conda_env=conda_env,
        rscript_path=rscript_path,
        guides_csv=guides_csv,
        output_tsv=output_tsv,
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
    # Heredoc upload
    ssh.exec(f"cat > {script_path} <<'PBS_EOF'\n{script_text}\nPBS_EOF")

    job = PBSJob(spec=spec, remote_script_path=script_path)
    qsub(ssh, job)
    poll_until_done(ssh, [job], max_hours=walltime_hours)

    # Verify success
    out, _err, _rc = ssh.exec(f"grep -c {OFFTARGET_SCORING_DONE} {log_path} || true")
    try:
        hits = int(out.strip().splitlines()[0]) if out.strip() else 0
    except ValueError:
        hits = 0
    if hits == 0:
        raise RuntimeError(f"off-target job did not emit {OFFTARGET_SCORING_DONE}")

    # Download the TSV
    with tempfile.NamedTemporaryFile(mode="r", suffix=".tsv", delete=False) as tmp:
        local_tsv = tmp.name
    ssh.download(output_tsv, local_tsv)
    with open(local_tsv) as fh:
        tsv_text = fh.read()
    return parse_offtarget_tsv(tsv_text)
