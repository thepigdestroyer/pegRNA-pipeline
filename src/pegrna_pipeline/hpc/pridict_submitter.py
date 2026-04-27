"""PRIDICT2 batch job builder + submitter.

This module renders the PBS body that runs PRIDICT 2.0 in ``batch`` mode on
one input CSV chunk, and exposes a higher-level :func:`submit_pridict_batches`
that partitions a large variant set, uploads each chunk, qsubs with the
concurrency cap, waits for all jobs to finish, and reports the output
directories.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Sequence

PRIDICT2_BATCH_DONE = "PRIDICT2_BATCH_DONE"


def build_pridict_pbs_body(
    *,
    pridict_install_path: str,
    conda_env: str,
    input_csv_path: str,
    output_dir: str,
    cores: int = 3,
) -> str:
    """Return the bash body (no ``#PBS`` directives) for one PRIDICT chunk.

    The caller wraps this in a :func:`render_pbs_script` call.
    """
    return (
        "set -euo pipefail\n"
        "module load miniforge/3\n"
        # Imperial CX3 ships only a miniforge bootstrap; conda lives in
        # ~/miniforge3 (created by setup_hpc.sh).
        'source "$HOME/miniforge3/etc/profile.d/conda.sh"\n'
        f"conda activate {conda_env}\n"
        f"cd {pridict_install_path}\n"
        f"python pridict2_pegRNA_design.py batch "
        f"--input-fname {input_csv_path} --output-dir {output_dir} "
        f"--nicking --cores {cores}\n"
        f"echo {PRIDICT2_BATCH_DONE}\n"
    )


# ------------------------------------------------------------------ T7
@dataclass
class BatchResult:
    """Outcome of :func:`submit_pridict_batches`."""

    submitted: List[object] = field(default_factory=list)
    failed: List[object] = field(default_factory=list)
    output_dirs: List[str] = field(default_factory=list)


def _write_chunk_csv_and_upload(ssh, variants, remote_csv_path: str) -> None:
    """Serialise variants as a PRIDICT input CSV locally then upload.

    The CSV schema (PRIDICT 2.0 batch mode) has two columns:
    ``sequence_name,editseq``. ``editseq`` is a bracket-notation edit
    string produced by :func:`format_bracket_insertion` upstream.

    This helper exists so tests can monkey-patch it; the fixture does
    not need a real filesystem.
    """
    import csv
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", newline="", suffix=".csv", delete=False) as tmp:
        writer = csv.writer(tmp)
        writer.writerow(["sequence_name", "editseq"])
        for v in variants:
            writer.writerow([v["sequence_name"], v["editseq"]])
        local_path = tmp.name
    ssh.upload(local_path, remote_csv_path)


def submit_pridict_batches(
    *,
    ssh,
    variants: Sequence,
    remote_run_dir: str,
    pridict_install_path: str,
    conda_env: str = "pridict2",
    max_per_chunk: int = 500_000,
    max_concurrent: int = 16,
    walltime_hours: int = 12,
    ncpus: int = 8,
    mem_gb: int = 32,
    cores: int = 3,
) -> BatchResult:
    """Partition ``variants`` into chunks, qsub each, wait, and verify.

    Each variant is a dict with at least ``{"sequence_name", "editseq"}``.
    The function:

    1. Partitions via :func:`partition` (T4).
    2. Uploads each chunk CSV to ``remote_run_dir/chunk_<i>.csv``.
    3. Renders & uploads the ``.pbs`` script.
    4. Submits every job through :func:`submit_with_concurrency_cap` (T5).
    5. Polls until finished via :func:`poll_until_done`.
    6. Verifies each chunk's log contains ``PRIDICT2_BATCH_DONE``.
    """
    # Deferred imports to avoid cycles.
    from pegrna_pipeline.hpc.chunking import partition
    from pegrna_pipeline.hpc.pbs import (
        PBSJob,
        PBSJobSpec,
        poll_until_done,
        render_pbs_script,
        submit_with_concurrency_cap,
    )

    chunks = list(partition(variants, max_per_chunk=max_per_chunk))
    jobs: List[PBSJob] = []
    output_dirs: List[str] = []

    # Stage per-chunk resources on the remote side.
    for i, chunk in enumerate(chunks):
        csv_path = f"{remote_run_dir}/chunk_{i}.csv"
        out_dir = f"{remote_run_dir}/out_{i}"
        log_path = f"{remote_run_dir}/chunk_{i}_pbs.log"
        script_path = f"{remote_run_dir}/chunk_{i}.pbs"

        _write_chunk_csv_and_upload(ssh, chunk, csv_path)
        ssh.exec(f"mkdir -p {out_dir}")

        body = build_pridict_pbs_body(
            pridict_install_path=pridict_install_path,
            conda_env=conda_env,
            input_csv_path=csv_path,
            output_dir=out_dir,
            cores=cores,
        )
        spec = PBSJobSpec(
            name=f"pridict_{i}",
            ncpus=ncpus,
            mem_gb=mem_gb,
            walltime_hours=walltime_hours,
            log_path=log_path,
            body=body,
        )
        script_text = render_pbs_script(spec)
        # Upload script via a heredoc exec; PBS scripts are tiny.
        escaped = script_text.replace("'", "'\\''")
        ssh.exec(f"cat > {script_path} <<'PBS_EOF'\n{escaped}\nPBS_EOF")
        jobs.append(PBSJob(spec=spec, remote_script_path=script_path))
        output_dirs.append(out_dir)

    # Fire jobs respecting the concurrency cap.
    submit_with_concurrency_cap(ssh, jobs, max_concurrent=max_concurrent)

    # Wait for all of them.
    poll_until_done(ssh, jobs, max_hours=walltime_hours)

    # Verify success by grepping the log for the marker.
    submitted: List[PBSJob] = []
    failed: List[PBSJob] = []
    for i, job in enumerate(jobs):
        log_path = f"{remote_run_dir}/chunk_{i}_pbs.log"
        out, _err, _rc = ssh.exec(f"grep -c {PRIDICT2_BATCH_DONE} {log_path} || true")
        try:
            hits = int(out.strip().splitlines()[0]) if out.strip() else 0
        except ValueError:
            hits = 0
        if hits > 0:
            submitted.append(job)
        else:
            failed.append(job)

    return BatchResult(submitted=submitted, failed=failed, output_dirs=output_dirs)
