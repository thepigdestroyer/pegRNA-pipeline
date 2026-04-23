"""PBS Pro job rendering + submission + polling.

This module is side-effect free in isolation: every operation that would
touch a real cluster goes through an :class:`SSHClient` (or compatible
mock). The main entry points are:

* :class:`PBSJobSpec` / :class:`PBSJob` — frozen data containers
* :func:`render_pbs_script` — deterministic ``#PBS ...`` script renderer
* :func:`parse_qstat_state` — parse the first matching row's state
* :func:`qsub` / :func:`qstat_one` / :func:`qdel` — one-shot wrappers
* :func:`poll_until_done` — batched polling until every job is finished
* :func:`retry_with_doubled_resources` — walltime/OOM retry helper (T3)
* :func:`submit_with_concurrency_cap` — rate-limited batch submitter (T5)
"""

from __future__ import annotations

import dataclasses
import re
import time
from dataclasses import dataclass
from typing import Callable, List, Literal, Optional, Sequence


# ---------------------------------------------------------------- constants


#: Upper cap on walltime for retry_with_doubled_resources (Plan 3 T3).
WALLTIME_CAP_HOURS: int = 24

#: Upper cap on memory for retry_with_doubled_resources (Plan 3 T3).
MEM_CAP_GB: int = 64


# ---------------------------------------------------------------- datatypes


@dataclass(frozen=True)
class PBSJobSpec:
    """Static description of a PBS job (script body + resources)."""

    name: str
    ncpus: int
    mem_gb: int
    walltime_hours: int
    log_path: str
    body: str
    ngpus: int = 0
    queue: Optional[str] = None


@dataclass
class PBSJob:
    """Mutable runtime handle for a submitted job."""

    spec: PBSJobSpec
    remote_script_path: str
    scheduler_id: str = ""
    state: str = ""


# ---------------------------------------------------------------- render


def render_pbs_script(spec: PBSJobSpec) -> str:
    """Render a PBS Pro script for ``spec`` as a single string.

    Layout::

        #!/bin/bash
        #PBS -N <name>
        #PBS -l select=1:ncpus=N:mem=Ngb[:ngpus=N]
        #PBS -l walltime=HH:00:00
        #PBS -o <log_path>
        #PBS -j oe
        [#PBS -q <queue>]

        <body>
    """
    select = f"select=1:ncpus={spec.ncpus}:mem={spec.mem_gb}gb"
    if spec.ngpus > 0:
        select = f"{select}:ngpus={spec.ngpus}"
    lines: List[str] = [
        "#!/bin/bash",
        f"#PBS -N {spec.name}",
        f"#PBS -l {select}",
        f"#PBS -l walltime={spec.walltime_hours:02d}:00:00",
        f"#PBS -o {spec.log_path}",
        "#PBS -j oe",
    ]
    if spec.queue:
        lines.append(f"#PBS -q {spec.queue}")
    lines.append("")
    lines.append(spec.body.rstrip("\n"))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------- parse


_QSTAT_ROW = re.compile(
    r"^\s*(?P<id>\S+)\s+\S+\s+\S+\s+\S+\s+(?P<state>[A-Z])\s+\S+\s*$",
    re.MULTILINE,
)


def parse_qstat_state(qstat_output: str, job_id_prefix: str) -> str:
    """Return the 1-letter state ('R', 'Q', 'F', ...) for the first row
    whose job id starts with ``job_id_prefix``.

    Returns an empty string if the job id is not present (in PBS Pro a
    completed job disappears from qstat after a brief retention window).
    """
    for match in _QSTAT_ROW.finditer(qstat_output):
        jid = match.group("id")
        if jid.startswith(job_id_prefix):
            return match.group("state")
    return ""


# ---------------------------------------------------------------- commands


def qsub(ssh, job: PBSJob) -> str:
    """Submit ``job`` via ``qsub <remote_script_path>`` on ``ssh``.

    The first whitespace-delimited token of stdout is treated as the
    scheduler id (PBS Pro prints ``1234567.hostname\\n``). The id is
    stored on ``job.scheduler_id`` and also returned.
    """
    out, err, rc = ssh.exec(f"qsub {job.remote_script_path}")
    if rc != 0:
        raise RuntimeError(f"qsub failed: rc={rc} stderr={err!r}")
    sid = out.strip().split()[0] if out.strip() else ""
    if not sid:
        raise RuntimeError(f"qsub returned no id: stdout={out!r} stderr={err!r}")
    job.scheduler_id = sid
    return sid


def qstat_one(ssh, job: PBSJob) -> str:
    """Return the current state of ``job`` ('' if it has aged out)."""
    out, _err, _rc = ssh.exec(f"qstat {job.scheduler_id}")
    return parse_qstat_state(out, job.scheduler_id.split(".")[0])


def qdel(ssh, job: PBSJob) -> None:
    """Best-effort ``qdel <scheduler_id>`` — does not raise on non-zero rc."""
    ssh.exec(f"qdel {job.scheduler_id}")


# ---------------------------------------------------------------- poll


def poll_until_done(
    ssh,
    jobs: Sequence[PBSJob],
    *,
    interval_seconds: int = 60,
    max_hours: int = 24,
    sleep_fn: Callable[[float], None] = time.sleep,
) -> List[PBSJob]:
    """Block until every job in ``jobs`` has disappeared from qstat.

    One batched ``qstat <id1> <id2> ...`` call per iteration minimises
    PBS daemon load. Each job's ``state`` is updated in place; jobs that
    leave the queue are marked ``'F'``. Raises :class:`TimeoutError` if
    ``max_hours`` elapses (measured in polling wall-clock seconds).

    Parameters
    ----------
    sleep_fn:
        Injectable so tests can pass a no-op (otherwise a test suite
        would need to wait a real minute per poll).
    """
    job_list = list(jobs)
    # Short-circuit on empty input.
    if not job_list:
        return job_list

    max_seconds = max_hours * 3600
    elapsed = 0
    while True:
        ids = " ".join(j.scheduler_id for j in job_list if j.scheduler_id)
        out, _err, _rc = ssh.exec(f"qstat {ids}")
        pending = 0
        for job in job_list:
            prefix = job.scheduler_id.split(".")[0]
            state = parse_qstat_state(out, prefix)
            if state == "":
                job.state = "F"
            else:
                job.state = state
                pending += 1
        if pending == 0:
            return job_list
        if elapsed >= max_seconds:
            raise TimeoutError(
                f"poll_until_done: {pending} job(s) still pending after {max_hours}h"
            )
        sleep_fn(interval_seconds)
        elapsed += interval_seconds


# ---------------------------------------------------------------- T3: retry


def retry_with_doubled_resources(
    spec: PBSJobSpec,
    *,
    reason: Literal["walltime", "oom"],
) -> PBSJobSpec:
    """Return a new spec with the failing dimension doubled, capped.

    ``reason='walltime'`` doubles ``walltime_hours`` (capped at
    :data:`WALLTIME_CAP_HOURS`); ``reason='oom'`` doubles ``mem_gb``
    (capped at :data:`MEM_CAP_GB`). Any other reason raises :class:`ValueError`.
    """
    if reason == "walltime":
        new = min(spec.walltime_hours * 2, WALLTIME_CAP_HOURS)
        return dataclasses.replace(spec, walltime_hours=new)
    if reason == "oom":
        new = min(spec.mem_gb * 2, MEM_CAP_GB)
        return dataclasses.replace(spec, mem_gb=new)
    raise ValueError(f"unknown reason {reason!r}; expected 'walltime' or 'oom'")


# ---------------------------------------------------------------- T5: cap


def _default_active_count_fn(ssh) -> int:
    """Return the number of this user's pending+running jobs."""
    out, _err, _rc = ssh.exec("qstat -u $USER")
    # Very simple count: any non-header line that contains an R or Q state.
    count = 0
    for line in out.splitlines():
        if re.search(r"\s[RQ]\s", line):
            count += 1
    return count


def submit_with_concurrency_cap(
    ssh,
    jobs: Sequence[PBSJob],
    *,
    max_concurrent: int = 16,
    backoff_seconds: int = 120,
    sleep_fn: Callable[[float], None] = time.sleep,
    active_count_fn: Optional[Callable[[object], int]] = None,
) -> List[PBSJob]:
    """Submit each job, blocking when active >= ``max_concurrent``.

    The active-count probe is pluggable so tests can drive the cap
    deterministically without spinning up a fake qstat parser.
    """
    _count = active_count_fn or _default_active_count_fn
    submitted: List[PBSJob] = []
    for job in jobs:
        # Wait until active count is below the cap.
        while _count(ssh) >= max_concurrent:
            sleep_fn(backoff_seconds)
        qsub(ssh, job)
        submitted.append(job)
    return submitted
