"""Laptop-side launcher for pegRNA-pipeline HPC runs.

:class:`HPCLauncher` orchestrates:

1. Login-node gene lookup via ``python -m pegrna_pipeline.hpc.sequence_fetcher``.
2. Creation of the remote run directory + ``inputs.json`` + ``construct.yaml``.
3. Submission of the driver PBS job via ``qsub -v RUN_DIR=...``.
4. Polling ``driver.log`` for the ``DRIVER_DONE`` marker.
5. Downloading the ``reports/`` directory to a local path.

All SSH interactions go through a :class:`SSHClient`-compatible object
produced by the ``ssh_factory`` — tests pass a lightweight fake; the
production default constructs a real :class:`SSHClient` from
``config['hostname']`` / ``config['username']``.
"""

from __future__ import annotations

import dataclasses
import json
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Optional

import yaml

from pegrna_pipeline.hpc.data_types import GeneRecord


# ---------------------------------------------------------------- result
@dataclass
class RunResult:
    """Return value of :meth:`HPCLauncher.launch_run`."""

    run_id: str
    remote_run_dir: str
    local_reports_dir: str
    driver_job_id: str
    success: bool


# ---------------------------------------------------------------- launcher
def _default_ssh_factory(*, hostname, username, password_provider, **kw):
    """Production ``ssh_factory`` default."""
    from pegrna_pipeline.hpc.ssh_client import SSHClient

    return SSHClient(
        hostname=hostname,
        username=username,
        password_provider=password_provider,
        **kw,
    )


class HPCLauncher:
    """Orchestrates a single end-to-end pegRNA-pipeline run."""

    def __init__(
        self,
        *,
        config: Dict[str, Any],
        password_provider: Optional[Callable[[], str]] = None,
        ssh_factory: Optional[Callable[..., Any]] = None,
    ) -> None:
        self.config = config
        self.password_provider = password_provider or (lambda: "")
        self.ssh_factory = ssh_factory or _default_ssh_factory

    # ---------------------------------------------------------- helpers
    def _build_ssh(self):
        """Build an SSH client per ``self.config`` and the injected factory."""
        return self.ssh_factory(
            hostname=self.config["hostname"],
            username=self.config["username"],
            password_provider=self.password_provider,
        )

    def _fetch_gene_record_on_login(self, ssh, gene: str) -> GeneRecord:
        """Run ``sequence_fetcher --json`` on the login node, parse the record."""
        repo_dir = self.config.get("repo_dir", "~/pegRNA-pipeline")
        email = self.config.get("entrez_email", "user@example.com")
        # Bootstrap conda — ssh exec runs in a non-login shell on Imperial
        # CX3, so neither the miniforge module nor ~/miniforge3 is on PATH
        # by default. See setup_hpc.sh / driver_pbs_wrapper.sh for the same
        # pattern.
        cmd = (
            "module load miniforge/3 >/dev/null 2>&1 || true; "
            'source "$HOME/miniforge3/etc/profile.d/conda.sh" && '
            "conda activate pegrna-pipeline && "
            f"cd {repo_dir} && "
            f"python -m pegrna_pipeline.hpc.sequence_fetcher "
            f"{gene} --email {email} --json"
        )
        out, err, rc = ssh.exec(cmd)
        if rc != 0:
            raise RuntimeError(f"login-node sequence_fetcher failed: rc={rc} err={err!r}")
        data = json.loads(out)
        return GeneRecord(**data)

    def _upload_json(self, ssh, remote_path: str, payload: Dict[str, Any]) -> None:
        """Upload ``payload`` as pretty-printed JSON via a heredoc."""
        text = json.dumps(payload, indent=2)
        escaped = text.replace("$", "\\$")
        ssh.exec(f"cat > {remote_path} <<'JSON_EOF'\n{escaped}\nJSON_EOF")

    def _upload_yaml(self, ssh, remote_path: str, payload: Dict[str, Any]) -> None:
        text = yaml.safe_dump(payload, sort_keys=False)
        ssh.exec(f"cat > {remote_path} <<'YAML_EOF'\n{text}\nYAML_EOF")

    def _qsub_driver(self, ssh, remote_run_dir: str) -> str:
        """qsub the wrapper with -v RUN_DIR=... and return the scheduler id."""
        repo_dir = self.config.get("repo_dir", "~/pegRNA-pipeline")
        wrapper = f"{repo_dir}/scripts/driver_pbs_wrapper.sh"
        out, err, rc = ssh.exec(f"/opt/pbs/bin/qsub -v RUN_DIR={remote_run_dir} {wrapper}")
        if rc != 0:
            raise RuntimeError(f"qsub failed: rc={rc} err={err!r}")
        sid = out.strip().split()[0] if out.strip() else ""
        if not sid:
            raise RuntimeError(f"qsub returned no id: {out!r}")
        return sid

    def _poll_driver(
        self,
        ssh,
        driver_id: str,
        remote_run_dir: str,
        interval: int = 60,
        max_hours: int = 48,
        sleep_fn: Callable[[float], None] = time.sleep,
    ) -> bool:
        """Poll until the driver's log emits ``DRIVER_DONE``.

        Per v1.1 corrections: greps ``<remote_run_dir>/driver.log`` — never
        a wildcard — to avoid colliding with other runs.
        """
        driver_log = f"{remote_run_dir}/driver.log"
        max_seconds = max_hours * 3600
        elapsed = 0
        while True:
            out, _err, _rc = ssh.exec(f"grep -c DRIVER_DONE {driver_log} || true")
            try:
                hits = int(out.strip().splitlines()[0]) if out.strip() else 0
            except ValueError:
                hits = 0
            if hits > 0:
                return True
            # Check if the job is still in the queue; if not, something went wrong.
            qstat_out, _q_err, _q_rc = ssh.exec(f"/opt/pbs/bin/qstat {driver_id}")
            if elapsed >= max_seconds:
                raise TimeoutError(f"driver {driver_id} not done after {max_hours}h")
            sleep_fn(interval)
            elapsed += interval

    def download_reports(self, *, remote_run_dir: str, local_reports_dir: str) -> None:
        """Download every file under ``remote_run_dir/reports/`` to local.

        Public per v1.1 (used by ``pegrna fetch``).
        """
        local = Path(local_reports_dir)
        local.mkdir(parents=True, exist_ok=True)
        with self._build_ssh() as ssh:
            # List files
            out, _err, _rc = ssh.exec(f"find {remote_run_dir}/reports -type f")
            for line in out.splitlines():
                remote_path = line.strip()
                if not remote_path:
                    continue
                name = Path(remote_path).name
                ssh.download(remote_path, str(local / name))

    # ---------------------------------------------------------- main entry
    def launch_run(
        self,
        *,
        gene: str,
        construct: str,
        mode: str,
        local_output_dir: str,
        cleanup: bool = False,
    ) -> RunResult:
        """Execute the full launch cycle and return a :class:`RunResult`."""
        run_id = f"{gene}_{mode}_{uuid.uuid4().hex[:8]}"
        remote_base = self.config.get("remote_runs_base", "/work/runs")
        remote_run_dir = f"{remote_base}/{run_id}"
        local_reports_dir = str(Path(local_output_dir) / run_id / "reports")

        with self._build_ssh() as ssh:
            # 1. Make remote directory + reports subdir.
            ssh.exec(f"mkdir -p {remote_run_dir}/reports")

            # 2. Fetch gene record on the login node.
            gene_record = self._fetch_gene_record_on_login(ssh, gene)

            # 3. Resolve the construct on the laptop side.
            resolved = self._resolve_construct_local(construct)

            # 4. Build + upload inputs.json.
            inputs = {
                "gene": dataclasses.asdict(gene_record),
                "mode": mode,
                "construct": construct,
                "cell_line": self.config.get("cell_line", "KOLF2_iPSC"),
                "pridict_install_path": self.config["pridict_install_path"],
                "rscript_path": self.config["rscript_path"],
                "conda_env": self.config.get("conda_env", "pridict2"),
                "compute_mode": self.config.get("compute_mode", "local"),
                "back_ssh": self.config.get("back_ssh", {}),
            }
            self._upload_json(ssh, f"{remote_run_dir}/inputs.json", inputs)

            # 5. Upload resolved construct as YAML.
            self._upload_yaml(
                ssh,
                f"{remote_run_dir}/construct.yaml",
                resolved,
            )

            # 6. qsub the driver wrapper.
            driver_id = self._qsub_driver(ssh, remote_run_dir)

            # 7. Poll for DRIVER_DONE.
            success = self._poll_driver(ssh, driver_id, remote_run_dir, interval=60)

        # 8. Download reports.
        self.download_reports(
            remote_run_dir=remote_run_dir,
            local_reports_dir=local_reports_dir,
        )

        # Optional cleanup (best-effort).
        if cleanup:
            try:
                with self._build_ssh() as ssh:
                    ssh.exec(f"rm -rf {remote_run_dir}")
            except Exception:  # pragma: no cover — best-effort only
                pass

        return RunResult(
            run_id=run_id,
            remote_run_dir=remote_run_dir,
            local_reports_dir=local_reports_dir,
            driver_job_id=driver_id,
            success=success,
        )

    # ---------------------------------------------------------- construct
    def _resolve_construct_local(self, construct: str) -> Dict[str, str]:
        """Resolve a construct token on the laptop using Plan 1's resolver.

        Returns a plain dict so it serialises nicely to YAML.
        """
        from pegrna_pipeline.agent.construct_parser import (
            parse_components,
            resolve_construct,
            tokenize,
        )
        from pegrna_pipeline.agent.motif_library import MotifLibrary

        library = MotifLibrary.load_full()
        comps = parse_components(tokenize(construct))
        resolved = resolve_construct(comps, library=library)
        return {
            "aa": resolved.aa,
            "dna": resolved.dna,
            "components": [
                {"name": c.name, "aa": c.aa, "dna": c.dna, "source": c.source}
                for c in resolved.components
            ],
        }
