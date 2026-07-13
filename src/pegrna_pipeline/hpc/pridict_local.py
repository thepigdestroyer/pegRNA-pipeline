"""Local (subprocess) PRIDICT2 runner for the DSI-box CPU path.

Unlike :func:`pegrna_pipeline.hpc.pridict_submitter.submit_pridict_batches` (which
qsubs onto Imperial CX3 over SSH/PBS), :class:`PridictLocalRunner` shells out to a
**local** ``python pridict2_pegRNA_design.py batch`` on the DSI box where PRIDICT2
is installed on CPU. It mirrors the ``run_batch -> {name: {pegrnas, nicking}}`` shape
so callers can swap the CX3 and local paths freely.

Two PRIDICT2 CLI quirks are honoured (verified on the DSI box; wiring spec section 3):

* **Quirk 1** — ``--input-fname`` is resolved relative to ``PRIDICT2/input/``, so the
  input CSV is written into ``<PRIDICT2_DIR>/input/`` and only its *basename* is passed.
* **Quirk 2** — ``--output-dir`` is relative to the repo dir and must **pre-exist**
  (PRIDICT does not create it), so the runner ``mkdir -p``s it before launching.

Failure is loud: a non-zero exit or a timeout raises :class:`RuntimeError` carrying the
captured stderr tail — there is no silent fallback. The temporary input CSV is always
removed afterwards. :class:`PridictRunnerFake` is a deterministic, subprocess-free
stand-in for CI, shaped exactly like the real parser output.
"""

from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path
from uuid import uuid4

from pegrna_pipeline.hpc.data_types import PegRNADesign
from pegrna_pipeline.hpc.pridict_parser import parse_pegrna_csv

__all__ = ["PridictLocalRunner", "PridictRunnerFake"]

#: PRIDICT2 batch entry point (run from the repo dir with ``cwd=PRIDICT2_DIR``).
_SCRIPT_NAME = "pridict2_pegRNA_design.py"
_DEFAULT_CORES = 8
_DEFAULT_TIMEOUT = 1800  # seconds (deep model ~10s/seq on CPU; batches amortise)
_STDERR_TAIL_CHARS = 4000


def _stderr_tail(value: object, limit: int = _STDERR_TAIL_CHARS) -> str:
    """Return the last ``limit`` characters of a str/bytes stderr blob (or "")."""
    if value is None:
        return ""
    if isinstance(value, bytes):
        value = value.decode("utf-8", "replace")
    return str(value)[-limit:]


class PridictLocalRunner:
    """Run PRIDICT2 as a local subprocess and parse its two output CSVs.

    Config comes from explicit constructor args, falling back to the environment
    (``PRIDICT2_DIR``, ``PRIDICT2_PYTHON``, ``PRIDICT2_CORES``). A missing dir or
    python path fails loud at construction rather than mid-run.
    """

    def __init__(
        self,
        *,
        pridict_dir: str | None = None,
        pridict_python: str | None = None,
        cores: int | None = None,
        timeout: int = _DEFAULT_TIMEOUT,
        script_name: str = _SCRIPT_NAME,
    ) -> None:
        self.pridict_dir = pridict_dir or os.environ.get("PRIDICT2_DIR")
        self.pridict_python = pridict_python or os.environ.get("PRIDICT2_PYTHON")
        if not self.pridict_dir:
            raise ValueError(
                "PRIDICT2_DIR is required (pass pridict_dir= or set the $PRIDICT2_DIR env var)"
            )
        if not self.pridict_python:
            raise ValueError(
                "PRIDICT2_PYTHON is required (pass pridict_python= or set $PRIDICT2_PYTHON)"
            )
        if cores is None:
            env_cores = os.environ.get("PRIDICT2_CORES")
            self.cores = int(env_cores) if env_cores else _DEFAULT_CORES
        else:
            self.cores = int(cores)
        self.timeout = timeout
        self.script_name = script_name

    def run_batch(self, rows: list[dict]) -> dict[str, dict]:
        """Design pegRNAs + nicking guides for every row in one PRIDICT batch.

        Each row is ``{"sequence_name", "editseq"}``. Returns a mapping
        ``sequence_name -> {"pegrnas": list[PegRNADesign], "nicking": list[dict]}``.
        """
        if not rows:
            raise ValueError("run_batch requires at least one row")

        uid = uuid4().hex
        repo = Path(self.pridict_dir)
        input_csv = repo / "input" / f"{uid}.csv"
        reldir = f"predictions/{uid}"  # repo-relative (quirk 2); POSIX for the Linux CLI
        outdir = repo / reldir

        try:
            input_csv.parent.mkdir(parents=True, exist_ok=True)
            self._write_input_csv(input_csv, rows)
            outdir.mkdir(parents=True, exist_ok=True)  # PRIDICT won't create it (quirk 2)
            # Quirk 1: --input-fname is basename-only, resolved under PRIDICT2/input/.
            self._run_subprocess(input_csv.name, reldir)
            return self._collect(outdir, rows)
        finally:
            input_csv.unlink(missing_ok=True)

    @staticmethod
    def _write_input_csv(path: Path, rows: list[dict]) -> None:
        with path.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["sequence_name", "editseq"])
            for row in rows:
                writer.writerow([row["sequence_name"], row["editseq"]])

    def _run_subprocess(self, input_basename: str, reldir: str) -> None:
        cmd = [
            self.pridict_python,
            self.script_name,
            "batch",
            "--input-fname",
            input_basename,
            "--output-dir",
            reldir,
            "--nicking",
            "--cores",
            str(self.cores),
        ]
        try:
            proc = subprocess.run(
                cmd,
                cwd=self.pridict_dir,
                capture_output=True,
                text=True,
                timeout=self.timeout,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"PRIDICT2 batch timed out after {self.timeout}s; "
                f"stderr tail:\n{_stderr_tail(exc.stderr)}"
            ) from exc
        if proc.returncode != 0:
            raise RuntimeError(
                f"PRIDICT2 batch failed (exit {proc.returncode}); "
                f"stderr tail:\n{_stderr_tail(proc.stderr)}"
            )

    def _collect(self, outdir: Path, rows: list[dict]) -> dict[str, dict]:
        result: dict[str, dict] = {}
        for row in rows:
            name = row["sequence_name"]
            peg_csv = outdir / f"{name}_pegRNA_Pridict_full.csv"
            if not peg_csv.is_file():
                # PRIDICT can exit 0 yet write no pegRNA CSV when it skips a row (e.g.
                # an insertion > 40 bp, or < 100 bp flanking each side). Fail loud with an
                # actionable message instead of the bare FileNotFoundError from the parser.
                raise RuntimeError(
                    f"PRIDICT2 produced no output for sequence {name!r} "
                    f"(expected {peg_csv.name}). PRIDICT ran but wrote no pegRNA CSV -- the "
                    "edit is likely outside its design domain (insertions must be <= 40 bp) "
                    "or the editseq had < 100 bp flanking sequence on each side."
                )
            pegrnas = parse_pegrna_csv(peg_csv)
            nicking = self._read_nicking(outdir / f"{name}_nicking_guides.csv")
            result[name] = {"pegrnas": pegrnas, "nicking": nicking}
        return result

    @staticmethod
    def _read_nicking(path: Path) -> list[dict]:
        """Read a ``*_nicking_guides.csv`` into raw ``csv.DictReader`` rows.

        Returns ``[]`` when the file is absent (PRIDICT can emit a pegRNA with no
        admissible nicking guide); the pegRNA CSV is the authoritative output and its
        absence surfaces loudly via :func:`parse_pegrna_csv`.
        """
        if not path.is_file():
            return []
        with path.open(newline="") as fh:
            return list(csv.DictReader(fh))


class PridictRunnerFake:
    """Deterministic, subprocess-free PRIDICT2 stand-in for CI.

    Shares :meth:`run_batch`'s signature with :class:`PridictLocalRunner` and returns
    the same shape (one :class:`PegRNADesign` + one ``csv.DictReader``-style nicking
    dict per row), so tests exercise the pipeline wiring without launching PRIDICT2.
    """

    def run_batch(self, rows: list[dict]) -> dict[str, dict]:
        result: dict[str, dict] = {}
        for row in rows:
            name = row["sequence_name"]
            peg = PegRNADesign(
                sequence_name=f"{name}_pegRNA1",
                spacer="GAATTCGTGTTCCTAAGCGT",
                pbs="TCGACGCTTAGGAACAC",
                rt="ACGTACGTACGTACGTACGTACGTACGTAC",
                pbs_length=13,
                rt_length=30,
                target_strand="Fw",
                editing_position=50,
                pe2_hek293t=0.682,
                pe2_k562=0.682,
                pegrna_full="",
            )
            nicking = [
                {
                    "Nicking-Protospacer": f"{name}_ng1",
                    "Nicking-Position-to-edit": "-30",
                    "Nicking-PAMdisrupt": "Nicking_PAM_disrupt",
                    "Target_Strand": "Fw",
                    "DeepCas9score": "0.55",
                }
            ]
            result[name] = {"pegrnas": [peg], "nicking": nicking}
        return result
