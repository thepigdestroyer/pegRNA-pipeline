"""Scaffold-nCas9 docking gate (decision 7b) — final + heaviest gate.

Boltz-2 (MIT, AF3-class) predicts the full pegRNA + nCas9 complex; we compare the
scaffold↔Cas9 interface confidence (ipTM / interface-pLDDT, + optional
binding-affinity proxy) of the candidate against a cached apo baseline and flag
if docking degrades. ``ComplexDockScorer`` is a Protocol so Boltz-2 (GPU/HPC) or
a fake are interchangeable; the apo baseline is computed once and reused across
candidates. Threshold has no published pegRNA-Cas9 value → starts advisory,
calibrate on real data. Fallbacks: Chai-1 (contact restraints) > Protenix >
RoseTTAFold2NA.
"""

from __future__ import annotations

from typing import Optional, Protocol, Tuple

__all__ = ["ComplexDockScorer", "DockFake", "Boltz2DockScorer", "passes_docking_gate"]


class ComplexDockScorer(Protocol):
    def interface_iptm(self, pegrna: str, cas9: str) -> float: ...
    def apo_interface_iptm(self, cas9: str) -> float: ...


class DockFake:
    """Deterministic fake for CI/tests."""

    def __init__(self, apo_iptm: float, candidate_iptm: float) -> None:
        self._apo = apo_iptm
        self._cand = candidate_iptm

    def apo_interface_iptm(self, cas9: str) -> float:  # noqa: ARG002
        return self._apo

    def interface_iptm(self, pegrna: str, cas9: str) -> float:  # noqa: ARG002
        return self._cand


class Boltz2DockScorer:
    """Real Boltz-2 adapter (CX3-only; big GPU + MSA).

    Predicts the pegRNA+nCas9 complex and reads the scaffold↔Cas9 interface ipTM.
    Heavy deps load lazily so the core package + CI never import them; raises a
    helpful error off-GPU — never silently skipped. The apo baseline is cached by
    the orchestrator (Task 13) and reused across candidates.
    """

    def __init__(
        self,
        *,
        cache_dir: str | None = None,
        use_msa_server: bool = False,
        offline_single_sequence: bool = True,
        use_kernels: bool = False,
        apo_scaffold: str | None = None,
    ) -> None:
        self.cache_dir = cache_dir
        self.use_msa_server = use_msa_server
        # Boltz-2 (>=2.2) accelerates its triangular-multiplicative updates with the
        # cuequivariance_torch CUDA kernels, which it hard-imports only when kernels
        # are enabled. That package is not part of the base install and needs a
        # matching CUDA toolchain, so by default we pass --no_kernels: Boltz then
        # uses the pure-PyTorch path (slower, but portable and offline-safe). Without
        # this the docking gate fails with "No module named 'cuequivariance_torch'"
        # (the real GPU-smoke failure). Set use_kernels=True if the CUDA kernels are
        # installed and you want the speed-up.
        self.use_kernels = use_kernels
        # Boltz-2 requires an MSA per protein/RNA chain. On offline CX3 compute
        # nodes the MSA server is unreachable, so by default we run in
        # SINGLE-SEQUENCE mode (`msa: empty` per chain) — the documented offline
        # path. Without this, boltz aborts with "Missing MSA's in input and
        # --use_msa_server flag not set" (the real GPU-smoke failure). Single-
        # sequence is lower-accuracy but deterministic and self-contained; set
        # use_msa_server=True to instead generate MSAs (login node / online only).
        self.offline_single_sequence = offline_single_sequence and not use_msa_server
        # The apo baseline is the STANDARD scaffold + Cas9 complex (an RNA chain
        # must be present or ipTM is undefined). Defaults to the SpCas9 sgRNA
        # scaffold so the baseline is a real two-chain interface.
        if apo_scaffold is None:
            from pegrna_pipeline.hpc.tjpe_workflow import (  # noqa: PLC0415
                STANDARD_SGRNA_SCAFFOLD,
            )

            apo_scaffold = STANDARD_SGRNA_SCAFFOLD
        self.apo_scaffold = apo_scaffold
        self._apo_cache: dict = {}

    def _build_predict_cmd(self, yaml_path: str, out_dir: str) -> list:
        """Build the ``boltz predict`` command line for this scorer's settings.

        Adds ``--cache`` (offline weights), ``--no_kernels`` (skip the
        cuequivariance CUDA kernels unless ``use_kernels``), and ``--use_msa_server``
        only when explicitly requested. Extracted for testability.
        """
        cmd = ["boltz", "predict", str(yaml_path), "--out_dir", str(out_dir)]
        if self.cache_dir:
            cmd += ["--cache", str(self.cache_dir)]
        if not self.use_kernels:
            cmd.append("--no_kernels")
        if self.use_msa_server:
            cmd.append("--use_msa_server")
        return cmd

    @staticmethod
    def _build_complex_yaml(rna: str, protein: str, *, single_sequence: bool) -> str:
        """Build the Boltz ``predict`` input YAML for an (RNA pegRNA + protein Cas9)
        complex. In ``single_sequence`` mode each chain gets ``msa: empty`` so Boltz
        runs without an MSA server (required on offline compute nodes)."""
        msa_line = "\n      msa: empty" if single_sequence else ""
        seqs = []
        if rna:
            seqs.append(f"  - rna:\n      id: P\n      sequence: {rna}{msa_line}")
        seqs.append(f"  - protein:\n      id: C\n      sequence: {protein}{msa_line}")
        return "version: 1\nsequences:\n" + "\n".join(seqs) + "\n"

    def _predict(self, components):  # pragma: no cover - exercised only on big GPU
        """Run Boltz-2 on a pegRNA(+Cas9) complex via the documented ``boltz
        predict`` CLI and return the scaffold<->Cas9 interface ipTM.

        Faithful to the public Boltz interface: write a YAML describing the
        chains (RNA pegRNA + protein nCas9), invoke ``boltz predict``, then parse
        the emitted ``confidence_*.json`` (Boltz writes ``iptm`` / ``complex_iptm``
        and per-chain-pair ipTM). VALIDATE ON CX3 — exact JSON key names and the
        chain-pair index for the scaffold<->Cas9 interface must be confirmed
        against the installed Boltz version before a production run (see the
        on-CX3 smoke harness, scripts/hpc/smoke_models.py).
        """
        import json
        import shutil
        import subprocess
        import tempfile
        from pathlib import Path

        if shutil.which("boltz") is None:
            raise ImportError(
                "Boltz2DockScorer requires the '[complex]' extra (boltz + torch) and "
                "a large GPU + MSA; the 'boltz' CLI was not found. Run on CX3."
            )
        rna = components.get("pegrna")
        protein = components.get("cas9", "")
        with tempfile.TemporaryDirectory(dir=self.cache_dir) as tmp:
            tmp_path = Path(tmp)
            yaml_path = tmp_path / "complex.yaml"
            yaml_path.write_text(
                self._build_complex_yaml(rna, protein, single_sequence=self.offline_single_sequence)
            )
            cmd = self._build_predict_cmd(str(yaml_path), str(tmp_path))
            proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
            if proc.returncode != 0:
                tail = (proc.stderr or proc.stdout or "").strip()[-800:]
                raise RuntimeError(
                    f"boltz predict failed (rc={proc.returncode}). "
                    f"cache={self.cache_dir!r}. stderr/stdout tail:\n{tail}"
                )
            # Boltz writes to <out_dir>/boltz_results_<stem>/predictions/<stem>/;
            # search broadly for the confidence JSON (name varies by version).
            conf_files = sorted(tmp_path.rglob("confidence_*.json")) or sorted(
                tmp_path.rglob("*confidence*.json")
            )
            if not conf_files:
                produced = [str(p.relative_to(tmp_path)) for p in tmp_path.rglob("*")][:40]
                tail = (proc.stderr or proc.stdout or "").strip()[-500:]
                raise RuntimeError(
                    "Boltz produced no confidence JSON. "
                    f"cache={self.cache_dir!r}; files under out_dir: {produced}; "
                    f"boltz output tail:\n{tail}"
                )
            data = json.loads(conf_files[0].read_text())
            # Prefer the chain-pair ipTM for the RNA<->protein interface; fall
            # back to the global complex ipTM. (Key names: confirm on CX3.)
            pair = data.get("pair_chains_iptm", {})
            try:
                return float(pair["P"]["C"])
            except (KeyError, TypeError):
                return float(data.get("iptm", data.get("complex_iptm", 0.0)))

    def apo_interface_iptm(self, cas9: str) -> float:  # pragma: no cover - CX3 only
        if cas9 in self._apo_cache:
            return self._apo_cache[cas9]
        # Baseline = scaffold-only pegRNA (no extension) + Cas9 — a real two-chain
        # interface, so the candidate's interface drop is meaningful.
        val = self._predict({"pegrna": self.apo_scaffold, "cas9": cas9})
        self._apo_cache[cas9] = val
        return val

    def interface_iptm(self, pegrna: str, cas9: str) -> float:  # pragma: no cover - CX3
        return self._predict({"pegrna": pegrna, "cas9": cas9})


def passes_docking_gate(
    pegrna: str,
    cas9: str,
    *,
    scorer: ComplexDockScorer,
    max_drop: float = 0.2,
    min_iptm: float = 0.6,
    metrics: Optional[dict] = None,
) -> Tuple[bool, str]:
    """HARD GATE on scaffold-Cas9 docking (interface ipTM).

    Rejects if the candidate interface ipTM either (a) falls below ``min_iptm``
    (AlphaFold-Multimer/AF3 "likely incorrect" floor, 0.6) or (b) drops by more
    than ``max_drop`` (one full confidence band, 0.2) relative to the apo
    scaffold-Cas9 baseline. Defaults are literature-derived (see
    gate_thresholds.THRESHOLDS: Evans 2022, Abramson 2024).

    If ``metrics`` is supplied it is populated with the raw values and cutoffs
    (interface ipTM vs the min-ipTM floor, plus apo baseline / drop vs max-drop)
    for the caller's progress display; it never changes the return value.
    """
    apo = scorer.apo_interface_iptm(cas9)
    cand = scorer.interface_iptm(pegrna, cas9)
    drop = apo - cand
    if metrics is not None:
        metrics.update({
            "metric": "interface ipTM",
            "value": round(cand, 3),
            "cutoff": min_iptm,
            "cmp": ">=",
            "extra": f"apo {apo:.2f}→{cand:.2f}, drop {drop:.2f} (max {max_drop})",
        })
    if cand < min_iptm - 1e-9:
        return False, (
            f"scaffold-Cas9 interface ipTM {cand:.2f} < {min_iptm} "
            "(AF-multimer 'likely incorrect' floor)"
        )
    if drop > max_drop + 1e-9:  # epsilon guards float artefacts at the exact boundary
        return False, (
            f"scaffold-Cas9 interface ipTM dropped {drop:.2f} "
            f"(apo {apo:.2f} -> {cand:.2f} > {max_drop})"
        )
    return True, f"interface preserved (apo {apo:.2f} -> {cand:.2f})"
