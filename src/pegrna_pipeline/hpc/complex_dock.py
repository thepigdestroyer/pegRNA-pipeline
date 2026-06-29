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

from typing import Protocol, Tuple

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

    def __init__(self, *, cache_dir: str | None = None, use_msa_server: bool = False) -> None:
        self.cache_dir = cache_dir
        self.use_msa_server = use_msa_server
        self._apo_cache: dict = {}

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
            seqs = []
            if rna:
                seqs.append(f"  - rna:\n      id: P\n      sequence: {rna}")
            seqs.append(f"  - protein:\n      id: C\n      sequence: {protein}")
            yaml_path.write_text("version: 1\nsequences:\n" + "\n".join(seqs) + "\n")
            cmd = ["boltz", "predict", str(yaml_path), "--out_dir", str(tmp_path)]
            if self.use_msa_server:
                cmd.append("--use_msa_server")
            subprocess.run(cmd, check=True, capture_output=True)
            conf_files = list(tmp_path.rglob("confidence_*.json"))
            if not conf_files:
                raise RuntimeError("Boltz produced no confidence JSON")
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
        return self._predict({"cas9": cas9})

    def interface_iptm(self, pegrna: str, cas9: str) -> float:  # pragma: no cover - CX3
        return self._predict({"pegrna": pegrna, "cas9": cas9})


def passes_docking_gate(
    pegrna: str,
    cas9: str,
    *,
    scorer: ComplexDockScorer,
    max_drop: float = 0.2,
    min_iptm: float = 0.6,
) -> Tuple[bool, str]:
    """HARD GATE on scaffold-Cas9 docking (interface ipTM).

    Rejects if the candidate interface ipTM either (a) falls below ``min_iptm``
    (AlphaFold-Multimer/AF3 "likely incorrect" floor, 0.6) or (b) drops by more
    than ``max_drop`` (one full confidence band, 0.2) relative to the apo
    scaffold-Cas9 baseline. Defaults are literature-derived (see
    gate_thresholds.THRESHOLDS: Evans 2022, Abramson 2024).
    """
    apo = scorer.apo_interface_iptm(cas9)
    cand = scorer.interface_iptm(pegrna, cas9)
    drop = apo - cand
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
