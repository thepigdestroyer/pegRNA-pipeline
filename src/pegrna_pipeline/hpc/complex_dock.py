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
        try:
            import boltz  # noqa: F401, PLC0415
            import torch  # noqa: F401, PLC0415
        except Exception as exc:
            raise ImportError(
                "Boltz2DockScorer requires the '[complex]' extra (boltz + torch) and "
                "a large GPU + MSA; run on CX3. "
                f"Original error: {exc}"
            ) from exc
        raise NotImplementedError(
            "Boltz-2 complex prediction + interface-ipTM extraction is wired in the "
            "follow-on plan"
        )

    def apo_interface_iptm(self, cas9: str) -> float:  # pragma: no cover - CX3 only
        if cas9 in self._apo_cache:
            return self._apo_cache[cas9]
        return self._predict({"cas9": cas9})

    def interface_iptm(self, pegrna: str, cas9: str) -> float:  # pragma: no cover - CX3
        return self._predict({"pegrna": pegrna, "cas9": cas9})


def passes_docking_gate(
    pegrna: str, cas9: str, *, scorer: ComplexDockScorer, max_drop: float = 0.15
) -> Tuple[bool, str]:
    """HARD GATE: reject if the scaffold-Cas9 interface ipTM drops by > ``max_drop``
    relative to the apo baseline."""
    apo = scorer.apo_interface_iptm(cas9)
    cand = scorer.interface_iptm(pegrna, cas9)
    drop = apo - cand
    if drop > max_drop + 1e-9:  # epsilon guards float artefacts at the exact boundary
        return False, (
            f"scaffold-Cas9 interface ipTM dropped {drop:.2f} "
            f"(apo {apo:.2f} -> {cand:.2f} > {max_drop})"
        )
    return True, f"interface preserved (apo {apo:.2f} -> {cand:.2f})"
