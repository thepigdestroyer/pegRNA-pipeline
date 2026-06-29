"""Factory that assembles the gate scorers for ``run_recoded_tjpe``.

Two modes:

* ``fake`` — deterministic in-repo fakes (CI, dry-run, smoke tests). No heavy
  deps, no GPU, no weights.
* ``real`` — the real model adapters (SpliceAI/Pangolin, AlphaFold-DB metrics,
  RiNALMo full structure, Boltz-2 docking, Cas-OFFinder off-target). These
  require the optional extras + staged weights + (mostly) a GPU; they raise a
  helpful error if those are missing. **There is no silent fallback to fakes in
  real mode** — a missing dependency fails loud so a production run never
  pretends a gate ran.

Returned bundle is a :class:`GateScorers` that maps 1:1 onto the keyword scorer
arguments of :func:`pegrna_pipeline.hpc.recoded_tjpe_workflow.run_recoded_tjpe`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Tuple

__all__ = ["GateScorers", "make_scorers"]


@dataclass
class GateScorers:
    splice_scorer: Any
    full_structure_scorer: Any
    importance_scorer: Any
    dock_scorer: Any
    offtarget_scorer: Callable[[List[str]], Dict[str, Tuple[float, int, bool]]]


def _fake_bundle() -> GateScorers:
    from pegrna_pipeline.hpc.complex_dock import DockFake
    from pegrna_pipeline.hpc.full_structure import FullStructureFake
    from pegrna_pipeline.hpc.splice_screen import SpliceFake
    from pegrna_pipeline.hpc.structure_importance import ResidueMetrics

    class _AllPassImportance:
        def metrics(self, index: int) -> ResidueMetrics:  # noqa: ARG002
            return ResidueMetrics(
                plddt=50.0,
                rsa=0.5,
                contact_number=4,
                is_catalytic=False,
                is_disulfide=False,
                is_positive_phi_gly=False,
                is_structural_pro=False,
            )

    def _clean_offtarget(protospacers: List[str]) -> Dict[str, Tuple[float, int, bool]]:
        return {p: (0.8, 0, False) for p in protospacers}

    return GateScorers(
        splice_scorer=SpliceFake(delta=0.0),
        full_structure_scorer=FullStructureFake(pairs={}),
        importance_scorer=_AllPassImportance(),
        dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.83),
        offtarget_scorer=_clean_offtarget,
    )


def _real_bundle(
    *,
    uniprot_acc: str = "",
    cas9: str = "nCas9",  # noqa: ARG001 — documented contract; used by caller
) -> GateScorers:
    from pegrna_pipeline.hpc.complex_dock import Boltz2DockScorer
    from pegrna_pipeline.hpc.full_structure import RiNALMoScorer
    from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer
    from pegrna_pipeline.hpc.structure_importance import AFDBImportanceScorer

    def _real_offtarget(protospacers: List[str]) -> Dict[str, Tuple[float, int, bool]]:
        # Real off-target scoring runs Cas-OFFinder / the offtarget_scorer HPC job
        # against the KOLF2.1J genome. Wired in the HPC follow-on; fail loud here.
        raise RuntimeError(
            "real off-target scoring requires the Cas-OFFinder / offtarget_scorer "
            "HPC job against the KOLF2.1J genome (configure on CX3)"
        )

    return GateScorers(
        splice_scorer=SpliceAIScorer(),
        full_structure_scorer=RiNALMoScorer(),
        importance_scorer=AFDBImportanceScorer(uniprot_acc or "UNKNOWN"),
        dock_scorer=Boltz2DockScorer(),
        offtarget_scorer=_real_offtarget,
    )


def make_scorers(*, mode: str = "fake", **kwargs: Any) -> GateScorers:
    """Return the gate-scorer bundle for the requested ``mode``.

    ``mode="fake"`` for CI/dry-run; ``mode="real"`` for an HPC production run
    (raises on missing deps/weights — never silently downgrades).
    """
    if mode == "fake":
        return _fake_bundle()
    if mode == "real":
        return _real_bundle(**kwargs)
    raise ValueError(f"mode must be 'fake' or 'real'; got {mode!r}")
