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

import os
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
    maxent_scorer: Any = None
    #: PRIDICT2 runner used by the auto-spacer design path (wiring spec 2026-07-10).
    #: ``PridictRunnerFake`` in fake mode; a local-subprocess ``PridictLocalRunner``
    #: in real/api mode when configured, else ``None`` (auto-spacer then fails loud
    #: only when requested — the fixed-spacer path never touches it).
    pridict_runner: Any = None
    #: Gates in this bundle served by a STAND-IN rather than the real model, e.g.
    #: ``("docking",)`` for an api run with no ``NVIDIA_API_KEY``. Threaded into the
    #: run manifest so a report can never pass a stand-in verdict off as a real one
    #: (a design is ordered as DNA off these reports — silence here is a safety bug).
    stub_gates: Tuple[str, ...] = ()


def _local_pridict_runner_or_none() -> Any:
    """Build the local (subprocess) PRIDICT2 runner for the auto-spacer path, or
    ``None`` when it is not configured.

    ``PridictLocalRunner`` fails loud at construction if ``PRIDICT2_DIR`` /
    ``PRIDICT2_PYTHON`` are unset. On a real/api box that only uses the fixed-spacer
    path (no auto-design), that would be a spurious failure, so a missing config
    degrades to ``None`` here; ``design_recoded`` still fails loud with an actionable
    message if auto-spacer is actually requested without a runner. The import is
    light (no torch/TF) and lives here so real-mode heavy-import errors surface first.
    """
    from pegrna_pipeline.hpc.pridict_local import PridictLocalRunner

    try:
        return PridictLocalRunner()
    except ValueError:
        return None


def _fake_bundle() -> GateScorers:
    from pegrna_pipeline.hpc.complex_dock import DockFake
    from pegrna_pipeline.hpc.cryptic_splice import MaxEntScanFake
    from pegrna_pipeline.hpc.full_structure import FullStructureFake
    from pegrna_pipeline.hpc.pridict_local import PridictRunnerFake
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
        maxent_scorer=MaxEntScanFake(),  # all motifs below threshold -> gate passes
        pridict_runner=PridictRunnerFake(),  # deterministic auto-spacer source for CI
        # Every gate here is a stand-in by construction — fake mode is a dry run.
        stub_gates=(
            "splice",
            "structure",
            "importance",
            "docking",
            "offtarget",
            "cryptic",
            "pridict",
        ),
    )


def _real_bundle(
    *,
    uniprot_acc: str = "",
    cas9: str = "nCas9",  # noqa: ARG001 — documented contract; used by caller
    genome: str | None = None,
    pam: str = "NGG",
    max_mismatches: int = 3,
    boltz_cache: str | None = None,
    rinalmo_weights: str | None = None,
    rnaformer_model_dir: str | None = None,
    **_ignored: Any,
) -> GateScorers:
    from pegrna_pipeline.hpc.complex_dock import Boltz2DockScorer
    from pegrna_pipeline.hpc.cryptic_splice import MaxEntScanScorer
    from pegrna_pipeline.hpc.full_structure import AutoStructureScorer
    from pegrna_pipeline.hpc.offtarget_scorer import CasOffinderScorer
    from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer
    from pegrna_pipeline.hpc.structure_importance import AFDBImportanceScorer

    if genome:
        # Real genome-wide specificity via cas-offinder (+ KOLF2.1J FASTA on CX3).
        offtarget: Callable[[List[str]], Dict[str, Tuple[float, int, bool]]] = CasOffinderScorer(
            genome=genome,
            pam=pam,
            max_mismatches=max_mismatches,
            # Per-guide results are cached under PEGRNA_CACHE_DIR when it is set: a scan costs
            # ~69 CPU-minutes, and the same guides recur across runs.
            cache_dir=os.environ.get("PEGRNA_CACHE_DIR"),
        )
    else:

        def offtarget(protospacers: List[str]) -> Dict[str, Tuple[float, int, bool]]:
            # Fail loud (never silently pass): real off-target needs a genome.
            raise RuntimeError(
                "real off-target scoring requires genome= (cas-offinder + the "
                "KOLF2.1J FASTA on CX3); pass genome=/path/to/genome.fa to make_scorers"
            )

    return GateScorers(
        splice_scorer=SpliceAIScorer(),
        # RiNALMo (Ampere+ GPU) with a transparent RNAformer (CPU) fallback so the
        # structure gate runs even without FlashAttention — both real models.
        full_structure_scorer=AutoStructureScorer(
            rinalmo_weights=rinalmo_weights, rnaformer_model_dir=rnaformer_model_dir
        ),
        importance_scorer=AFDBImportanceScorer(uniprot_acc or "UNKNOWN"),
        # --cache points Boltz at pre-staged weights (offline CX3 compute nodes).
        dock_scorer=Boltz2DockScorer(cache_dir=boltz_cache),
        offtarget_scorer=offtarget,
        # MaxEntScan cryptic-splice gate (maxentpy staged in the [splice] extra).
        maxent_scorer=MaxEntScanScorer(),
        # Local-subprocess PRIDICT2 for the auto-spacer path (None if unconfigured).
        pridict_runner=_local_pridict_runner_or_none(),
    )


class _PaddedSpliceAIScorer:
    """SpliceAI on short local design windows (api/CPU path).

    The real ``SpliceAIScorer`` fail-louds unless it is handed >= 2*context nt
    (>= 5 kb each side) — the HPC path supplies the true GRCh38 flanks. On a
    GPU-less web host we usually only have the local design window (hundreds of
    nt), so we N-pad both REF and ALT up to the required context before scoring.
    N is SpliceAI's own end-padding symbol (zero one-hot, no splice signal), and
    the prefix/suffix alignment in ``delta`` keeps the edit region intact, so the
    delta over the real bases is preserved. CAVEAT: without the true genomic
    neighbourhood this is an approximation of the full-context HPC gate.
    """

    def __init__(self) -> None:
        from pegrna_pipeline.hpc.splice_screen import SpliceAIScorer

        self._inner = SpliceAIScorer()

    def delta(self, ref_seq: str, alt_seq: str) -> float:
        pad = self._inner._context  # 5000 nt each side
        need = 2 * pad
        if len(ref_seq) >= need and len(alt_seq) >= need:
            return self._inner.delta(ref_seq, alt_seq)
        flank = "N" * pad
        return self._inner.delta(flank + ref_seq + flank, flank + alt_seq + flank)


def _api_bundle(
    *,
    uniprot_acc: str = "",
    cas9: str = "nCas9",  # noqa: ARG001 — documented contract; used by caller
    genome: str | None = None,
    pam: str = "NGG",
    max_mismatches: int = 3,
    nvidia_api_key: str | None = None,
    **_ignored: Any,
) -> GateScorers:
    """CPU/API gate bundle for a GPU-less host (the web deployment).

    Every gate is REAL except that the two GPU-hard models are handled without a
    local GPU:

    * splice      -> _PaddedSpliceAIScorer      (SpliceAIScorer on CPU, N-padded so
                                                 the advisory cross-check runs on
                                                 short local design windows)
    * cryptic     -> MaxEntScanScorer          (pure-Python matrices, CPU)
    * structure   -> ViennaRNAStructureScorer  (partition-function bpp, CPU)
    * importance  -> AFDBImportanceScorer       (AlphaFold-DB web API + biotite, CPU)
    * docking     -> Boltz-2 via NVIDIA-hosted NIM API (needs NVIDIA_API_KEY);
                     if no key is configured, a transparent stand-in is used so a
                     design still runs, and this is flagged in the report manifest.
    * off-target  -> cas-offinder if a genome FASTA is supplied, else clean-pass.
    """
    from pegrna_pipeline.hpc.complex_dock import DockFake
    from pegrna_pipeline.hpc.complex_dock_nim import Boltz2NIMDockScorer
    from pegrna_pipeline.hpc.cryptic_splice import MaxEntScanScorer
    from pegrna_pipeline.hpc.structure_importance import AFDBImportanceScorer
    from pegrna_pipeline.hpc.vienna_structure import ViennaRNAStructureScorer

    stubs: List[str] = []
    key = nvidia_api_key or os.environ.get("NVIDIA_API_KEY")
    if key:
        dock: Any = Boltz2NIMDockScorer(api_key=key)
    else:
        # No key -> cannot call the hosted GPU model. Run every other gate for real
        # and use a neutral stand-in for docking rather than failing the whole
        # design. Recorded in ``stub_gates`` so the report says so in plain text.
        dock = DockFake(apo_iptm=0.85, candidate_iptm=0.83)
        stubs.append("docking")

    if genome:
        from pegrna_pipeline.hpc.offtarget_scorer import CasOffinderScorer

        offtarget: Callable[[List[str]], Dict[str, Tuple[float, int, bool]]] = CasOffinderScorer(
            genome=genome,
            pam=pam,
            max_mismatches=max_mismatches,
            # Per-guide results are cached under PEGRNA_CACHE_DIR when it is set: a scan costs
            # ~69 CPU-minutes, and the same guides recur across runs.
            cache_dir=os.environ.get("PEGRNA_CACHE_DIR"),
        )
    else:
        # No genome FASTA -> cas-offinder cannot run. Clean-pass stand-in, declared.
        stubs.append("offtarget")

        def offtarget(protospacers: List[str]) -> Dict[str, Tuple[float, int, bool]]:
            return {p: (0.8, 0, False) for p in protospacers}

    return GateScorers(
        # N-pad short web-host design windows up to SpliceAI's >=10 kb context so the
        # advisory cross-check actually runs here; the raw SpliceAIScorer fail-louds
        # below that and the workflow would silently record it as an "advisory skip".
        splice_scorer=_PaddedSpliceAIScorer(),
        full_structure_scorer=ViennaRNAStructureScorer(),
        importance_scorer=AFDBImportanceScorer(uniprot_acc or "UNKNOWN"),
        dock_scorer=dock,
        offtarget_scorer=offtarget,
        maxent_scorer=MaxEntScanScorer(),
        # Local-subprocess PRIDICT2 for the auto-spacer path (None if unconfigured).
        pridict_runner=_local_pridict_runner_or_none(),
        stub_gates=tuple(stubs),
    )


def make_scorers(*, mode: str = "fake", **kwargs: Any) -> GateScorers:
    """Return the gate-scorer bundle for the requested ``mode``.

    ``mode="fake"`` for CI/dry-run; ``mode="real"`` for an HPC production run
    (raises on missing deps/weights — never silently downgrades); ``mode="api"``
    for a GPU-less host — CPU-real gates + Boltz-2 via the NVIDIA-hosted NIM API.
    """
    if mode == "fake":
        return _fake_bundle()
    if mode == "real":
        return _real_bundle(**kwargs)
    if mode == "api":
        return _api_bundle(**kwargs)
    raise ValueError(f"mode must be 'fake', 'real', or 'api'; got {mode!r}")
