"""Scaffold-nCas9 docking gate (decision 7b) — final + heaviest gate.

Boltz-2 (MIT, AF3-class) predicts the full pegRNA + nCas9 complex; we compare the
scaffold↔Cas9 interface confidence (ipTM / interface-pLDDT, + optional
binding-affinity proxy) of the candidate against a cached apo baseline and flag
if docking degrades. ``ComplexDockScorer`` is a Protocol so Boltz-2 (GPU/HPC) or
a fake are interchangeable; the apo baseline is computed once and reused across
candidates. Fallbacks: Chai-1 (contact restraints) > Protenix > RoseTTAFold2NA.

Calibration. This module originally noted that no published pegRNA-Cas9 threshold exists
and that it should be calibrated on real data. It now has been, and the result is that the
AlphaFold-Multimer absolute floor (ipTM >= 0.6) does not transfer: interface ipTM falls with
total RNA length regardless of design quality, so no tjPE-sized pegRNA, neutral or not, can
reach 0.6 (see ``gate_thresholds.DOCKING_REFERENCE``). Scorers that return per-sample ipTMs
are therefore judged by :func:`passes_calibrated_docking_gate`, against references measured
at the candidate's own length. The absolute-floor gate remains for scorers that cannot.
"""

from __future__ import annotations

import statistics
from dataclasses import dataclass
from typing import Optional, Protocol, Sequence, Tuple

from pegrna_pipeline.hpc.gate_thresholds import (
    DOCKING_CALIBRATED_RANGE,
    DOCKING_REFERENCE,
    THRESHOLDS,
)

__all__ = [
    "ComplexDockScorer",
    "DockFake",
    "Boltz2DockScorer",
    "DockingReference",
    "DockingServiceError",
    "docking_reference",
    "passes_calibrated_docking_gate",
    "passes_docking_gate",
]


class DockingServiceError(RuntimeError):
    """The docking SERVICE failed (rate limit, outage, bad credentials, timeout).

    This is not a verdict on a design and must never be recorded as one. Only a
    successfully obtained ipTM may reject a candidate; a rate limit that gets
    reported as "docking scorer error" and skipped turns infrastructure trouble into
    a confident-looking "no valid pegRNA", which is the mirror image of the failure the
    stand-in disclosure guards against (there a placeholder is mistaken for a pass,
    here an outage is mistaken for a fail). The workflow lets this propagate so the
    run fails loudly as what it is.
    """


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


#: AF-Multimer/AF3 boundary between the "high" and "medium" interface-ipTM bands
#: (Evans et al. 2022; Abramson et al. 2024). The medium/likely-incorrect boundary is
#: the credible-interface floor ``THRESHOLDS.docking_min_iptm`` (0.6), so it is not
#: hard-coded here.
_IPTM_HIGH_BAND = 0.8


def _iptm_confidence_band(iptm: float, *, floor: float) -> str:
    """Label the AlphaFold-Multimer / AF3 interface-ipTM confidence band an interface
    falls in (Evans et al. 2022; Abramson et al. 2024 Nature): ``ipTM >= 0.8`` high,
    ``floor <= ipTM < 0.8`` medium, ``ipTM < floor`` likely-incorrect. ``floor`` is the
    credible-interface cutoff (``THRESHOLDS.docking_min_iptm``), so a run using a custom
    floor annotates against that same floor.
    """
    if iptm >= _IPTM_HIGH_BAND:
        return f"high (ipTM>={_IPTM_HIGH_BAND:g})"
    if iptm >= floor:
        return f"medium ({floor:g}-{_IPTM_HIGH_BAND:g})"
    return f"likely-incorrect (<{floor:g})"


#: Appended to the advisory message so a reader sees WHY a long design is not judged, not only that
#: it is not.
_WHY_NOT_JUDGED = (
    " (past about 230 nt the model's score for a normal RNA and a non-binder converge until the "
    "gap is smaller than its own noise)"
)


@dataclass(frozen=True)
class DockingReference:
    """What interface ipTM looks like at one total-RNA length."""

    length: int
    neutral: float  # a pegRNA with a neutral extension of this length
    nonbinder: float  # an RNA of this length that has lost Cas9
    calibrated: bool  # False when ``length`` is outside the measured range

    @property
    def span(self) -> float:
        return self.neutral - self.nonbinder

    def cutoff(self, fraction: float) -> float:
        """The ipTM that sits ``fraction`` of the way from non-binder to neutral."""
        return self.nonbinder + fraction * self.span


def docking_reference(
    length: int, table: Optional[Sequence[Tuple[int, float, float]]] = None
) -> DockingReference:
    """Neutral and non-binder interface ipTM at ``length`` nt (piecewise-linear).

    Interpolates ``DOCKING_REFERENCE``; the endpoints are held flat just beyond the table
    (96 nt is the shortest possible pegRNA; the dynamic range barely changes past 232 nt),
    and ``calibrated`` is False outside ``DOCKING_CALIBRATED_RANGE`` so callers can refuse
    to treat an extrapolated number as a verdict. ``table`` substitutes another
    ``(length, neutral, non-binder)`` table; the calibration analysis uses it to test the
    gate's own interpolation against references built with a sequence left out.
    """
    lo, hi = DOCKING_CALIBRATED_RANGE
    pts = DOCKING_REFERENCE if table is None else table
    x = min(max(length, pts[0][0]), pts[-1][0])
    for (l0, n0, b0), (l1, n1, b1) in zip(pts, pts[1:]):
        if l0 <= x <= l1:
            t = (x - l0) / (l1 - l0)
            return DockingReference(
                length, n0 + t * (n1 - n0), b0 + t * (b1 - b0), lo <= length <= hi
            )
    raise AssertionError("unreachable: DOCKING_REFERENCE must cover its own range")


def passes_calibrated_docking_gate(
    pegrna: str,
    cas9: str,
    *,
    scorer,
    pass_fraction: float = THRESHOLDS.docking_pass_fraction,
    borderline_margin: float = THRESHOLDS.docking_borderline_margin,
    metrics: Optional[dict] = None,
) -> Tuple[bool, str]:
    """Is this pegRNA's Cas9 interface more like a neutral extension of ITS length, or like
    an RNA that has lost Cas9?

    ``scorer.interface_iptm_samples`` returns the interface ipTM of several independent
    diffusion samples; their mean is compared with the references measured at this
    pegRNA's own length (:func:`docking_reference`). The design passes when it lies at least
    ``pass_fraction`` (default 0.5, midway) of the way from the non-binder reference up to
    the neutral one.

    Why not an absolute floor: ipTM falls with RNA length independently of binding quality,
    so 0.6 is unreachable for tjPE-sized pegRNAs (neutral 194/232-nt extensions average
    0.50/0.42). Why not a drop from the apo scaffold: any long RNA "drops". Why keep the gate
    at all: it separates a neutral extension from an RNA that has lost Cas9 (a scrambled
    backbone) with no error in 106 leave-one-out sequences and 12 held-out ones, and it rejects
    an extension built to tie up most of the scaffold (a 40-60 nt complement scores like a
    scrambled backbone). What it does NOT do is detect localized pairing: a 24-nt complement
    to almost any part of the scaffold mostly passes (ADR 0024). It is a check that Cas9 can
    still bind, not a check of the scaffold's folding.

    A verdict within ``borderline_margin`` (as a fraction of the neutral-to-non-binder span)
    of the cutoff is still decided by the rule, but flagged BORDERLINE in the message and in
    ``metrics`` because it is not robust: the closest calls in the calibration data fall in
    that band.

    Outside the calibrated length range there is no valid reference, so the result is
    ADVISORY ONLY: it passes, is flagged in ``metrics['advisory']`` and in the message, and
    the caller surfaces a warning. A safety gate that cannot be validated must say so rather
    than silently accept or reject. The range is where the data support a verdict, not a
    place where measurement stopped: past about 230 nt Boltz-2's score for a normal RNA and for
    a non-binder converge until the gap is inside the model's own sequence-to-sequence noise
    (2.4 sd at 336 nt), and leave-one-out starts to err.
    """
    samples = [float(v) for v in scorer.interface_iptm_samples(pegrna, cas9)]
    if not samples:
        raise DockingServiceError("docking scorer returned no ipTM samples")
    mean = statistics.fmean(samples)
    n = len(samples)
    ref = docking_reference(len(pegrna))
    cutoff = ref.cutoff(pass_fraction)
    frac = (mean - ref.nonbinder) / ref.span
    borderline = ref.calibrated and abs(frac - pass_fraction) < borderline_margin
    borderline_text = (
        f"ipTM {mean:.2f} is {100 * frac:.0f}% of the way from non-binder to neutral "
        f"at {len(pegrna)} nt, within {100 * borderline_margin:.0f} points of the "
        f"{100 * pass_fraction:.0f}% cutoff, so this verdict is not robust"
    )
    borderline_note = "borderline: " + borderline_text
    detail = (
        f"interface ipTM {mean:.2f} (mean of {n}, range {min(samples):.2f}-{max(samples):.2f}) "
        f"at {len(pegrna)} nt; neutral extension {ref.neutral:.2f}, non-binder "
        f"{ref.nonbinder:.2f}; {100 * frac:.0f}% of the way from non-binder to neutral "
        f"(needs >= {100 * pass_fraction:.0f}%, i.e. ipTM >= {cutoff:.2f})"
    )
    if metrics is not None:
        if ref.calibrated:
            metrics.update(
                {
                    "metric": f"interface ipTM (mean of {n})",
                    "value": round(mean, 3),
                    "cutoff": round(cutoff, 3),
                    "cmp": ">=",
                    "extra": f"neutral {ref.neutral:.2f} \u00b7 non-binder {ref.nonbinder:.2f} at {len(pegrna)} nt",
                }
            )
            if borderline:
                metrics.update({"borderline": True, "note": borderline_note})
        else:
            metrics.update(
                {
                    "advisory": True,
                    "note": (
                        f"advisory: ipTM {mean:.2f} at {len(pegrna)} nt is outside the calibrated "
                        f"{DOCKING_CALIBRATED_RANGE[0]}-{DOCKING_CALIBRATED_RANGE[1]} nt range, so it "
                        f"is reported but does not accept or reject{_WHY_NOT_JUDGED}"
                    ),
                }
            )
    if not ref.calibrated:
        return True, (
            f"ADVISORY ONLY - {len(pegrna)} nt is outside the calibrated "
            f"{DOCKING_CALIBRATED_RANGE[0]}-{DOCKING_CALIBRATED_RANGE[1]} nt range, so docking "
            f"cannot be judged here{_WHY_NOT_JUDGED}: {detail}"
        )
    flag = f" [BORDERLINE: {borderline_text}]" if borderline else ""
    if mean < cutoff - 1e-9:
        return (
            False,
            f"scaffold-Cas9 interface looks like a non-binder at this length: {detail}{flag}",
        )
    return True, f"scaffold-Cas9 interface consistent with a neutral extension: {detail}{flag}"


def passes_docking_gate(
    pegrna: str,
    cas9: str,
    *,
    scorer: ComplexDockScorer,
    max_drop: float = THRESHOLDS.docking_max_drop,
    min_iptm: float = THRESHOLDS.docking_min_iptm,
    metrics: Optional[dict] = None,
) -> Tuple[bool, str]:
    """HARD GATE on scaffold-Cas9 docking (interface ipTM).

    Rejects if the candidate interface ipTM either (a) falls below ``min_iptm``
    (AlphaFold-Multimer/AF3 "likely incorrect" floor, 0.6) or (b) drops by more
    than ``max_drop`` (one full confidence band, 0.2) relative to the apo
    scaffold-Cas9 baseline.

    The floor and band annotation are single-sourced from
    :data:`pegrna_pipeline.hpc.gate_thresholds.THRESHOLDS` (``docking_min_iptm`` /
    ``docking_max_drop``), whose citation records the AF-Multimer/AF3 interpretation:
    ipTM ``>= 0.8`` high, ``0.6-0.8`` medium, ``< 0.6`` likely incorrect (Evans et al.
    2022; Abramson et al. 2024 Nature). Every candidate's ipTM is annotated with the band
    it lands in -- so a borderline pass (e.g. 0.64 = medium) and a rejection (e.g. 0.45 =
    likely-incorrect) are both legible in the report, and a near-floor result reads as an
    AF-multimer confidence band rather than an opaque threshold miss.

    If ``metrics`` is supplied it is populated with the raw values and cutoffs
    (interface ipTM vs the min-ipTM floor, plus apo baseline / drop vs max-drop, and the
    AF-multimer band) for the caller's progress display; it never changes the return
    value.

    A scorer that exposes ``interface_iptm_samples`` is judged by
    :func:`passes_calibrated_docking_gate` instead; ``max_drop`` / ``min_iptm`` then do not
    apply (see its docstring for why the absolute floor cannot).
    """
    if hasattr(scorer, "interface_iptm_samples"):
        return passes_calibrated_docking_gate(pegrna, cas9, scorer=scorer, metrics=metrics)
    apo = scorer.apo_interface_iptm(cas9)
    cand = scorer.interface_iptm(pegrna, cas9)
    drop = apo - cand
    band = _iptm_confidence_band(cand, floor=min_iptm)
    if metrics is not None:
        metrics.update(
            {
                "metric": "interface ipTM",
                "value": round(cand, 3),
                "cutoff": min_iptm,
                "cmp": ">=",
                "extra": (
                    f"apo {apo:.2f}->{cand:.2f}, drop {drop:.2f} (max {max_drop}); "
                    f"AF-multimer band: {band}"
                ),
            }
        )
    if cand < min_iptm - 1e-9:
        return False, (
            f"scaffold-Cas9 interface ipTM {cand:.2f} < {min_iptm:g} floor -- AF-multimer "
            f"'{band}' band (Evans 2022; Abramson 2024: >={_IPTM_HIGH_BAND:g} high, "
            f"{min_iptm:g}-{_IPTM_HIGH_BAND:g} medium, <{min_iptm:g} likely incorrect)"
        )
    if drop > max_drop + 1e-9:  # epsilon guards float artefacts at the exact boundary
        return False, (
            f"scaffold-Cas9 interface ipTM dropped {drop:.2f} "
            f"(apo {apo:.2f} -> {cand:.2f} > {max_drop}; {band} band)"
        )
    return True, f"interface preserved (apo {apo:.2f} -> {cand:.2f}; AF-multimer band: {band})"
