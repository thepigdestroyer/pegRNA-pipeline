"""Boltz-2 docking gate via the NVIDIA-hosted NIM API (no local GPU needed).

Drop-in ``ComplexDockScorer`` (see ``complex_dock.py``) that offloads the heavy
Boltz-2 pegRNA+Cas9 complex prediction to NVIDIA's hosted Boltz-2 NIM
(https://health.api.nvidia.com/v1/biology/mit/boltz2/predict) instead of a local
GPU. Implements ``interface_iptm`` / ``apo_interface_iptm`` with the same
semantics as ``Boltz2DockScorer`` so ``passes_docking_gate`` is unchanged.

Two things the upstream local adapter left as "validate on CX3" and this one
resolves:

* ``cas9`` arrives as the NAME ``"nCas9"`` (not a sequence). We resolve it to the
  real S. pyogenes Cas9 sequence (UniProt Q99ZW2), applying the D10A nickase
  mutation for ``nCas9``. Fetched once from UniProt and cached; override with the
  ``PEGRNA_CAS9_SEQ`` env var or by passing a raw sequence as ``cas9``.
* ipTM extraction: the NIM response reliably carries ``confidence_scores`` and
  the predicted mmCIF; interface ipTM (``pair_chains_iptm`` / ``iptm``) may live
  under ``metrics`` depending on NIM version. We look for a real ipTM first and
  fall back to the aggregate ``confidence_scores`` value as a proxy (the gate is
  a RELATIVE candidate-vs-apo comparison, so a consistent proxy still functions);
  the proxy path is flagged via ``self.last_metric_source``.

Cost/latency: a Cas9(~1368 aa)+pegRNA(~700 nt) complex is large — each call can
take minutes and consume NIM credits. The apo baseline is cached and reused.

Waiting model (measured against the live endpoint, not assumed). A plain request is held
open until the result is ready: 25-150 s for the sizes used here. Do NOT send NVIDIA's
``NVCF-POLL-SECONDS`` header, although the cloud-function docs describe it as returning a
202 plus a request id to poll. On this hosted endpoint a short window returns HTTP 504
(``Nvcf-Status: errored``) when it elapses, and ``GET /v2/nvcf/pexec/status/<id>`` for that
id is 404, so the job is gone, not resumable. A 504 therefore means "resubmit", which is
how it is handled. Because a shared trial endpoint also stalls (a request once hung for
over 30 minutes), every call has a per-attempt socket timeout and an overall deadline.
"""

from __future__ import annotations

import json
import os
import re
import time
import urllib.error
import urllib.request
from typing import Optional

from pegrna_pipeline.hpc.complex_dock import DockingServiceError
from pegrna_pipeline.hpc.gate_thresholds import THRESHOLDS

__all__ = ["Boltz2NIMDockScorer"]

#: HTTP statuses worth retrying: rate limiting and transient server trouble.
_TRANSIENT_HTTP = frozenset({429, 500, 502, 503, 504})
_MAX_BACKOFF_S = 60.0
DEFAULT_ATTEMPT_TIMEOUT_S = 300.0
DEFAULT_DEADLINE_S = 600.0

DEFAULT_NIM_URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
_SPCAS9_UNIPROT = "Q99ZW2"  # S. pyogenes Cas9


class Boltz2NIMDockScorer:
    def __init__(
        self,
        *,
        api_key: Optional[str] = None,
        url: Optional[str] = None,
        apo_scaffold: Optional[str] = None,
        recycling_steps: int = 3,
        sampling_steps: int = 50,
        timeout: float = DEFAULT_ATTEMPT_TIMEOUT_S,
        retries: int = 6,
        samples: int = THRESHOLDS.docking_samples,
        deadline: float = DEFAULT_DEADLINE_S,
    ) -> None:
        self.api_key = api_key or os.environ.get("NVIDIA_API_KEY")
        self.url = url or os.environ.get("BOLTZ_NIM_URL") or DEFAULT_NIM_URL
        self.recycling_steps = recycling_steps
        self.sampling_steps = sampling_steps
        #: Socket timeout of ONE attempt. Real calls take 25-150 s; 300 s detects a hung
        #: connection in minutes instead of the quarter of an hour it used to take.
        self.timeout = timeout
        self.retries = retries
        #: Wall-clock budget for one call INCLUDING retries and backoff. Without it
        #: ``retries`` x ``timeout`` could hold a web request for the better part of an hour.
        self.deadline = deadline
        #: Diffusion samples per request. A single draw is too noisy to threshold (sd
        #: 0.03-0.11; a scrambled non-binder drew 0.634 once), so the calibrated gate
        #: averages several, which cost one request and about the same wall time.
        self.samples = max(1, int(samples))
        if apo_scaffold is None:
            from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD  # noqa: PLC0415

            apo_scaffold = STANDARD_SGRNA_SCAFFOLD
        self.apo_scaffold = apo_scaffold
        self._apo_cache: dict = {}
        self._cas9_seq_cache: dict = {}
        self.last_metric_source: Optional[str] = None

    # -- Cas9 sequence resolution -------------------------------------------
    def _resolve_cas9(self, cas9: str) -> str:
        """Map a Cas9 NAME to a protein sequence; pass a raw AA sequence through."""
        env_seq = os.environ.get("PEGRNA_CAS9_SEQ")
        if env_seq:
            return env_seq.strip()
        # A raw amino-acid sequence (long, only AA letters) is used as-is.
        if len(cas9) > 40 and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", cas9.upper()):
            return cas9.upper()
        name = cas9.lower()
        if name in self._cas9_seq_cache:
            return self._cas9_seq_cache[name]
        seq = self._fetch_uniprot(_SPCAS9_UNIPROT)
        if name in ("ncas9", "nickase", "cas9n"):
            # D10A nickase (RuvC-dead): position 10 (1-based) Asp -> Ala.
            if len(seq) >= 10 and seq[9] == "D":
                seq = seq[:9] + "A" + seq[10:]
        self._cas9_seq_cache[name] = seq
        return seq

    @staticmethod
    def _fetch_uniprot(acc: str) -> str:
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        with urllib.request.urlopen(url, timeout=60) as r:  # noqa: S310
            text = r.read().decode()
        return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))

    # -- NIM call ------------------------------------------------------------
    def _predict_samples(self, rna: str, protein: str) -> list:
        if not self.api_key:
            raise DockingServiceError(
                "Boltz2NIMDockScorer needs an NVIDIA API key: set NVIDIA_API_KEY "
                "(nvapi-...) from build.nvidia.com. No key -> cannot run the docking gate."
            )
        # Boltz RNA chains must be ACGU; the pegRNA/scaffold are stored in DNA
        # letters, so transcribe T->U (a raw T here 422s the NIM).
        rna_acgu = rna.upper().replace("T", "U")
        payload = {
            "polymers": [
                {"id": "P", "molecule_type": "rna", "sequence": rna_acgu},
                {"id": "C", "molecule_type": "protein", "sequence": protein.upper()},
            ],
            "recycling_steps": self.recycling_steps,
            "sampling_steps": self.sampling_steps,
            "diffusion_samples": self.samples,
            "step_scale": 1.638,
            "output_format": "mmcif",
            "concatenate_msas": False,
        }
        data = json.dumps(payload).encode()
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}",
        }
        last_err: Optional[BaseException] = None
        started = time.monotonic()
        attempts = 0
        for attempt in range(self.retries + 1):
            remaining = self.deadline - (time.monotonic() - started)
            if remaining <= 0:
                break
            attempts += 1
            try:
                req = urllib.request.Request(self.url, data=data, headers=headers, method="POST")
                # An attempt may not outlive the call's overall budget.
                with urllib.request.urlopen(  # noqa: S310
                    req, timeout=min(self.timeout, remaining)
                ) as resp:
                    body = json.loads(resp.read().decode())
                return self._extract_samples(body)
            except urllib.error.HTTPError as exc:
                last_err = exc
                if exc.code not in _TRANSIENT_HTTP:
                    # 401/403 (key), 422 (payload), ...: retrying cannot help.
                    detail = ""
                    try:
                        detail = exc.read().decode(errors="replace")[:300]
                    except Exception:  # noqa: BLE001
                        pass
                    raise DockingServiceError(
                        f"Boltz-2 NIM rejected the request (HTTP {exc.code}): {detail}"
                    ) from exc
                retry_after = exc.headers.get("Retry-After") if exc.headers else None
                try:
                    wait = float(retry_after) if retry_after else 0.0
                except ValueError:
                    wait = 0.0
                wait = wait or min(5.0 * 2**attempt, _MAX_BACKOFF_S)
            except (urllib.error.URLError, TimeoutError, OSError, ValueError) as exc:
                # ValueError covers a non-JSON body (json.JSONDecodeError): a flaky proxy
                # or a truncated response, retried like a dropped connection.
                last_err = exc
                wait = min(5.0 * 2**attempt, _MAX_BACKOFF_S)
            if attempt < self.retries:
                pause = min(wait, _MAX_BACKOFF_S)
                if pause >= self.deadline - (time.monotonic() - started):
                    break  # sleeping would use up what is left; fail now rather than later
                time.sleep(pause)
        elapsed = time.monotonic() - started
        raise DockingServiceError(
            f"Boltz-2 NIM unavailable after {attempts} attempt(s) in {elapsed:.0f} s "
            f"(budget {self.deadline:.0f} s): {last_err or 'no attempt was possible'}. "
            "This is the docking service failing, not a verdict on the design."
        )

    @staticmethod
    def _interface_value(matrix: dict) -> Optional[float]:
        """The pegRNA<->Cas9 entry of one sample's pair-chains matrix (max off-diagonal)."""
        off = [
            float(v)
            for i, row in matrix.items()
            if isinstance(row, dict)
            for j, v in row.items()
            if i != j
        ]
        return max(off) if off else None

    def _extract_samples(self, body: dict) -> list:
        """The interface ipTM of EVERY diffusion sample in the NIM response.

        Verified against the live NVIDIA Boltz-2 NIM response schema: the interface
        metrics are TOP-LEVEL lists with one entry per diffusion sample:

        * ``pair_chains_iptm_scores``: ``[{"0": {"0": .., "1": ..}, "1": {..}}, ...]``
          keyed by CHAIN INDEX (order in the polymers list). We send P (rna) then C
          (protein), so the pegRNA<->Cas9 interface is the off-diagonal entry.
        * ``iptm_scores``: ``[float, ...]`` -- the global interface ipTM (fallback).
        * ``confidence_scores``: ``[float, ...]`` -- aggregate confidence (last-resort
          proxy; flagged via ``last_metric_source``).
        """
        pcs = body.get("pair_chains_iptm_scores")
        if isinstance(pcs, list) and pcs and all(isinstance(m, dict) for m in pcs):
            vals = [self._interface_value(m) for m in pcs]
            if vals and all(v is not None for v in vals):
                self.last_metric_source = "pair_chains_iptm_scores"
                return vals
        for key, source in (
            ("iptm_scores", "iptm_scores"),
            ("confidence_scores", "confidence_scores(proxy)"),
        ):
            scores = body.get(key)
            if isinstance(scores, list) and scores:
                self.last_metric_source = source
                return [float(v) for v in scores]
        raise DockingServiceError(
            "Boltz-2 NIM response had no usable confidence metric " f"(keys: {sorted(body.keys())})"
        )

    def _extract_iptm(self, body: dict) -> float:
        """First sample's interface ipTM (kept for callers that want a single number)."""
        return self._extract_samples(body)[0]

    def _predict_iptm(self, rna: str, protein: str) -> float:
        """Mean interface ipTM over the diffusion samples (legacy single-number view)."""
        vals = self._predict_samples(rna, protein)
        return sum(vals) / len(vals)

    # -- ComplexDockScorer interface ----------------------------------------
    def apo_interface_iptm(self, cas9: str) -> float:
        if cas9 in self._apo_cache:
            return self._apo_cache[cas9]
        val = self._predict_iptm(self.apo_scaffold, self._resolve_cas9(cas9))
        self._apo_cache[cas9] = val
        return val

    def interface_iptm(self, pegrna: str, cas9: str) -> float:
        return self._predict_iptm(pegrna, self._resolve_cas9(cas9))

    def interface_iptm_samples(self, pegrna: str, cas9: str) -> list:
        """Per-sample interface ipTM, for the length-calibrated gate.

        The calibration in ``gate_thresholds.DOCKING_REFERENCE`` was measured on
        ``pair_chains_iptm_scores``. If the service stops returning it and the scorer falls
        back to the global ipTM or an aggregate-confidence proxy, comparing that against
        references for a different quantity would be a silent miscalibration, so refuse.
        """
        vals = self._predict_samples(pegrna, self._resolve_cas9(cas9))
        if self.last_metric_source != "pair_chains_iptm_scores":
            raise DockingServiceError(
                f"Boltz-2 NIM returned '{self.last_metric_source}' instead of "
                "pair_chains_iptm_scores; the docking calibration was measured on the latter "
                "and does not apply to this metric."
            )
        return vals
