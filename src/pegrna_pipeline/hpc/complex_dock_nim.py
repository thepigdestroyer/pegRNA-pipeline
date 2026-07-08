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
"""

from __future__ import annotations

import json
import os
import re
import time
import urllib.request
from typing import Optional

__all__ = ["Boltz2NIMDockScorer"]

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
        timeout: int = 900,
        retries: int = 2,
    ) -> None:
        self.api_key = api_key or os.environ.get("NVIDIA_API_KEY")
        self.url = url or os.environ.get("BOLTZ_NIM_URL") or DEFAULT_NIM_URL
        self.recycling_steps = recycling_steps
        self.sampling_steps = sampling_steps
        self.timeout = timeout
        self.retries = retries
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
    def _predict_iptm(self, rna: str, protein: str) -> float:
        if not self.api_key:
            raise RuntimeError(
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
            "diffusion_samples": 1,
            "step_scale": 1.638,
            "output_format": "mmcif",
            "concatenate_msas": False,
        }
        data = json.dumps(payload).encode()
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}",
        }
        last_err = None
        for attempt in range(self.retries + 1):
            try:
                req = urllib.request.Request(self.url, data=data, headers=headers, method="POST")
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:  # noqa: S310
                    body = json.loads(resp.read().decode())
                return self._extract_iptm(body)
            except Exception as exc:  # noqa: BLE001 - retry transient 5xx/timeouts
                last_err = exc
                if attempt < self.retries:
                    time.sleep(5 * (attempt + 1))
        raise RuntimeError(f"Boltz-2 NIM request failed after {self.retries + 1} tries: {last_err}")

    def _extract_iptm(self, body: dict) -> float:
        """Pull the chain-pair interface ipTM from the NIM response.

        Verified against the live NVIDIA Boltz-2 NIM response schema: the
        interface metrics are TOP-LEVEL lists (one entry per diffusion sample):

        * ``pair_chains_iptm_scores``: ``[{"0": {"0": .., "1": ..}, "1": {..}}]``
          keyed by CHAIN INDEX (order in the polymers list). We send P (rna) then
          C (protein), so the pegRNA<->Cas9 interface is ["0"]["1"] — taken as the
          max off-diagonal to be robust to ordering.
        * ``iptm_scores``: ``[float]`` — the global interface ipTM (fallback).
        * ``confidence_scores``: ``[float]`` — aggregate confidence (last-resort
          proxy; flagged via ``last_metric_source``).
        """
        pcs = body.get("pair_chains_iptm_scores")
        if isinstance(pcs, list) and pcs and isinstance(pcs[0], dict):
            mat = pcs[0]
            off = [
                float(v)
                for i, row in mat.items()
                if isinstance(row, dict)
                for j, v in row.items()
                if i != j
            ]
            if off:
                self.last_metric_source = "pair_chains_iptm_scores"
                return max(off)
        iptm = body.get("iptm_scores")
        if isinstance(iptm, list) and iptm:
            self.last_metric_source = "iptm_scores"
            return float(iptm[0])
        scores = body.get("confidence_scores")
        if isinstance(scores, list) and scores:
            self.last_metric_source = "confidence_scores(proxy)"
            return float(scores[0])
        raise RuntimeError(
            "Boltz-2 NIM response had no usable confidence metric "
            f"(keys: {sorted(body.keys())})"
        )

    # -- ComplexDockScorer interface ----------------------------------------
    def apo_interface_iptm(self, cas9: str) -> float:
        if cas9 in self._apo_cache:
            return self._apo_cache[cas9]
        val = self._predict_iptm(self.apo_scaffold, self._resolve_cas9(cas9))
        self._apo_cache[cas9] = val
        return val

    def interface_iptm(self, pegrna: str, cas9: str) -> float:
        return self._predict_iptm(pegrna, self._resolve_cas9(cas9))
