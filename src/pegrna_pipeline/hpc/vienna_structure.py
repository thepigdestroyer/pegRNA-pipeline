"""CPU base-pair-probability scorer via ViennaRNA (for the `api`/CPU gate bundle).

RiNALMo needs an Ampere+ GPU and RNAformer needs a staged checkpoint; neither is
available on a GPU-less host. ViennaRNA's partition-function base-pair
probabilities give a real, fast, CPU-only ``pair_probs`` implementing the same
``FullStructureScorer`` Protocol (``full_structure.py``), so the 2D-structure
gate runs for real without a GPU. This is a physics-based (Turner nearest-
neighbour) model rather than a learned one — a legitimate cross-check the repo
already names as the "cheap inner-loop pre-screen", promoted here to the gate
scorer when learned models can't run.

Indexing matches the RiNALMo/RNAformer adapters: 0-based sequence positions,
``{(i, j): p}`` for i < j, full map above a small floor (so region_accessibility
sees the total paired mass).
"""

from __future__ import annotations

from typing import Dict, Tuple

__all__ = ["ViennaRNAStructureScorer"]


class ViennaRNAStructureScorer:
    """CPU base-pair-probability scorer using ViennaRNA's partition function.

    Requires the ``ViennaRNA`` package (``pip install ViennaRNA`` or
    ``conda install -c bioconda viennarna``). Loads lazily and raises a helpful
    ImportError if missing (fail-loud, consistent with the other real adapters).
    """

    def __init__(self, *, floor: float = 0.01) -> None:
        self.floor = floor

    def pair_probs(self, seq: str) -> Dict[Tuple[int, int], float]:
        try:
            import RNA  # noqa: PLC0415
        except Exception as exc:  # pragma: no cover - env-dependent
            raise ImportError(
                "ViennaRNAStructureScorer requires the ViennaRNA python package "
                "(pip install ViennaRNA / conda install -c bioconda viennarna). "
                f"Original error: {exc}"
            ) from exc

        # Fold the RNA form: uppercase, DNA T -> U (the emitted pegRNA mixes DNA
        # letters with a PolyU tail; the Turner model is RNA).
        rna = seq.upper().replace("T", "U")
        fc = RNA.fold_compound(rna)
        fc.pf()  # compute the partition function (populates the bpp matrix)
        bpp = fc.bpp()  # 1-based: bpp[i][j] = P(pair i,j) for 1 <= i < j <= n

        n = len(rna)
        out: Dict[Tuple[int, int], float] = {}
        for i in range(1, n + 1):
            row = bpp[i]
            for j in range(i + 1, n + 1):
                p = float(row[j])
                if p >= self.floor:
                    out[(i - 1, j - 1)] = p  # -> 0-based to match the other scorers
        return out
