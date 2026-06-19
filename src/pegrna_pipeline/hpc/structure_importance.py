"""AlphaFold structural-importance hard gate (per-residue).

Conservative OR across hard-block axes (when in doubt, block). Metrics come from
an AlphaFold model: pLDDT/PAE from AlphaFold-DB, RSA via DSSP, contacts via
biotite, side-chain flags via geometry + UniProt/M-CSA annotation. The
``ImportanceScorer`` Protocol lets that heavy extraction be injected/faked; the
pure decision rules (:func:`residue_blocked`, :func:`indel_blocked`) are fully
testable without any model.

Fail-loud (handled by the caller): no AFDB entry, footprint identity < 98 %, or
ligand/partner present (monomer insufficient → route to AF-multimer).

ROOT RULE: AlphaFold predicts ONE static monomer + a CONFIDENCE score; it does
NOT predict ΔΔG. pLDDT is a gate-ENABLER (must pair with burial), never a
standalone block. Biggest false-PASS = intermolecular interface residues absent
from the monomer — mitigate by OR-ing UniProt BINDING/METAL/DNA_BIND.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Protocol

__all__ = [
    "ResidueMetrics",
    "ImportanceScorer",
    "residue_blocked",
    "indel_blocked",
    "span_blocked",
    "AFDBImportanceScorer",
]

#: DSSP secondary-structure codes where an in-frame indel disrupts the backbone.
_SSE_BLOCK_CODES = frozenset({"H", "G", "I", "E", "B"})


@dataclass(frozen=True)
class ResidueMetrics:
    plddt: float
    rsa: float
    contact_number: int
    is_catalytic: bool
    is_disulfide: bool
    is_positive_phi_gly: bool
    is_structural_pro: bool


def residue_blocked(
    m: ResidueMetrics,
    *,
    plddt_min: float = 70.0,
    rsa_buried: float = 0.15,
    contact_max: int = 18,
) -> bool:
    """Return True if a MISSENSE at this residue should block the editing pair.

    Block (any true): catalytic / disulfide-Cys / positive-φ Gly / structural
    Pro; OR (pLDDT ≥ ``plddt_min`` AND RSA < ``rsa_buried``) buried core; OR
    contact number > ``contact_max``. High pLDDT alone never blocks.
    """
    if m.is_catalytic or m.is_disulfide or m.is_positive_phi_gly or m.is_structural_pro:
        return True
    if m.plddt >= plddt_min and m.rsa < rsa_buried:
        return True
    if m.contact_number > contact_max:
        return True
    return False


def indel_blocked(dssp_code: str) -> bool:
    """Return True if an in-frame indel landing on this DSSP code should block.

    Indels disrupt the backbone, so any helix/strand element (H/G/I/E/B) blocks;
    loops/turns/coil (S/T/C/' ') are tolerated.
    """
    return dssp_code.strip().upper()[:1] in _SSE_BLOCK_CODES


class ImportanceScorer(Protocol):
    def metrics(self, residue_index: int) -> ResidueMetrics: ...


def span_blocked(scorer: ImportanceScorer, residue_indices: Iterable[int]) -> bool:
    """Block if ANY residue in the footprint is structurally important (OR)."""
    return any(residue_blocked(scorer.metrics(i)) for i in residue_indices)


class AFDBImportanceScorer:
    """Real AlphaFold-DB-backed scorer (lazy heavy imports; AFDB + DSSP + biotite).

    Fetches the AFDB model for a UniProt accession, computes RSA (DSSP), contact
    number (biotite), and side-chain flags (geometry + UniProt/M-CSA), then maps
    transcript→residue through an alignment (≥98 % identity or fail-loud). Heavy
    deps load lazily so the core package + CI never import them. Real runs happen
    where biotite/DSSP + the AFDB model are available (CX3 or a configured box).
    """

    def __init__(self, uniprot_acc: str, *, afdb_dir: str | None = None) -> None:
        self.uniprot_acc = uniprot_acc
        self.afdb_dir = afdb_dir
        self._struct = None

    def _load(self):  # pragma: no cover - exercised only with biotite/DSSP + model
        if self._struct is not None:
            return
        try:
            import biotite.structure.io.pdbx as pdbx  # noqa: F401, PLC0415
            from biotite.database import rcsb  # noqa: F401, PLC0415
        except Exception as exc:
            raise ImportError(
                "AFDBImportanceScorer requires the '[importance]' extra "
                "(biotite + a DSSP binary) and an AFDB model; run on CX3 or a "
                f"configured host. Original error: {exc}"
            ) from exc
        raise NotImplementedError(
            "AFDB fetch + DSSP/biotite metric extraction is wired in the follow-on plan"
        )

    def metrics(self, residue_index: int) -> ResidueMetrics:  # pragma: no cover - CX3
        self._load()
        raise NotImplementedError
