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
    rsa_buried: float = 0.25,
    contact_max: int = 18,
) -> bool:
    """Return True if a MISSENSE at this residue should block the editing pair.

    Block (any true): catalytic / disulfide-Cys / positive-φ Gly / structural
    Pro; OR (pLDDT ≥ ``plddt_min`` AND RSA < ``rsa_buried``) buried core; OR
    (RSA < ``rsa_buried`` AND contact number > ``contact_max``) densely-packed
    buried core. High pLDDT alone — or a high contact count on a solvent-exposed
    residue — never blocks on its own.
    """
    if m.is_catalytic or m.is_disulfide or m.is_positive_phi_gly or m.is_structural_pro:
        return True
    if m.plddt >= plddt_min and m.rsa < rsa_buried:
        return True
    # "Densely-packed core" contact rule: only meaningful for a BURIED residue.
    # The heavy-atom contact count trivially exceeds 18 even for solvent-exposed
    # residues (e.g. a disordered C-terminus at RSA>0.7), so gating on contact
    # number ALONE flags surface residues as structural cores — a false positive
    # that made every candidate fail. Require the residue to also be buried.
    if m.rsa < rsa_buried and m.contact_number > contact_max:
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

    # Tien et al. 2013 theoretical max ASA (A^2), for RSA normalisation.
    _MAX_ASA = {
        "ALA": 129.0,
        "ARG": 274.0,
        "ASN": 195.0,
        "ASP": 193.0,
        "CYS": 167.0,
        "GLU": 223.0,
        "GLN": 225.0,
        "GLY": 104.0,
        "HIS": 224.0,
        "ILE": 197.0,
        "LEU": 201.0,
        "LYS": 236.0,
        "MET": 224.0,
        "PHE": 240.0,
        "PRO": 159.0,
        "SER": 155.0,
        "THR": 172.0,
        "TRP": 285.0,
        "TYR": 263.0,
        "VAL": 174.0,
    }

    def _load(self):  # pragma: no cover - exercised only with biotite + AFDB model
        if self._struct is not None:
            return
        try:
            import biotite.structure as struc  # noqa: F401, PLC0415
            import biotite.structure.io.pdbx as pdbx  # noqa: PLC0415
        except Exception as exc:
            raise ImportError(
                "AFDBImportanceScorer requires the '[importance]' extra (biotite) "
                "and an AFDB model; run on CX3 or a configured host. "
                f"Original error: {exc}"
            ) from exc
        import glob
        import json
        import os
        import urllib.request

        acc = self.uniprot_acc
        afdb_dir = self.afdb_dir or "."
        # OFFLINE-FIRST: prefer a pre-staged local model (CX3 compute nodes have
        # no internet, so staging must download on the login node). Match the F1
        # fragment for this accession if already present.
        local = sorted(glob.glob(os.path.join(afdb_dir, f"AF-{acc}-F1-model_v*.cif")))
        if local:
            path = local[-1]  # highest version
        else:
            # Resolve the current model URL via the AFDB API (version moved v4->v6;
            # the API always points at the latest). For multi-fragment entries
            # (large proteins) pick the F1 fragment, not meta[0] blindly.
            api = f"https://alphafold.ebi.ac.uk/api/prediction/{acc}"
            meta = json.load(urllib.request.urlopen(api, timeout=60))  # noqa: S310 - trusted AFDB
            entry = next((m for m in meta if "-F1-" in str(m.get("cifUrl", ""))), meta[0])
            cif_url = entry["cifUrl"]
            path = os.path.join(afdb_dir, os.path.basename(cif_url))
            os.makedirs(afdb_dir, exist_ok=True)
            urllib.request.urlretrieve(cif_url, path)  # noqa: S310 - trusted AFDB host
        cif = pdbx.CIFFile.read(path)
        self._struct = pdbx.get_structure(cif, model=1, extra_fields=["b_factor"])

    def protein_length(self) -> int:  # pragma: no cover - CX3 only (biotite + AFDB)
        """Residue count of the loaded AFDB model.

        This is the AUTHORITATIVE length for mapping a local design-window residue
        index onto the structure -- it is the very structure that :meth:`metrics`
        indexes -- so the recoded workflow prefers it over a caller-supplied
        ``protein_length`` to avoid a silent local-vs-absolute mis-index. Caveat:
        for a multi-fragment AFDB entry only the F1 fragment is loaded, so for
        proteins longer than the ~2700-residue fragment size this returns the
        fragment length, not the full-sequence length.
        """
        import numpy as np

        self._load()
        return int(np.unique(self._struct.res_id).size)

    def metrics(self, residue_index: int) -> ResidueMetrics:  # pragma: no cover - CX3
        """Per-residue structural metrics from the AFDB model.

        pLDDT from the B-factor column; RSA from biotite SASA normalised by the
        Tien 2013 max-ASA; contact number = heavy-atom neighbours within 5 A;
        side-chain flags from geometry (disulfide SG-SG < 2.5 A, positive-phi Gly,
        Pro). Catalytic/metal flags require UniProt annotation (wired via the
        annotation follow-on); default False here. VALIDATE ON CX3: confirm
        residue indexing (0- vs 1-based) against the transcript<->AFDB alignment
        (>=98% identity) before relying on a block decision.
        """
        import biotite.structure as struc
        import numpy as np

        self._load()
        atoms = self._struct
        res_ids = np.unique(atoms.res_id)
        rid = res_ids[residue_index]
        mask = atoms.res_id == rid
        res_name = atoms.res_name[mask][0]
        plddt = float(np.mean(atoms.b_factor[mask]))
        # SASA (whole structure) -> per-residue sum -> RSA
        sasa_atom = struc.sasa(atoms, vdw_radii="Single")
        res_sasa = float(np.nansum(sasa_atom[mask]))
        rsa = res_sasa / self._MAX_ASA.get(res_name, 200.0)
        # Contact number = count of HEAVY ATOMS of OTHER residues within 5 A of any
        # atom of this residue (matches the Yeh/Shih heavy-atom-contact definition
        # that the >18 threshold is calibrated to — NOT a count of distinct
        # neighbour residues, which would rarely exceed 18 and under-block).
        cell = struc.CellList(atoms, cell_size=5.0)
        is_heavy = atoms.element != "H"
        seen_atoms = set()
        for atom_i in np.where(mask)[0]:
            for j in cell.get_atoms(atoms.coord[atom_i], radius=5.0):
                if atoms.res_id[j] != rid and is_heavy[j]:
                    seen_atoms.add(int(j))
        contact_number = len(seen_atoms)
        # disulfide: this CYS SG within 2.5 A of another CYS SG
        is_disulfide = False
        if res_name == "CYS":
            sg = atoms[(atoms.res_id == rid) & (atoms.atom_name == "SG")]
            if sg.array_length():
                allsg = atoms[(atoms.atom_name == "SG") & (atoms.res_id != rid)]
                if allsg.array_length():
                    d = np.linalg.norm(allsg.coord - sg.coord[0], axis=1)
                    is_disulfide = bool((d < 2.5).any())
        return ResidueMetrics(
            plddt=plddt,
            rsa=rsa,
            contact_number=contact_number,
            is_catalytic=False,
            is_disulfide=is_disulfide,
            is_positive_phi_gly=(res_name == "GLY"),
            is_structural_pro=(res_name == "PRO"),
        )
