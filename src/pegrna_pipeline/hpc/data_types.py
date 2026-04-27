"""Shared frozen dataclasses used across Plan-2 workflows.

This module intentionally has no runtime dependencies beyond the standard
library so it can be imported by any layer of the pipeline (tests, HPC
jobs, the Click CLI) without dragging in Biopython or pandas.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Literal, Optional


@dataclass(frozen=True)
class GeneRecord:
    """Immutable snapshot of a gene's mRNA context.

    Attributes
    ----------
    gene_name:
        The HGNC-style symbol the caller asked for (case-preserved).
    refseq_accession:
        The RefSeq accession that was resolved (e.g. ``NM_006565.4``).
    mrna_sequence:
        Full mRNA sequence as an uppercase ACGT string (no Us).
    cds_start, cds_end:
        Half-open indices of the CDS within ``mrna_sequence``
        (``[cds_start, cds_end)``), 0-indexed.
    stop_codon_start:
        Index of the first nucleotide of the stop codon, i.e.
        ``cds_end - 3``.
    stop_codon:
        3-letter stop codon (``TAA`` | ``TAG`` | ``TGA``).
    upstream_flank:
        ``flanking_bp`` nucleotides immediately preceding the stop codon
        (ends at ``stop_codon_start - 1``).
    downstream_flank:
        ``flanking_bp`` nucleotides immediately after the stop codon
        (starts at ``stop_codon_start + 3``).
    flanking_bp:
        Length of the flanks actually returned.
    """

    gene_name: str
    refseq_accession: str
    mrna_sequence: str
    cds_start: int
    cds_end: int
    stop_codon_start: int
    stop_codon: str
    upstream_flank: str
    downstream_flank: str
    flanking_bp: int


@dataclass(frozen=True)
class PegRNADesign:
    """A single pegRNA row parsed from a PRIDICT 2.0 ``_pegRNA_Pridict_full.csv``.

    ``target_strand`` is either ``"Fw"`` (forward, same strand as mRNA) or
    ``"Rv"`` (reverse complement). ``editing_position`` is the 0-indexed
    position within the input window where the edit lies, matching PRIDICT's
    output convention (see ``tests/test_pridict_index_convention.py``).
    """

    sequence_name: str
    spacer: str
    pbs: str
    rt: str
    pbs_length: int
    rt_length: int
    target_strand: Literal["Fw", "Rv"]
    editing_position: int
    pe2_hek293t: float
    pe2_k562: float
    pegrna_full: str


@dataclass(frozen=True)
class NgRNADesign:
    """A nicking-guide row parsed from a PRIDICT ``_nicking_guides.csv``.

    ``distance_to_edit`` is signed: positive = downstream of the edit on the
    mRNA, negative = upstream. Strand convention mirrors ``PegRNADesign``.
    """

    spacer: str
    distance_to_edit: int
    target_strand: Literal["Fw", "Rv"]
    deepcas9_score: float
    pe3b: bool
    pam_disrupt: bool


#: Default cell line whose PRIDICT 2.0 PE2 score is used to rank pegRNAs.
#: KOLF2 iPSC are closer to K562 than HEK293T in chromatin / PE3 efficiency
#: data, so we default to ``"k562"``. Override per-run via the
#: ``score_cell_line`` field on :class:`VariantScore` or the matching
#: ``score_cell_line`` argument on each workflow's ``run_*_workflow`` entry.
DEFAULT_SCORE_CELL_LINE = "k562"


def _score_attr(cell_line: str) -> str:
    """Translate a cell-line tag to the matching :class:`PegRNADesign` field."""
    cl = cell_line.lower()
    if cl == "k562":
        return "pe2_k562"
    if cl == "hek293t":
        return "pe2_hek293t"
    raise ValueError(f"score_cell_line must be 'k562' or 'hek293t'; got {cell_line!r}")


@dataclass
class VariantScore:
    """Aggregate container: one codon variant, plus all its PRIDICT hits.

    Not frozen — callers routinely build a :class:`VariantScore` incrementally
    (add designs as PRIDICT jobs finish). The ``best_pegRNA`` / ``top_score``
    helpers let sort-by-quality be a one-liner.

    ``score_cell_line`` controls which PRIDICT PE2 column is used for ranking:
    ``"k562"`` (default, recommended for KOLF2 / iPSC work) or ``"hek293t"``.
    """

    sequence_name: str
    codon_variant_dna: str
    pegRNA_designs: List[PegRNADesign] = field(default_factory=list)
    ngRNA_candidates: List[NgRNADesign] = field(default_factory=list)
    score_cell_line: str = DEFAULT_SCORE_CELL_LINE

    @property
    def best_pegRNA(self) -> Optional[PegRNADesign]:
        """Return the highest-scoring design (per ``score_cell_line``); None if empty."""
        if not self.pegRNA_designs:
            return None
        attr = _score_attr(self.score_cell_line)
        return max(self.pegRNA_designs, key=lambda d: getattr(d, attr))

    @property
    def top_score(self) -> float:
        """Convenience: best PE2 score per ``score_cell_line``, or 0.0 if empty.

        Returns a real number (never ``None``) so callers can sort without
        guarding against missing designs.
        """
        best = self.best_pegRNA
        if best is None:
            return 0.0
        return float(getattr(best, _score_attr(self.score_cell_line)))
