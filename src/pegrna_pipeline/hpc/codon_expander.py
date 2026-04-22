"""Codon expansion and optimisation.

This module grows over tasks 16-25. Early tasks only rely on ``optimise_bulk``
to codon-optimise raw aa: payloads inside the construct parser; later tasks
add enumerate_variants, metadata, etc.
"""

from __future__ import annotations

# A deliberately compact initial codon map: preferred human codon per amino
# acid. Task 16 replaces this with full ``CODONS_BY_AA`` / frequency-driven
# selection; this placeholder keeps T11 honest without duplicating a motif-
# library codon table inside the agent package.
_SIMPLE_PREFERRED = {
    "A": "GCC", "R": "CGC", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TGA",
}


def optimise_bulk(aa_sequence: str) -> str:
    """Return a DNA string encoding ``aa_sequence`` using preferred codons.

    Later tasks extend this to apply a weighted formula that penalises CpG
    dinucleotides at junctions, restriction sites, and homopolymer runs.
    """
    if not aa_sequence:
        return ""
    out = []
    for aa in aa_sequence:
        upper = aa.upper()
        if upper not in _SIMPLE_PREFERRED:
            raise ValueError(f"unknown amino acid {aa!r}")
        out.append(_SIMPLE_PREFERRED[upper])
    return "".join(out)
