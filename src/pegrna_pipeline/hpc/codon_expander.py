"""Codon expansion and optimisation.

Core building blocks used across the pipeline:
- :data:`CODONS_BY_AA`  -- canonical codon map, cytosine-first per aa.
- :data:`HUMAN_CODON_FREQUENCIES` -- flat dict of {codon: fraction} from the
  ``h_sapiens_9606`` table shipped by ``python-codon-tables``.
- :func:`codons_for_aa` / :func:`codon_frequency` -- convenience accessors.
- :func:`optimise_bulk` -- pick one codon per residue using a weighted score.
- (added in later tasks) :func:`enumerate_variants`,
  :func:`enumerate_with_metadata`, quality filters, etc.
"""

from __future__ import annotations

from typing import Dict, List

import python_codon_tables as _pct


# ---------------------------------------------------------------------- tables


def _load_human_codon_table() -> Dict[str, Dict[str, float]]:
    return _pct.get_codons_table("h_sapiens_9606")


def _build_codons_by_aa() -> Dict[str, List[str]]:
    """Canonical codon lists per amino acid, sorted cytosine-first then lex.

    Cytosine-first ordering matters for downstream enumeration where the
    pipeline prefers C-rich codons (GC-content stability in the C-terminal
    edit window).
    """
    table = _load_human_codon_table()
    out: Dict[str, List[str]] = {}
    for aa, codon_map in table.items():
        codons = list(codon_map.keys())
        # Stable, cytosine-first, then lexicographic.
        codons.sort(key=lambda c: (0 if c.startswith("C") else 1, c))
        out[aa] = codons
    return out


CODONS_BY_AA: Dict[str, List[str]] = _build_codons_by_aa()


def _build_flat_frequency_table() -> Dict[str, float]:
    table = _load_human_codon_table()
    flat: Dict[str, float] = {}
    for aa, codon_map in table.items():
        for codon, freq in codon_map.items():
            flat[codon] = float(freq)
    return flat


HUMAN_CODON_FREQUENCIES: Dict[str, float] = _build_flat_frequency_table()


# --------------------------------------------------------------------- lookups


def codons_for_aa(aa: str) -> List[str]:
    """Return the cytosine-first list of codons encoding ``aa``.

    Case-insensitive. Raises :class:`KeyError` for unknown aa letters.
    """
    key = aa.upper()
    if key not in CODONS_BY_AA:
        raise KeyError(aa)
    return list(CODONS_BY_AA[key])


def codon_frequency(aa: str, codon: str) -> float:
    """Return the human frequency of ``codon`` within the synonymous group for ``aa``."""
    key = aa.upper()
    table = _load_human_codon_table()
    try:
        return float(table[key][codon.upper()])
    except KeyError as exc:
        raise KeyError(f"codon {codon!r} not listed for aa {aa!r}") from exc


# ---------------------------------------------------------- bulk optimisation


def optimise_bulk(aa_sequence: str) -> str:
    """Preliminary greedy codon optimiser.

    Picks the most-frequent human codon per residue; ties are broken by the
    cytosine-first ordering of :data:`CODONS_BY_AA`. Task 25 upgrades this to
    a weighted-score policy that penalises CpG junctions, restriction sites,
    and homopolymer runs.

    An empty string returns ``""``. Unknown amino acids raise
    :class:`ValueError`.
    """
    if not aa_sequence:
        return ""
    out: List[str] = []
    for aa in aa_sequence:
        key = aa.upper()
        if key not in CODONS_BY_AA:
            raise ValueError(f"unknown amino acid {aa!r}")
        codons = CODONS_BY_AA[key]
        # Pick highest-frequency codon; stable tie-break favours cytosine-first.
        best = max(codons, key=lambda c: (codon_frequency(key, c), codons.index(c) * -1))
        out.append(best)
    return "".join(out)
