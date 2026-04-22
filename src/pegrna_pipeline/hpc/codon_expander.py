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

import heapq
import itertools
import re
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from typing import Dict, Iterator, List, Tuple

import python_codon_tables as _pct
import yaml


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


# ---------------------------------------------------------------- CodonVariant


@dataclass(frozen=True)
class CodonVariant:
    """A candidate DNA variant for a given aa sequence.

    Attributes
    ----------
    dna:
        Uppercase ACGT string encoding the aa sequence.
    c_count:
        Number of cytosine residues in the DNA (higher -> preferred).
    gc_content:
        GC fraction in [0, 1].
    rank:
        Position in a sorted variant list (0 = best).
    """

    dna: str
    c_count: int
    gc_content: float
    rank: int


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


# --------------------------------------------------------- variant enumeration


def _surviving_variants(aa_sequence: str) -> Iterator[Tuple[str, int, float]]:
    """Yield (dna, c_count, gc) for every synonymous DNA combination that
    passes all quality filters. Intentionally a generator so the caller can
    stream into heapq for top-K selection.
    """
    codon_choices = [CODONS_BY_AA[aa.upper()] for aa in aa_sequence]
    for combo in itertools.product(*codon_choices):
        dna = "".join(combo)
        if _has_internal_stop(dna):
            continue
        if _has_cpg(dna):
            continue
        if _has_restriction_site(dna):
            continue
        if _has_homopolymer(dna):
            continue
        yield dna, dna.count("C"), _gc_content(dna)


def enumerate_variants(aa_sequence: str, max_variants: int = 100) -> List[CodonVariant]:
    """Return top ``max_variants`` surviving variants, globally ranked.

    Ranking: (c_count DESC, abs(gc_content - 0.55) ASC). Uses :func:`heapq.nlargest`
    with a key that mirrors this ordering so the top-K is globally correct
    even when only a subset of combinations is kept in memory.
    """
    if not aa_sequence:
        return []

    # heapq.nlargest keeps the K biggest items by ``key``.  Our desired order
    # is (c_count DESC, |gc - 0.55| ASC) so a larger key = (c_count, -|gc-0.55|).
    top = heapq.nlargest(
        max_variants,
        _surviving_variants(aa_sequence),
        key=lambda t: (t[1], -abs(t[2] - 0.55)),
    )

    # nlargest produces items in descending key order; but ties may land in
    # arbitrary order.  Re-sort exactly by the canonical ordering before
    # assigning ranks so the result is deterministic.
    top.sort(key=lambda t: (-t[1], abs(t[2] - 0.55), t[0]))

    return [
        CodonVariant(dna=dna, c_count=c, gc_content=gc, rank=i)
        for i, (dna, c, gc) in enumerate(top)
    ]


# ---------------------------------------------------------- bulk optimisation


# ---------------------------------------------------------------- filters


_STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})


@lru_cache(maxsize=1)
def load_restriction_sites() -> List[str]:
    """Return the sorted list of recognition sequences bundled with the package."""
    raw = (
        resources.files("pegrna_pipeline.data")
        .joinpath("restriction_sites.yaml")
        .read_text()
    )
    data = yaml.safe_load(raw) or {}
    sites = sorted({e["recognition"].upper() for e in data.get("enzymes", [])})
    return sites


def _has_restriction_site(dna: str) -> bool:
    """Return True if ``dna`` contains any bundled restriction site as a substring."""
    upper = dna.upper()
    return any(site in upper for site in load_restriction_sites())


def _gc_content(dna: str) -> float:
    """Return the GC fraction of ``dna`` in [0.0, 1.0]. Empty input -> 0.0."""
    if not dna:
        return 0.0
    upper = dna.upper()
    return (upper.count("G") + upper.count("C")) / len(upper)


def _has_homopolymer(dna: str, min_run: int = 5) -> bool:
    """Return True if ``dna`` contains a homopolymer run of length >= min_run."""
    if min_run <= 0:
        return bool(dna)
    pattern = re.compile(rf"(.)\1{{{min_run - 1},}}")
    return pattern.search(dna) is not None


def _has_cpg(dna: str) -> bool:
    """Return True if ``dna`` contains a CpG dinucleotide (literal 'CG')."""
    return "CG" in dna


def _has_internal_stop(dna: str) -> bool:
    """Return True if ``dna`` contains any in-frame stop codon.

    Walks indices ``range(0, len(dna) - 2, 3)`` so every full codon is
    inspected -- including the last one (off-by-one fix).
    """
    for i in range(0, len(dna) - 2, 3):
        if dna[i : i + 3] in _STOP_CODONS:
            return True
    return False


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
