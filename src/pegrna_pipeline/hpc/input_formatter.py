"""Format a :class:`GeneRecord` + insertion DNA as PRIDICT 2.0 input.

PRIDICT's bracket notation for an insertion between two bases is::

    AAAAA(+GGG)TTTTT
          ^^^^
          '+' prefix → this is an insertion edit (not a substitution).

This module exposes :func:`format_bracket_insertion` which assembles that
string from a :class:`GeneRecord` and validates the insertion payload.
"""

from __future__ import annotations

import re

from pegrna_pipeline.hpc.data_types import GeneRecord


_ACGT_ONLY = re.compile(r"^[ACGT]+$")


def format_bracket_insertion(gene: GeneRecord, insertion: str) -> str:
    """Return ``f"{gene.upstream_flank}(+{insertion}){gene.downstream_flank}"``.

    Parameters
    ----------
    gene:
        A :class:`GeneRecord` whose ``upstream_flank`` / ``downstream_flank``
        fields will surround the bracketed insertion.
    insertion:
        Uppercase ACGT-only DNA to splice in. Empty strings, lowercase,
        whitespace, ambiguity codes (``N``, ``U``, etc.) all raise
        :class:`ValueError`.

    Raises
    ------
    ValueError
        If ``insertion`` is empty or contains any character outside ACGT.
    """
    if not insertion:
        raise ValueError("insertion must be non-empty")
    if not _ACGT_ONLY.match(insertion):
        raise ValueError(
            f"insertion must be uppercase ACGT only; got {insertion!r}"
        )
    return f"{gene.upstream_flank}(+{insertion}){gene.downstream_flank}"
