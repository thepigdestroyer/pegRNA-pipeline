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

# PRIDICT2 editseq bracket conventions (batch mode): a deletion is written
# ``FLANK(-DELETEDSEQ)FLANK``, a substitution ``FLANK(REF/ALT)FLANK``, and an
# INSERTION ``FLANK(+INSERTEDSEQ)FLANK``.
#
# CONFIRMED against a live PRIDICT2 run on the DSI box (2026-07-10): ``(+SEQ)``
# produces a valid pegRNA; ``(/SEQ)`` errors ("string index out of range"). This
# matches the legacy ``format_bracket_insertion`` marker. Isolated in the two
# constants below so a single edit corrects every editseq the pipeline emits.
_INSERTION_OPEN = "(+"
_INSERTION_CLOSE = ")"


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
        raise ValueError(f"insertion must be uppercase ACGT only; got {insertion!r}")
    return f"{gene.upstream_flank}(+{insertion}){gene.downstream_flank}"


def format_tjpe_insertion_editseq(
    genome: str,
    site_s: int,
    insert: str,
    *,
    min_flank: int = 100,
) -> str:
    """Return a PRIDICT2 insertion editseq for a tj-PE insertion at ``site_s``.

    The result is ``genome[:site_s] + "(+" + insert + ")" + genome[site_s:]``,
    uppercased — i.e. the ``(+INSERT)`` insertion bracket spliced into the
    local genome window at ``site_s`` (a half-open index: the insert lands
    *between* ``genome[site_s - 1]`` and ``genome[site_s]``).

    Both flanks (the ``site_s`` bases before the bracket and the
    ``len(genome) - site_s`` bases after it) must be at least ``min_flank`` bp
    so PRIDICT has adequate sequence context on each side.

    Parameters
    ----------
    genome:
        Local genome window (ACGT; case-insensitive, uppercased on return).
    site_s:
        Half-open insertion index within ``[0, len(genome)]``.
    insert:
        Uppercase ACGT-only DNA to insert (case-insensitive, uppercased on
        return). Empty strings and any non-ACGT character raise
        :class:`ValueError`.
    min_flank:
        Minimum required flank length on *each* side of the insertion
        (default 100 bp).

    Raises
    ------
    ValueError
        If ``insert`` is empty or non-ACGT, ``site_s`` is outside
        ``[0, len(genome)]``, or either flank is shorter than ``min_flank``.
    """
    if not insert:
        raise ValueError("insert must be non-empty")
    insert_up = insert.upper()
    if not _ACGT_ONLY.match(insert_up):
        raise ValueError(f"insert must be ACGT only; got {insert!r}")

    genome_up = genome.upper()
    genome_len = len(genome_up)
    if site_s < 0 or site_s > genome_len:
        raise ValueError(f"site_s must be within [0, {genome_len}]; got {site_s}")

    left_flank = site_s
    right_flank = genome_len - site_s
    if left_flank < min_flank or right_flank < min_flank:
        raise ValueError(
            f"genome window too short for PRIDICT context: need >= {min_flank} bp "
            f"each side of site_s={site_s}, got left={left_flank} right={right_flank}"
        )

    return (
        f"{genome_up[:site_s]}{_INSERTION_OPEN}{insert_up}{_INSERTION_CLOSE}"
        f"{genome_up[site_s:]}"
    )
