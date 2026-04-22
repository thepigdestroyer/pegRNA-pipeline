"""Prime-Assembly (PA) workflow — two-pegRNA construct insertion.

Entry point: :func:`run_pa_workflow`. The donor DNA consists of three
concatenated pieces:

    upstream_flap  +  bulk  +  downstream_flap

Each flap is the portion of the insert that a pegRNA's RT extension lays
down. Plan 2 v1.1 sets the default flap length to 10 amino acids.

Tasks in this module:

* :func:`compute_flap_boundaries` — T15
* :func:`assemble_linear_donor`   — T17
* :class:`PAResult`               — T18-T22
* :func:`run_pa_workflow`         — T18-T22
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union


_ACGT_ONLY = re.compile(r"^[ACGT]+$")


def _check_acgt(label: str, value: str) -> None:
    if not value or not _ACGT_ONLY.match(value):
        raise ValueError(f"{label} must be non-empty uppercase ACGT; got {value!r}")


# ------------------------------------------------------------------ T15


def compute_flap_boundaries(
    full_aa: str,
    *,
    flap_length_aa: int = 10,
) -> Tuple[str, str, str]:
    """Split ``full_aa`` into (upstream_flap, bulk, downstream_flap).

    Each flap has exactly ``flap_length_aa`` characters. Raises
    :class:`ValueError` if ``len(full_aa) < 2 * flap_length_aa`` — you
    can't have two flaps of the requested size otherwise.
    """
    required = 2 * flap_length_aa
    if len(full_aa) < required:
        raise ValueError(
            f"full_aa is {len(full_aa)} aa; need at least {required} "
            f"(2 * flap_length_aa={flap_length_aa})"
        )
    upstream = full_aa[:flap_length_aa]
    downstream = full_aa[-flap_length_aa:]
    bulk = full_aa[flap_length_aa : len(full_aa) - flap_length_aa]
    return upstream, bulk, downstream


# ------------------------------------------------------------------ T17


def assemble_linear_donor(
    *,
    upstream_flap_dna: str,
    bulk_dna: str,
    downstream_flap_dna: str,
) -> str:
    """Return ``upstream_flap_dna + bulk_dna + downstream_flap_dna``.

    Every component must be non-empty uppercase ACGT — mixed case, U, or
    ambiguity codes raise :class:`ValueError`.
    """
    _check_acgt("upstream_flap_dna", upstream_flap_dna)
    _check_acgt("bulk_dna", bulk_dna)
    _check_acgt("downstream_flap_dna", downstream_flap_dna)
    return upstream_flap_dna + bulk_dna + downstream_flap_dna
