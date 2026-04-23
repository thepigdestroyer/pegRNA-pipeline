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

from pegrna_pipeline.hpc.data_types import (
    GeneRecord,
    PegRNADesign,
    VariantScore,
)
from pegrna_pipeline.hpc.pair_scorer import (
    PegRNAPair,
    filter_compatible_pairs,
)
from pegrna_pipeline.hpc.pridict_parser import (
    parse_ngrna_csv,
    parse_pegrna_csv,
)


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


# -------------------------------------------------------------- T18-T22


@dataclass
class PAResult:
    """Final output of :func:`run_pa_workflow`.

    Fields follow Plan 2 v1.1 exactly. ``winning_pair`` is the top pair
    by joint score; ``top_pairs`` carries the full ranked list so the
    caller can pick a backup. Off-target fields are reserved for Plan 3
    (CRISPR specificity scoring) — left as empty strings here.
    """

    winning_pair: Optional[PegRNAPair]
    top_pairs: List[PegRNAPair]
    linear_donor_dna: str
    upstream_pegRNA_sequence: str
    downstream_pegRNA_sequence: str
    upstream_flap_dna: str
    downstream_flap_dna: str
    bulk_dna: str
    full_integrated_product: str
    upstream_pegRNA_off_target: str
    downstream_pegRNA_off_target: str
    total_pairs_evaluated: int
    pairs_after_compat_filter: int


def _load_pegrna_designs_from_fixtures(fixtures_dir: Path) -> List[PegRNADesign]:
    """Flatten every ``*_pegRNA_Pridict_full.csv`` in ``fixtures_dir``
    into one list of :class:`PegRNADesign`.
    """
    fixtures_dir = Path(fixtures_dir)
    designs: List[PegRNADesign] = []
    for peg_path in sorted(fixtures_dir.glob("*_pegRNA_Pridict_full.csv")):
        designs.extend(parse_pegrna_csv(peg_path))
    return designs


def _resolve_flap_designs(
    *,
    flap_name: str,
    fixtures_dir: Optional[Union[str, Path]],
    variant_scores: Optional[Sequence[VariantScore]],
) -> List[PegRNADesign]:
    """Return a flat list of pegRNA designs for one flap, enforcing the
    "exactly one source" contract."""
    if (fixtures_dir is None) == (variant_scores is None):
        raise ValueError(
            f"exactly one of {flap_name}_fixtures_dir or {flap_name}_variant_scores "
            "must be provided"
        )
    if fixtures_dir is not None:
        return _load_pegrna_designs_from_fixtures(Path(fixtures_dir))
    # variant_scores path: flatten per-variant design lists.
    designs: List[PegRNADesign] = []
    for vs in variant_scores or []:
        designs.extend(vs.pegRNA_designs)
    return designs


def run_pa_workflow(
    *,
    gene: GeneRecord,
    upstream_fixtures_dir: Optional[Union[str, Path]] = None,
    downstream_fixtures_dir: Optional[Union[str, Path]] = None,
    upstream_variant_scores: Optional[Sequence[VariantScore]] = None,
    downstream_variant_scores: Optional[Sequence[VariantScore]] = None,
    upstream_flap_dna: Optional[str] = None,
    downstream_flap_dna: Optional[str] = None,
    bulk_dna: Optional[str] = None,
    donor_dna: Optional[str] = None,
    flap_length_aa: int = 10,  # noqa: ARG001 — passed through for Plan 3
    min_pair_distance: int = 30,
    max_pair_distance: int = 80,
) -> PAResult:
    """Run the Prime Assembly workflow for one gene + flap pair.

    Exactly one of ``{upstream,downstream}_fixtures_dir`` or
    ``{upstream,downstream}_variant_scores`` must be supplied **per flap**
    (so you can mix fixture-driven and programmatic paths in Plan 3 tests).

    Donor DNA policy: if ``donor_dna`` is supplied, it is used verbatim
    (nothing else is assembled). Otherwise, ``upstream_flap_dna +
    bulk_dna + downstream_flap_dna`` are concatenated via
    :func:`assemble_linear_donor`.
    """
    upstream_designs = _resolve_flap_designs(
        flap_name="upstream",
        fixtures_dir=upstream_fixtures_dir,
        variant_scores=upstream_variant_scores,
    )
    downstream_designs = _resolve_flap_designs(
        flap_name="downstream",
        fixtures_dir=downstream_fixtures_dir,
        variant_scores=downstream_variant_scores,
    )

    total_evaluated = len(upstream_designs) * len(downstream_designs)
    compat_pairs = filter_compatible_pairs(
        upstream=upstream_designs,
        downstream=downstream_designs,
        stop_codon_position=gene.stop_codon_start,
        min_distance=min_pair_distance,
        max_distance=max_pair_distance,
    )
    winning = compat_pairs[0] if compat_pairs else None

    # Donor assembly
    if donor_dna is not None:
        linear_donor = donor_dna
        final_upstream_flap = upstream_flap_dna or ""
        final_bulk = bulk_dna or ""
        final_downstream_flap = downstream_flap_dna or ""
    else:
        if not (upstream_flap_dna and bulk_dna and downstream_flap_dna):
            raise ValueError(
                "when donor_dna is not supplied, upstream_flap_dna, bulk_dna, "
                "and downstream_flap_dna must all be provided"
            )
        linear_donor = assemble_linear_donor(
            upstream_flap_dna=upstream_flap_dna,
            bulk_dna=bulk_dna,
            downstream_flap_dna=downstream_flap_dna,
        )
        final_upstream_flap = upstream_flap_dna
        final_bulk = bulk_dna
        final_downstream_flap = downstream_flap_dna

    upstream_full = winning.upstream.pegrna_full if winning else ""
    downstream_full = winning.downstream.pegrna_full if winning else ""

    # Represent the "full integrated product" as the donor dropped into the
    # mRNA flanks at the stop codon — a convenient integration check for
    # downstream users / tests.
    full_integrated = gene.upstream_flank + linear_donor + gene.downstream_flank

    return PAResult(
        winning_pair=winning,
        top_pairs=list(compat_pairs),
        linear_donor_dna=linear_donor,
        upstream_pegRNA_sequence=upstream_full,
        downstream_pegRNA_sequence=downstream_full,
        upstream_flap_dna=final_upstream_flap,
        downstream_flap_dna=final_downstream_flap,
        bulk_dna=final_bulk,
        full_integrated_product=full_integrated,
        upstream_pegRNA_off_target="",
        downstream_pegRNA_off_target="",
        total_pairs_evaluated=total_evaluated,
        pairs_after_compat_filter=len(compat_pairs),
    )
