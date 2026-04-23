"""Fetch mRNA context for a gene via NCBI Entrez.

Public API
----------
:func:`fetch_gene_record` — resolve a gene symbol to a :class:`GeneRecord`.

Private helpers (patchable in tests)
------------------------------------
:func:`_search_refseq_accession` — `esearch` for a RefSeq mRNA.
:func:`_fetch_genbank_record` — `efetch` the GenBank record.
:func:`_extract_cds_coordinates` — parse CDS feature + stop codon.

All three are module-level so unit tests can monkey-patch them without
hitting the network. Network is only contacted when the real helpers run
(i.e., in integration tests or a user's live invocation).
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple

from pegrna_pipeline.hpc.data_types import GeneRecord

# ------------------------------------------------------------------ errors


class SequenceFetchError(RuntimeError):
    """Base class for all :mod:`sequence_fetcher` failures."""


class GeneNotFoundError(SequenceFetchError):
    """The Entrez gene database returned no hit for the requested symbol."""


class RefSeqNotFoundError(SequenceFetchError):
    """No RefSeq mRNA accession could be resolved for the gene."""


class InsufficientFlankingError(SequenceFetchError):
    """The mRNA is too short to provide the requested flank length."""


# ------------------------------------------------------------ private helpers


def _search_refseq_accession(
    gene_name: str,
    *,
    organism: str = "Homo sapiens",
    api_key: Optional[str] = None,
) -> str:
    """Resolve ``gene_name`` to a RefSeq mRNA accession via esearch.

    Returns the first accession ordered by NCBI relevance. Raises
    :class:`RefSeqNotFoundError` if the search returns no ids.
    """
    from Bio import Entrez  # deferred to keep import cost low

    term = f"{gene_name}[Gene Name] AND {organism}[Organism] AND refseq[Filter] AND mrna[Filter]"
    kwargs = dict(db="nuccore", term=term, retmax=5, sort="relevance")
    if api_key:
        kwargs["api_key"] = api_key
    handle = Entrez.esearch(**kwargs)
    try:
        record = Entrez.read(handle)
    finally:
        handle.close()
    id_list = record.get("IdList", [])
    if not id_list:
        raise RefSeqNotFoundError(f"no RefSeq mRNA for {gene_name!r} in {organism!r}")
    # Resolve first id to its accession.
    sum_kwargs = dict(db="nuccore", id=id_list[0])
    if api_key:
        sum_kwargs["api_key"] = api_key
    sh = Entrez.esummary(**sum_kwargs)
    try:
        summary = Entrez.read(sh)
    finally:
        sh.close()
    if not summary:
        raise RefSeqNotFoundError(f"no summary for id {id_list[0]!r}")
    acc = summary[0].get("AccessionVersion") or summary[0].get("Caption")
    if not acc:
        raise RefSeqNotFoundError(f"no accession field in summary for {gene_name!r}")
    return str(acc)


def _fetch_genbank_record(
    accession: str,
    *,
    email: str,
    api_key: Optional[str] = None,
) -> Dict[str, Any]:
    """Download a GenBank record and return a plain dict:

    ``{"accession": str, "sequence": str, "cds_features": list[(start, end)]}``.

    Indices are 0-based, end-exclusive (Python conventions), NOT GenBank's.
    """
    from Bio import Entrez, SeqIO

    Entrez.email = email
    kwargs = dict(db="nuccore", id=accession, rettype="gb", retmode="text")
    if api_key:
        kwargs["api_key"] = api_key
    handle = Entrez.efetch(**kwargs)
    try:
        record = SeqIO.read(handle, "genbank")
    finally:
        handle.close()
    sequence = str(record.seq).upper().replace("U", "T")
    cds_features: List[Tuple[int, int]] = []
    for feat in record.features:
        if feat.type == "CDS":
            cds_features.append((int(feat.location.start), int(feat.location.end)))
    return {"accession": accession, "sequence": sequence, "cds_features": cds_features}


def _extract_cds_coordinates(record: Dict[str, Any]) -> Tuple[int, int, int, str]:
    """From a record dict, return ``(cds_start, cds_end, stop_codon_start, stop_codon)``.

    Picks the first CDS feature (RefSeq mRNAs have exactly one) and computes
    the stop-codon window as the final 3 nt of the CDS. Raises
    :class:`SequenceFetchError` if no CDS features are present or the stop
    codon is not a canonical TAA/TAG/TGA.
    """
    cds_list = record.get("cds_features", [])
    if not cds_list:
        raise SequenceFetchError("record has no CDS features")
    cds_start, cds_end = cds_list[0]
    sequence = record["sequence"]
    stop_codon_start = cds_end - 3
    stop_codon = sequence[stop_codon_start:cds_end].upper()
    if stop_codon not in {"TAA", "TAG", "TGA"}:
        raise SequenceFetchError(f"final CDS codon is {stop_codon!r}, expected TAA/TAG/TGA")
    return cds_start, cds_end, stop_codon_start, stop_codon


# ------------------------------------------------------------------ public


def fetch_gene_record(
    gene_name: str,
    *,
    email: str,
    organism: str = "Homo sapiens",
    flanking_bp: int = 100,
    api_key: Optional[str] = None,
) -> GeneRecord:
    """Resolve ``gene_name`` to a :class:`GeneRecord` with flanking context.

    Raises
    ------
    RefSeqNotFoundError
        No RefSeq mRNA matched.
    GeneNotFoundError
        Explicit upstream "no such gene" from Entrez.
    InsufficientFlankingError
        The mRNA is too short for the requested flank length.
    SequenceFetchError
        Any other problem (bad CDS, network failure, etc.).
    """
    accession = _search_refseq_accession(gene_name, organism=organism, api_key=api_key)
    record = _fetch_genbank_record(accession, email=email, api_key=api_key)
    cds_start, cds_end, stop_codon_start, stop_codon = _extract_cds_coordinates(record)
    seq = record["sequence"]

    if stop_codon_start - flanking_bp < 0:
        raise InsufficientFlankingError(
            f"upstream flank ({flanking_bp} bp) exceeds available sequence "
            f"before stop codon at {stop_codon_start}"
        )
    if stop_codon_start + 3 + flanking_bp > len(seq):
        raise InsufficientFlankingError(
            f"downstream flank ({flanking_bp} bp) exceeds available sequence "
            f"after stop codon end at {stop_codon_start + 3}"
        )

    upstream = seq[stop_codon_start - flanking_bp : stop_codon_start]
    downstream = seq[stop_codon_start + 3 : stop_codon_start + 3 + flanking_bp]

    return GeneRecord(
        gene_name=gene_name,
        refseq_accession=accession,
        mrna_sequence=seq,
        cds_start=cds_start,
        cds_end=cds_end,
        stop_codon_start=stop_codon_start,
        stop_codon=stop_codon,
        upstream_flank=upstream,
        downstream_flank=downstream,
        flanking_bp=flanking_bp,
    )


# ---------------------------------------------------------------- __main__


def _main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m pegrna_pipeline.hpc.sequence_fetcher",
        description="Fetch mRNA context for a gene via NCBI Entrez.",
    )
    parser.add_argument("gene", help="Gene symbol (e.g. CTCF).")
    parser.add_argument("--email", required=True, help="Email for Entrez (required by NCBI).")
    parser.add_argument(
        "--organism",
        default="Homo sapiens",
        help="Organism name passed to Entrez (default: Homo sapiens).",
    )
    parser.add_argument("--flanking-bp", type=int, default=100, help="Flank length (default: 100).")
    parser.add_argument("--api-key", default=None, help="Optional NCBI API key.")
    parser.add_argument("--json", action="store_true", help="Emit the record as JSON on stdout.")
    args = parser.parse_args(argv)

    record = fetch_gene_record(
        args.gene,
        email=args.email,
        organism=args.organism,
        flanking_bp=args.flanking_bp,
        api_key=args.api_key,
    )
    if args.json:
        sys.stdout.write(json.dumps(dataclasses.asdict(record), indent=2))
        sys.stdout.write("\n")
    else:
        for field, value in dataclasses.asdict(record).items():
            sys.stdout.write(f"{field}: {value}\n")
    return 0


if __name__ == "__main__":  # pragma: no cover - exercised via subprocess
    raise SystemExit(_main())
