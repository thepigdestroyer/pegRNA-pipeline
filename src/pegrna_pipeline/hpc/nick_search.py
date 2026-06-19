"""Enumerate candidate nicking-guide sites for the TJ-PE second nick.

Candidates lie in the ``[distance_min, distance_max]`` window (default 40-70 bp)
on the strand OPPOSITE the pegRNA, each anchored to an ``NGG`` PAM with a 20-nt
protospacer. Off-target / on-target scoring is applied separately
(:func:`screen_nick_candidates`, Task 6) to the protospacers returned here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple

from Bio.Seq import Seq

__all__ = ["enumerate_nick_sites", "screen_nick_candidates", "NickCandidate"]


@dataclass(frozen=True)
class NickCandidate:
    protospacer: str       # 20 nt
    pam: str               # 3 nt, NGG
    nick_position: int     # 0-indexed on the sense genome
    distance_to_edit: int  # |nick - S|
    target_strand: str     # "Fw" or "Rv" (opposite the pegRNA)


def _opposite(strand: str) -> str:
    return {"Fw": "Rv", "Rv": "Fw"}[strand]


def enumerate_nick_sites(
    genome: str, *, site_s: int, pegrna_strand: str,
    distance_min: int, distance_max: int,
) -> List[NickCandidate]:
    """Return PAM-anchored nick candidates within the distance window on the
    strand opposite ``pegrna_strand``.

    SpCas9 nicks ~3 nt 5' of the PAM; the protospacer is the 20 nt immediately
    5' of the PAM. Positions are reported in sense-genome coordinates.
    """
    genome = genome.upper()
    want = _opposite(pegrna_strand)
    n = len(genome)
    out: List[NickCandidate] = []

    def scan(seq: str, strand: str, to_sense: Callable[[int], int]) -> None:
        for i in range(len(seq) - 2):
            if seq[i + 1 : i + 3] != "GG":  # NGG PAM (positions i+1, i+2)
                continue
            proto_end = i                  # protospacer is the 20 nt 5' of the PAM
            proto_start = proto_end - 20
            if proto_start < 0:
                continue
            nick = to_sense(proto_end - 3)
            dist = abs(nick - site_s)
            if distance_min <= dist <= distance_max:
                out.append(
                    NickCandidate(
                        protospacer=seq[proto_start:proto_end],
                        pam=seq[i : i + 3],
                        nick_position=nick,
                        distance_to_edit=dist,
                        target_strand=strand,
                    )
                )

    if want == "Fw":
        scan(genome, "Fw", lambda x: x)
    else:
        rc = str(Seq(genome).reverse_complement())
        scan(rc, "Rv", lambda x: n - 1 - x)
    # Deterministic ordering: nearest first, then protospacer.
    out.sort(key=lambda c: (c.distance_to_edit, c.protospacer))
    return out


def screen_nick_candidates(
    candidates: List[NickCandidate], *,
    scorer: Callable[[List[str]], Dict[str, Tuple[float, int, bool]]],
    min_ontarget: float = 0.5,
    max_offtargets: int = 0,
    other_allele_seqs: Dict[str, str] | None = None,
) -> List[NickCandidate]:
    """Keep PAM-valid candidates that are off-target-clean and active.

    ``scorer`` maps protospacer -> ``(ontarget_score, offtarget_count,
    has_offtarget)``; inject the real Cas9 off-target tool
    (``offtarget_scorer.submit_offtarget_job``) or a fake in tests.

    ``other_allele_seqs`` (decision 13) maps protospacer -> the homologous
    sequence on the OTHER haplotype; when supplied, a candidate whose
    protospacer still matches the other allele exactly is treated as having an
    on-locus off-target and dropped (an allele-specific guide must NOT cut the
    other allele). Returns survivors sorted by descending on-target score.
    """
    scores = scorer([c.protospacer for c in candidates])
    other = other_allele_seqs or {}
    kept: List[Tuple[float, NickCandidate]] = []
    for c in candidates:
        on, off, has_off = scores.get(c.protospacer, (0.0, 10**9, True))
        if has_off or off > max_offtargets or on < min_ontarget:
            continue
        if other.get(c.protospacer) == c.protospacer:
            continue  # other allele identical => not allele-discriminating, would cut both
        kept.append((on, c))
    kept.sort(key=lambda t: (-t[0], t[1].distance_to_edit, t[1].protospacer))
    return [c for _, c in kept]
