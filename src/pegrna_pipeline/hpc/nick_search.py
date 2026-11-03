"""Enumerate candidate nicking-guide sites for the TJ-PE second nick.

Candidates lie in the ``[distance_min, distance_max]`` window (default 40-70 bp)
on the strand OPPOSITE the pegRNA, each anchored to an ``NGG`` PAM with a 20-nt
protospacer. Off-target / on-target scoring is applied separately
(:func:`screen_nick_candidates`, Task 6) to the protospacers returned here.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable, Dict, List, Optional, Tuple

from Bio.Seq import Seq

from pegrna_pipeline.hpc.guide_sanity import assess_spacer

__all__ = [
    "enumerate_nick_sites",
    "filter_nick_sanity",
    "offtarget_rejection",
    "screen_nick_candidates",
    "nick_usable_for_arm",
    "NickCandidate",
]

ARM_SIDES = ("3prime", "5prime", "auto")


def nick_usable_for_arm(nick_position: int, site_s: int, arm_side: str) -> bool:
    """Can a nick at ``nick_position`` serve a design whose arm sits on ``arm_side``?

    ``arm_side`` selects which side of the insertion site S the recode/replacement
    arm occupies, and the second nick must lie on that same side: ``"5prime"`` needs
    the nick 5' of S (the arm is carved out of the CDS upstream of the stop codon),
    ``"3prime"`` needs it 3' of S, and ``"auto"`` accepts either. A nick exactly at S
    has no arm and is never usable.

    This is the single definition of that rule. It has to be applied wherever a
    candidate pool is BUILT, not just where it is walked: a pool whose every member
    the walk will discard is not a pool, and must not suppress the fallback to
    other candidate sources.
    """
    if arm_side not in ARM_SIDES:
        raise ValueError(f"arm_side must be one of {ARM_SIDES}; got {arm_side!r}")
    if nick_position == site_s:
        return False
    if arm_side == "5prime":
        return nick_position < site_s
    if arm_side == "3prime":
        return nick_position > site_s
    return True


@dataclass(frozen=True)
class NickCandidate:
    protospacer: str  # 20 nt
    pam: str  # 3 nt, NGG
    nick_position: int  # 0-indexed on the sense genome
    distance_to_edit: int  # |nick - S|
    target_strand: str  # "Fw" or "Rv" (opposite the pegRNA)
    #: Off-target sites the screen counted for this guide (extra exact copies plus near
    #: matches within the mismatch tier). ``None`` until the candidate has been screened.
    offtargets: Optional[int] = None


def offtarget_rejection(
    score: Tuple[float, int, bool], *, max_offtargets: int, min_ontarget: float
) -> Optional[str]:
    """Why a guide's ``(on-target, off-target count, flagged)`` triple disqualifies it, or
    ``None`` when it is acceptable.

    The count is what is compared with ``max_offtargets``. The ``flagged`` bit is redundant
    with it (the scorers set it to ``count > 0``), and it must not be consulted on its own:
    doing so made ``max_offtargets`` a no-op for every value above 0, because any guide with
    a single off-target was already flagged. A guide flagged with a count of zero is
    contradictory, and is rejected rather than trusted.

    Exact duplicates of the guide elsewhere in the genome are never tolerated, whatever the
    limit: ``max_offtargets`` is a tolerance for near-matches (a cutting site that is weak and
    rare), and a perfect copy is neither. Scorers that report the duplicate count do so on the
    score (``OffTargetScore.duplicates``); one that returns a plain tuple has none known.

    This is the one definition of the rule; both the genome-enumerated and the PRIDICT paths
    call it, so they cannot disagree about the same guide.
    """
    on, off, flagged = score
    if on < min_ontarget:
        return f"on-target score {on:.2f} is below {min_ontarget:g}"
    duplicates = getattr(score, "duplicates", 0)
    if duplicates > 0:
        return (
            f"{duplicates} exact copy/copies of this guide elsewhere in the genome: it is not "
            "unique to this locus (never tolerated)"
        )
    if off > max_offtargets:
        return f"{off} off-target site(s), more than the {max_offtargets} allowed"
    if flagged and off <= 0:
        return "the scorer flagged an off-target but reported no count"
    return None


def filter_nick_sanity(
    candidates: List["NickCandidate"],
) -> Tuple[List["NickCandidate"], List[Tuple["NickCandidate", Tuple[str, ...]]]]:
    """Split candidates into those whose spacer breaks no HARD sequence rule (TTTT, GC) and
    those that do, with the reasons. Applied before the genome scan, which is the expensive
    step: there is no point scanning a guide that cannot be used. This is the same check the
    PRIDICT path applies (``guide_sanity``), so a guide is judged identically on both."""
    kept: List["NickCandidate"] = []
    dropped: List[Tuple["NickCandidate", Tuple[str, ...]]] = []
    for c in candidates:
        verdict = assess_spacer(c.protospacer)
        if verdict.ok:
            kept.append(c)
        else:
            dropped.append((c, verdict.failures))
    return kept, dropped


def _opposite(strand: str) -> str:
    return {"Fw": "Rv", "Rv": "Fw"}[strand]


def enumerate_nick_sites(
    genome: str,
    *,
    site_s: int,
    pegrna_strand: str,
    distance_min: int,
    distance_max: int,
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
            proto_end = i  # protospacer is the 20 nt 5' of the PAM
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
    candidates: List[NickCandidate],
    *,
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
        score = scores.get(c.protospacer, (0.0, 10**9, True))
        if offtarget_rejection(score, max_offtargets=max_offtargets, min_ontarget=min_ontarget):
            continue
        if other.get(c.protospacer) == c.protospacer:
            continue  # other allele identical => not allele-discriminating, would cut both
        kept.append((score[0], replace(c, offtargets=int(score[1]))))
    kept.sort(key=lambda t: (-t[0], t[1].distance_to_edit, t[1].protospacer))
    return [c for _, c in kept]
