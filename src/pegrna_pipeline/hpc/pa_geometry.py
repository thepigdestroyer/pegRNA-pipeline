"""Pure geometry of a Prime Assembly design: guides, cuts, arms, the donor, the pegRNAs. No models.

Prime Assembly (Liu et al., Nature 2026) writes two ~35-nt flaps with two pegRNAs on opposite strands, and the
ends of a separate linear donor match those flaps. The DNA between the two nicks is replaced by the donor. In
this project (decision A+) the donor also RESTORES that DNA, the way the Replacement Arm does in tj-PE: coding
stretches are recoded synonymously, non-coding stretches are byte-exact, and the net deletion is zero. The
whole picture is drawn in ``docs/STRUCTURE.md`` and specified in
``docs/superpowers/specs/2026-11-04-prime-assembly-design.md``.

Conventions (``docs/STRUCTURE.md`` section 0). The window is the SENSE strand, 0-indexed; everything is written
5'->3'; antiparallel relationships use ``reverse_complement``. A **cut** ``b`` is the index of the first base
to the RIGHT of the cut, ``genome[:b] | genome[b:]``, for a guide on either strand. SpCas9 nicks the PAM strand
3 bp 5' of the PAM and the flap grows 3' from the nick, toward the PAM, so for two flaps to fill the gap the
top-strand guide **A** must sit left of the bottom-strand guide **B** ("PAM-in"): both PAMs lie between the nicks.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

from Bio.Seq import Seq

from pegrna_pipeline.hpc.codon_expander import codon_frequency
from pegrna_pipeline.hpc.microhomology import longest_common_run
from pegrna_pipeline.hpc.recoded_tjpe_workflow import _frame_at, _partition_arm
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD

__all__ = [
    "MIN_NICK_SEPARATION_BP",
    "PAM_FOOTPRINT_CLASH_BP",
    "Guide",
    "GuidePair",
    "Arm",
    "PamOutcome",
    "Extension",
    "find_guides",
    "footprints_clear",
    "pam_in_pairs",
    "build_arm",
    "break_pams",
    "assemble_donor",
    "product_allele",
    "flap_lengths",
    "derive_extension",
    "choose_extension",
    "genomic_run_under_flap",
    "assemble_pa_pegrna",
    "pa_regions",
    "nick_continuity_pairs",
    "SPACER_BASES_RIGHT_OF_NICK",
]

#: Hard floor on the distance between the two nicks (the project's value since the original plan: closer
#: than this and the two editing events interfere). The geometric minimum is smaller, see below.
MIN_NICK_SEPARATION_BP = 23
#: Two guides' PAM-bearing footprints overlap below this span: the two Cas9s cannot both bind.
PAM_FOOTPRINT_CLASH_BP = 12

SPACER_LEN = 20
_DNA = frozenset("ACGT")


def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _require_dna(label: str, seq: str, *, allow_empty: bool = False) -> None:
    if (not seq and not allow_empty) or set(seq.upper()) - _DNA:
        raise ValueError(f"{label} must be uppercase ACGT; got {seq!r}")


# ---- guides, cuts, pairs ----------------------------------------------------------------------


@dataclass(frozen=True)
class Guide:
    """One NGG-bearing guide and where it cuts (``cut`` is the first base right of the cut)."""

    strand: str  # "top" | "bottom"
    spacer: str  # 20 nt, 5'->3', as written in the guide RNA
    pam: str  # 3 nt, 5'->3' on the PAM strand
    cut: int
    pam_span: Tuple[int, int]  # sense coordinates [lo, hi)
    protospacer_span: Tuple[int, int]  # sense coordinates [lo, hi)

    @property
    def footprint(self) -> Tuple[int, int]:
        lo = min(self.protospacer_span[0], self.pam_span[0])
        hi = max(self.protospacer_span[1], self.pam_span[1])
        return lo, hi

    @property
    def critical_bases(self) -> Tuple[Tuple[int, int], str]:
        """The two sense positions whose identity makes this a PAM (the GG of NGG), and the base
        they must hold in the sense window: G for a top guide, C for a bottom guide (the complement).
        """
        lo = self.pam_span[0]
        return ((lo + 1, lo + 2), "G") if self.strand == "top" else ((lo, lo + 1), "C")


def find_guides(genome: str) -> List[Guide]:
    """Every NGG-bearing guide on both strands, in a stable order (by cut, then strand, then spacer).

    A top-strand guide has its PAM ``genome[i:i+3] = NGG`` and protospacer ``genome[i-20:i]``; its cut is
    ``i - 3``. A bottom-strand guide shows ``CCN`` at ``genome[j:j+3]``; its protospacer is
    ``RC(genome[j+3:j+23])`` and its cut ``j + 6``. A protospacer or PAM with an ambiguous base is skipped.
    """
    g = genome.upper()
    out: List[Guide] = []
    for i in range(SPACER_LEN, len(g) - 2):
        if g[i + 1 : i + 3] == "GG":
            spacer, pam = g[i - SPACER_LEN : i], g[i : i + 3]
            if not set(spacer + pam) - _DNA:
                out.append(Guide("top", spacer, pam, i - 3, (i, i + 3), (i - SPACER_LEN, i)))
    for j in range(0, len(g) - SPACER_LEN - 2):
        if g[j : j + 2] == "CC":
            proto, pam_sense = g[j + 3 : j + 3 + SPACER_LEN], g[j : j + 3]
            if not set(proto + pam_sense) - _DNA:
                out.append(
                    Guide(
                        "bottom",
                        _rc(proto),
                        _rc(pam_sense),
                        j + 6,
                        (j, j + 3),
                        (j + 3, j + 3 + SPACER_LEN),
                    )
                )
    out.sort(key=lambda x: (x.cut, x.strand, x.spacer))
    return out


def footprints_clear(a: Guide, b: Guide) -> bool:
    """True when the two guides' PAM-bearing footprints do not overlap (the span is at least 12)."""
    return a.footprint[1] <= b.footprint[0]


@dataclass(frozen=True)
class GuidePair:
    a: Guide  # top strand, 5' cut n1
    b: Guide  # bottom strand, 3' cut n2

    @property
    def span(self) -> int:
        return self.b.cut - self.a.cut


def pam_in_pairs(
    guides: Sequence[Guide], *, site_s: int, min_span: int, max_span: int
) -> List[GuidePair]:
    """Every PAM-in pair with S between the nicks and a span in ``[min_span, max_span]``, smallest span first.

    Guide A is on the top strand with ``cut <= S``; guide B on the bottom strand with ``cut >= S``. Left-A,
    right-B is the only arrangement whose flaps fill the gap, so it is the only one built.
    """
    if min_span < MIN_NICK_SEPARATION_BP:
        raise ValueError(
            f"min_span must be at least {MIN_NICK_SEPARATION_BP} bp (the project floor); got {min_span}"
        )
    tops = [g for g in guides if g.strand == "top" and g.cut <= site_s]
    bottoms = [g for g in guides if g.strand == "bottom" and g.cut >= site_s]
    pairs = [
        GuidePair(a, b)
        for a in tops
        for b in bottoms
        if min_span <= b.cut - a.cut <= max_span and footprints_clear(a, b)
    ]
    pairs.sort(key=lambda p: (p.span, p.a.cut, p.b.cut))
    return pairs


# ---- arms --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Arm:
    """The restored copy of ``genome[lo:hi)``: recoded where coding, byte-exact where not."""

    lo: int
    hi: int
    sense: str
    original: str
    recode_span: Optional[Tuple[int, int]]
    locked_span: Optional[Tuple[int, int]]
    locked_expected: Optional[str]
    n_edits: int
    max_run_after: int
    microhomology_flag: bool  # a non-coding (locked) stretch is present: it cannot be recoded


@lru_cache(maxsize=4096)
def _cached_partition(genome: str, lo: int, hi: int, mask: Tuple[bool, ...], max_run: int):
    return _partition_arm(genome, lo=lo, hi=hi, cds_mask=mask, max_run=max_run)


def build_arm(
    genome: str, lo: int, hi: int, cds_mask: Sequence[bool], max_run: int
) -> Optional[Arm]:
    """The Replacement-Arm logic of tj-PE (``_partition_arm``) applied to ``genome[lo:hi)``: coding runs
    synonymously recoded with no identical stretch over ``max_run``, non-coding runs byte-exact. An empty
    span is an empty arm. ``None`` if a coding run has no synonymous solution or the span crosses more than
    one coding boundary."""
    if hi < lo:
        raise ValueError(f"arm bounds out of order: [{lo}, {hi})")
    if hi == lo:
        return Arm(lo, hi, "", "", None, None, None, 0, 0, False)
    built = _cached_partition(genome, lo, hi, tuple(bool(x) for x in cds_mask), max_run)
    if built is None:
        return None
    return Arm(
        lo=lo,
        hi=hi,
        sense=built.sense_arm,
        original=genome[lo:hi].upper(),
        recode_span=built.recode_span,
        locked_span=built.locked_span,
        locked_expected=built.locked_expected,
        n_edits=built.n_edits,
        max_run_after=built.max_run_after,
        microhomology_flag=built.microhomology_flag,
    )


@dataclass(frozen=True)
class PamOutcome:
    side: str  # "A" | "B"
    #: "broken-by-recode"            the recode already changed a base of the PAM, and no PAM is left
    #: "broken"                      no PAM is left, by one synonymous base change (``position``)
    #: "weakened"                    only a weak PAM (NAG or NGA) is left: no synonymous change gives none, so
    #:                               SpCas9 can still cut it, at a fraction of its NGG rate
    #: "intact-noncoding"            the PAM is regenerated: non-coding DNA is byte-exact (Decision 15)
    #: "intact-no-synonymous-change" the PAM is regenerated: no synonymous change exists there
    status: str
    position: Optional[int]
    #: The PAM of the product as the guide reads it (3 nt, on its own strand).
    pam: str = ""


#: The last two bases of the PAMs that wild-type SpCas9 still cuts, at a fraction of its NGG rate. NAG: about a fifth
#: of NGG (Hsu et al. 2013). NGA: 0.08-0.14 of the cleaved PAM reads in a human-cell PAM-library screen, against
#: 0.21-0.23 for NAG (Tang et al. 2019, three targets; a share of the cleavage events, not a rate relative to NGG), and
#: over 20% indels at 3 of 19 NGA sites against 17 of 17 NGG sites (Nishimasu et al. 2018). Every other change to a G
#: of NGG leaves no PAM: that second study found no such indels at 0 of 18 NGT and 0 of 15 NGC sites. So when a
#: synonymous change can give a C or a T, it is preferred over an A (ADR 0030).
_WEAK_PAM_TAILS = frozenset({"AG", "GA"})


def _product_pam(guide: Guide, restored: Sequence[str], n1: int, n2: int, genome: str) -> str:
    """The guide's PAM in the restored DNA, 3 nt as the guide reads it (``restored`` is aligned to ``genome[n1:n2]``)."""
    lo, hi = guide.pam_span
    sense = "".join(restored[q - n1] if n1 <= q < n2 else genome[q] for q in range(lo, hi))
    return sense if guide.strand == "top" else _rc(sense)


def _with_edits(arm: Arm, sense: str) -> Arm:
    return replace(
        arm,
        sense=sense,
        n_edits=sum(1 for a, b in zip(sense, arm.original) if a != b),
        max_run_after=longest_common_run(sense, arm.original) if sense else 0,
    )


def break_pams(
    genome: str,
    arm5: Arm,
    arm3: Arm,
    guides: Sequence[Tuple[str, Guide]],
    cds_mask: Sequence[bool],
    max_run: int,
) -> Tuple[Arm, Arm, List[PamOutcome]]:
    """Break each guide's PAM in the restored DNA by one synonymous base change, where one exists.

    The paper's replacement removes both PAMs; a restoring donor brings them back, so the guides could
    re-engage the product. A PAM can be broken here only inside a coding (recoded) stretch, by changing one
    of the two G's of NGG (C's on the bottom strand) so that the codon still encodes the same amino acid
    and the arm still has no identical stretch over ``max_run``. Either guide's PAM may lie in either arm
    (in the CTCF best pair B's PAM is inside arm5), so this is a pair-level edit. Everything else is
    reported, not forced: non-coding DNA stays byte-exact.

    A change that leaves NAG or NGA does not remove the PAM (see ``_WEAK_PAM_TAILS``), so a change to a C or a T
    is always preferred, and a PAM that can only be taken to NAG or NGA is reported as ``weakened``, never as
    broken. A change that merely undoes the recode's own edit (back to NGG) is not a break.
    """
    n1, s, n2 = arm5.lo, arm5.hi, arm3.hi
    if arm3.lo != s:
        raise ValueError("the two arms must meet at S")
    mask = tuple(bool(x) for x in cds_mask)
    restored = list(arm5.sense + arm3.sense)  # aligned with genome[n1:n2)
    outcomes: List[PamOutcome] = []

    def arm_of(pos: int) -> Arm:
        return arm5 if pos < s else arm3

    def coding(pos: int) -> bool:
        span = arm_of(pos).recode_span
        return span is not None and span[0] <= pos < span[1]

    for side, guide in guides:
        crit, _want = guide.critical_bases
        if any(not n1 <= p < n2 for p in crit):
            raise ValueError("a PAM lies outside the replaced span: not a PAM-in pair")
        tail = _product_pam(guide, restored, n1, n2, genome)[1:]
        if tail != "GG" and tail not in _WEAK_PAM_TAILS:
            outcomes.append(
                PamOutcome(
                    side, "broken-by-recode", None, _product_pam(guide, restored, n1, n2, genome)
                )
            )
            continue
        # The PAM is NGG, or the recode left only a weak one: look for the best single synonymous change.
        best: Optional[Tuple[Tuple[bool, float, int, str], int, str]] = None
        for p in crit:
            if not coding(p):
                continue
            arm = arm_of(p)
            c = p - _frame_at(mask, p)  # first base of the codon holding p
            if c < arm.recode_span[0] or c + 3 > arm.recode_span[1]:
                continue
            old = "".join(restored[c - n1 : c - n1 + 3])
            for alt in "ACGT":
                if alt == restored[p - n1]:
                    continue
                new = old[: p - c] + alt + old[p - c + 1 :]
                aa_old, aa_new = str(Seq(old).translate()), str(Seq(new).translate())
                if aa_old != aa_new or "*" in aa_new:
                    continue
                trial = restored[:]
                trial[p - n1] = alt
                after = _product_pam(guide, trial, n1, n2, genome)[1:]
                if after == "GG":
                    continue  # this only undoes the recode's own change: the PAM is back
                lo, hi = (0, s - n1) if arm is arm5 else (s - n1, n2 - n1)
                # The edit may only break identity, never extend it past what the arm already allowed.
                allowed = max(max_run, longest_common_run("".join(restored[lo:hi]), arm.original))
                if longest_common_run("".join(trial[lo:hi]), arm.original) > allowed:
                    continue
                key = (after in _WEAK_PAM_TAILS, -codon_frequency(aa_new, new), p, alt)
                if best is None or key < best[0]:
                    best = (key, p, alt)
        if best is not None:
            key, p, alt = best
            restored[p - n1] = alt
            status = "weakened" if key[0] else "broken"
            outcomes.append(
                PamOutcome(side, status, p, _product_pam(guide, restored, n1, n2, genome))
            )
            continue
        now = _product_pam(guide, restored, n1, n2, genome)
        if tail != "GG":
            status = "weakened"  # the recode left NAG or NGA and nothing synonymous improves on it
        elif any(coding(p) for p in crit):
            status = "intact-no-synonymous-change"
        else:
            status = "intact-noncoding"
        outcomes.append(PamOutcome(side, status, None, now))

    new5, new3 = "".join(restored[: s - n1]), "".join(restored[s - n1 :])
    return (
        arm5 if new5 == arm5.sense else _with_edits(arm5, new5),
        arm3 if new3 == arm3.sense else _with_edits(arm3, new3),
        outcomes,
    )


# ---- donor and product -------------------------------------------------------------------------


def assemble_donor(arm5: str, cargo: str, arm3: str) -> str:
    """``arm5 + cargo + arm3``: the restored DNA either side of the cargo, which sits at S."""
    _require_dna("arm5", arm5, allow_empty=True)
    _require_dna("cargo", cargo)
    _require_dna("arm3", arm3, allow_empty=True)
    return arm5 + cargo + arm3


def product_allele(genome: str, n1: int, n2: int, donor: str) -> str:
    """The genome after assembly: the donor stands where ``genome[n1:n2)`` was."""
    return genome[:n1] + donor + genome[n2:]


# ---- flaps, extensions, the pegRNA -------------------------------------------------------------


def flap_lengths(default: int = 35, lo: int = 30, hi: int = 50) -> List[int]:
    """Flap lengths in the order they are tried: the default, then outward (+1, -1, +2, -2, ...) while
    inside ``[lo, hi]``, then the rest upward."""
    out = [default]
    for k in range(1, hi - lo + 2):
        for cand in (default + k, default - k):
            if lo <= cand <= hi and cand not in out:
                out.append(cand)
    return out


@dataclass(frozen=True)
class Extension:
    side: str  # "A" | "B"
    flap_len: int
    pbs_len: int
    pbs: str
    rtt: str


def derive_extension(
    side: str, genome: str, cut: int, donor: str, flap_len: int, pbs_len: int
) -> Optional[Extension]:
    """The PBS and RTT of one guide. Guide A (top strand, cut ``n1``): ``PBS = RC(genome[n1-L:n1])`` pairs
    with the 3' end of the nicked top strand and ``RTT = RC(donor[:F])`` writes the donor's start. Guide B
    (bottom strand, cut ``n2``): ``PBS = genome[n2:n2+L]`` and ``RTT = donor[-F:]``. ``None`` when the PBS
    window leaves the genome, the flap outruns the donor, or a base is not ACGT."""
    if side not in ("A", "B"):
        raise ValueError(f"side must be 'A' or 'B'; got {side!r}")
    if flap_len > len(donor):
        return None
    if side == "A":
        if cut - pbs_len < 0:
            return None
        pbs, rtt = _rc(genome[cut - pbs_len : cut]), _rc(donor[:flap_len])
    else:
        if cut + pbs_len > len(genome):
            return None
        pbs, rtt = genome[cut : cut + pbs_len], donor[-flap_len:]
    pbs, rtt = pbs.upper(), rtt.upper()
    if set(pbs + rtt) - _DNA:
        return None
    return Extension(side, flap_len, pbs_len, pbs, rtt)


def choose_extension(
    side: str,
    genome: str,
    cut: int,
    donor: str,
    pbs_len: int,
    *,
    default: int = 35,
    lo: int = 30,
    hi: int = 50,
) -> Optional[Extension]:
    """Scan the flap length from the default outward: the first length whose extension (the RTT followed by the PBS,
    as transcribed) holds no TTTT (a Pol III terminator) and whose RTT does not start with C (a first-base C pairs
    with the scaffold); failing that, the first length with no TTTT; failing that, ``None``.

    The two are read as one string because the polymerase reads one transcript: an RTT ending in TT followed by a
    PBS starting TT terminates it, though neither part holds TTTT on its own. Note that the RTT's 3' end is fixed by
    the donor's end, so a run formed there is the same at every flap length and no length can avoid it.
    """
    fallback: Optional[Extension] = None
    for flap_len in flap_lengths(default, lo, hi):
        ext = derive_extension(side, genome, cut, donor, flap_len, pbs_len)
        if ext is None or "TTTT" in ext.rtt + ext.pbs:
            continue
        if ext.rtt[0] != "C":
            return ext
        if fallback is None:
            fallback = ext
    return fallback


def genomic_run_under_flap(
    side: str, genome: str, n1: int, n2: int, donor: str, flap_len: int
) -> int:
    """The longest stretch the flap shares with the strand it displaces. A flap that is a copy of the
    genome can re-anneal where it came from instead of on the donor, which is what recoding prevents.
    """
    if side == "A":
        flap, displaced = donor[:flap_len], genome[n1 : min(n1 + flap_len, n2)]
    else:
        flap, displaced = donor[-flap_len:], genome[max(n2 - flap_len, n1) : n2]
    return longest_common_run(flap, displaced)


def assemble_pa_pegrna(
    *,
    spacer: str,
    rtt: str,
    pbs: str,
    tail: str = "UUU",
    scaffold: str = STANDARD_SGRNA_SCAFFOLD,
) -> str:
    """``Spacer · Scaffold · RTT · PBS · tail``, 5'->3': the standard layout, not tj-PE's (there is no
    rcPBS2: the second strand is captured by the donor)."""
    for label, value in (("spacer", spacer), ("scaffold", scaffold), ("rtt", rtt), ("pbs", pbs)):
        _require_dna(label, value)
    if not tail or set(tail.upper()) - {"U"} or not 2 <= len(tail) <= 3:
        raise ValueError(f"tail must be 2-3 uridines (U only); got {tail!r}")
    return spacer + scaffold + rtt + pbs + tail


#: SpCas9 nicks 3 bp from its PAM, so the spacer's last three bases lie to the right of the nick.
SPACER_BASES_RIGHT_OF_NICK = 3


def nick_continuity_pairs(regions: Dict[str, Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    """The position pairs ``(i, j)``, i < j, at which the RTT's 3' end meets the spacer's 3' end.

    The donor copies the genome at each nick, so the RTT's 3' end is the reverse complement of the bases just right
    of the nick. The protospacer is continuous across the nick, so those bases are the spacer's own last three
    (``SPACER_BASES_RIGHT_OF_NICK``), and the spacer-PBS helix runs on into the RTT: spacer base ``17 + m`` against
    the RTT's last-but-``m`` base, for ``m`` = 0, 1, 2. The folding gate already treats the spacer-PBS helix as
    expected; this is the same helix, up to three base pairs longer, and it is inherent to a donor that restores
    the DNA between the nicks (tj-PE starts its insert at the nick, so never meets it).
    """
    first = regions["spacer"][1] - SPACER_BASES_RIGHT_OF_NICK
    rtt_end = regions["rtt"][1]
    return tuple((first + m, rtt_end - 1 - m) for m in range(SPACER_BASES_RIGHT_OF_NICK))


def pa_regions(spacer: str, rtt: str, pbs: str, tail: str = "UUU") -> Dict[str, Tuple[int, int]]:
    """Region offsets in the assembled pegRNA, for the folding gate. The PBS is named ``pbs1`` because the
    shared gate treats the spacer-PBS pair as the intrinsic complementarity it is; the tail is left out, as
    in tj-PE's map."""
    scaffold = len(STANDARD_SGRNA_SCAFFOLD)
    i1 = len(spacer)
    i2 = i1 + scaffold
    i3 = i2 + len(rtt)
    i4 = i3 + len(pbs)
    return {"spacer": (0, i1), "scaffold": (i1, i2), "rtt": (i2, i3), "pbs1": (i3, i4)}
