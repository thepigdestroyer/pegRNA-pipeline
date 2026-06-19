"""Het-SNP advisory + allele-specific pairing + allele-hunt nick movement.

Decisions 12 + 13 + 15. All advisory — never gates. A het SNP overlapping a
DISCRIMINATING position lets a pegRNA-ngRNA pairing edit only one allele
(homogeneous all-het population, often preferable). We classify the SNP's
position (pam/seed/pbs1/distal), estimate discrimination strength, and ALWAYS
surface >=1 allele-specific pairing when available.

Both pegRNA-side and ngRNA-side are valid. ngRNA-side: allele B sees only an SSB
(the ngRNA doesn't cut it); a lone nick is faithfully repaired and MMR reverts
the heteroduplex toward WT, and the large cargo cannot install without the 2nd
nick, so insertion is cleanly allele-specific. Residual = a low-frequency PE2-like
byproduct on allele B (raised under PE4/PE5 MMR-inhibition) — a REPORT NOTE, not
a scar warning.

Specificity cannot be manufactured (the <=2 silent PAM edits write to both
alleles); it must come from a pre-existing PHASED het, so it is capped by natural
heterozygosity (~5-10 % of targets) — compute the exact count from the diploid
genome and report it.

``allele_hunt`` nick movement (decision 15): only when hunting, the pegRNA NICK
may shift (the insert site S NEVER moves) to capture a discriminating het. The
re-templated S<->nick span partitions into a CDS ``recode_span`` (synonymous
recoding allowed) and a non-coding ``locked_span`` (byte-exact to the targeted
allele; microhomology there is an advisory flag).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

__all__ = [
    "classify_snp_position",
    "discrimination_strength",
    "find_allele_specific",
    "enumerate_nick_movements",
    "AlleleSpecificOption",
    "NickMovement",
]

_SEED_NT = 10  # PAM-proximal window length that discriminates


def classify_snp_position(*, element: str, offset: int, length: int) -> str:
    """Classify a SNP's position within a functional element."""
    if element.endswith("pam") or element == "pam":
        return "pam"
    if element == "pbs1":
        return "pbs1"
    if element == "pbs2":
        return "pbs2"
    if element == "spacer":
        # seed = PAM-proximal _SEED_NT (3' end of a Fw spacer)
        return "seed" if offset >= length - _SEED_NT else "distal"
    return "distal"


def discrimination_strength(position_class: str, side: str) -> float:
    """Predicted allele-discrimination strength. PAM-destroying > seed > pbs1 >
    distal(0). ngRNA-side is valid but modestly down-weighted (byproduct note)."""
    base = {"pam": 1.0, "seed": 0.7, "pbs1": 0.4, "pbs2": 0.3, "distal": 0.0}.get(
        position_class, 0.0
    )
    return base * (0.85 if side == "ngrna" else 1.0)


@dataclass(frozen=True)
class AlleleSpecificOption:
    allele_specific: bool
    side: str            # "pegrna" | "ngrna"
    position_class: str
    strength: float
    haplotype: int
    note: str


def find_allele_specific(het_snps: Sequence[dict]) -> List[AlleleSpecificOption]:
    """Return allele-specific pairing options (possibly empty) from PHASED het SNPs."""
    opts: List[AlleleSpecificOption] = []
    for snp in het_snps:
        if not snp.get("phased", False):
            continue
        element = str(snp["element"])
        side = "ngrna" if element.startswith("ngrna") else "pegrna"
        pos = classify_snp_position(
            element=element, offset=int(snp["offset"]), length=int(snp.get("length", 20))
        )
        strength = discrimination_strength(pos, side)
        if strength <= 0.0:
            continue
        if side == "ngrna":
            note = (
                "ngRNA-side: allele B sees only an SSB (nick) -> faithfully repaired, "
                "MMR reverts toward WT, cargo cannot install without the 2nd nick; "
                "residual low-frequency PE2 byproduct on allele B (raised under PE4/PE5)"
            )
        else:
            note = "pegRNA-side: only allele A is nicked -> clean single-allele edit"
        opts.append(
            AlleleSpecificOption(
                allele_specific=True, side=side, position_class=pos,
                strength=strength, haplotype=int(snp["haplotype"]), note=note,
            )
        )
    opts.sort(key=lambda o: -o.strength)
    return opts


@dataclass(frozen=True)
class NickMovement:
    shift: int                                  # signed nt shift of the nick (S fixed)
    shift_abs: int
    direction: str                              # "5prime" | "3prime"
    captured_snp: Optional[dict]
    position_class: Optional[str]
    strength: float
    recode_span: Optional[Tuple[int, int]]      # CDS portion of S<->nick (recodable)
    locked_span: Optional[Tuple[int, int]]      # non-coding portion (byte-exact)
    microhomology_advisory: bool                # True if a non-coding locked span exists


def enumerate_nick_movements(
    *, s_site: int, genome: str, cds_mask: Sequence[bool],
    max_shift: int, het_snps: Sequence[dict],
) -> List[NickMovement]:
    """Allele-hunt only: yield candidate pegRNA-nick shifts (S never moves).

    For each shift the re-templated span between S and the moved nick is
    partitioned via ``cds_mask`` into a CDS ``recode_span`` (recodable) and a
    non-coding ``locked_span`` (byte-exact to the targeted allele). A non-coding
    locked span sets ``microhomology_advisory`` (it cannot be recoded away).
    Candidates are ranked by smallest absolute shift first (minimise locked span
    + RTT length). Distal-only SNPs do not yield discrimination.
    """
    snp_by_pos = {int(s["pos"]): s for s in het_snps if "pos" in s}
    out: List[NickMovement] = []

    for shift in range(1, max_shift + 1):
        for direction, nick in (("5prime", s_site - shift), ("3prime", s_site + shift)):
            if nick < 0 or nick >= len(genome):
                continue
            lo, hi = sorted((s_site, nick))
            # Partition [lo, hi) into CDS vs non-coding using the mask.
            cds_positions = [i for i in range(lo, hi) if i < len(cds_mask) and cds_mask[i]]
            nc_positions = [i for i in range(lo, hi) if i < len(cds_mask) and not cds_mask[i]]
            recode_span = (min(cds_positions), max(cds_positions) + 1) if cds_positions else None
            locked_span = (min(nc_positions), max(nc_positions) + 1) if nc_positions else None
            # Did this shift bring a discriminating SNP into the moved region?
            captured = None
            pos_class = None
            strength = 0.0
            for p, snp in snp_by_pos.items():
                if lo <= p < hi or p == nick:
                    pc = classify_snp_position(
                        element=str(snp.get("element", "spacer")),
                        offset=int(snp.get("offset", 0)),
                        length=int(snp.get("length", 20)),
                    )
                    st = discrimination_strength(pc, "pegrna")
                    if st > strength:
                        captured, pos_class, strength = snp, pc, st
            out.append(
                NickMovement(
                    shift=(-shift if direction == "5prime" else shift),
                    shift_abs=shift,
                    direction=direction,
                    captured_snp=captured,
                    position_class=pos_class,
                    strength=strength,
                    recode_span=recode_span,
                    locked_span=locked_span,
                    microhomology_advisory=locked_span is not None,
                )
            )
    out.sort(key=lambda m: (m.shift_abs, m.direction))
    return out
