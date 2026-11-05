"""What a cargo may change, and what changed (ADR 0029).

Only bases that code for protein are ever recoded. Everything else in a cargo is kept EXACTLY: a promoter's
consensus cannot change without breaking it, a stop codon is not a choice, a linker the user typed is theirs.
This module is the one place that decides which bases those are, and the one place that reports what moved.

**The lock model.** A cargo is a sequence plus a set of *locked* bases. A *coding part* is a stretch read in
frame from its own first base; a codon of it is **free** to change only if it is not a stop codon and touches no
locked base (a codon that touches one is kept whole). The recoder (the INSERT lever, ADR 0022) changes free
codons to other synonymous codons and nothing else. Nothing outside the free codons can move.

**Where the locks come from**, strongest first:

1. ``cargo_parts``: the parts the user (or the web builder) assembled, each with a ``locked`` flag. When given,
   nothing is guessed.
2. ``cargo_locked``: stretches the user says must never change, as 1-based inclusive ``[first, last]``
   nucleotide ranges. They apply on top of 1 or 3.
3. Known elements that are not protein, found by their sequence in either orientation wherever they sit (the
   T7 promoter, with its G run, or the C run of its reverse complement). Always locked, even when the user said
   nothing and even if they sit in what looks like an ORF.
4. Otherwise the cargo is read from its first base in frame 0 and only the ORF up to the first stop codon is
   coding. That ORF counts when it is the whole cargo (a clean CDS, or a CDS and one terminal stop) or at
   least :data:`MIN_ORF_CODONS` long, which a wrong reading frame reaches by chance rarely. Everything after the
   stop is locked.

This is deliberately conservative: an unannotated cargo gets less recoding, never more.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from pegrna_pipeline.hpc.insert_elements import T7_MIN_PROMOTER

__all__ = [
    "MIN_ORF_CODONS",
    "CODE",
    "CargoPart",
    "CargoSpec",
    "BaseChange",
    "ChangeReport",
    "build_cargo_spec",
    "free_codons",
    "cargo_variants",
    "diff_cargo",
    "assert_locks_hold",
    "arm_changes",
    "build_change_report",
    "library_part_is_coding",
]

#: A prefix ORF that is not the whole cargo must be at least this long to be taken for protein. In a wrong
#: reading frame a stop appears about every 21 codons, so 30 or more without one is unlikely by chance.
MIN_ORF_CODONS = 30

_BASES = "TCAG"
_AAS = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"  # standard code, TCAG order
)
CODE: Dict[str, str] = {
    a + b + c: _AAS[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
_ACGT = frozenset("ACGT")
_COMP = str.maketrans("ACGT", "TGCA")

#: Library categories whose members are not protein even when they carry an amino-acid string.
_NON_CODING_CATEGORIES = frozenset({"regulatory", "amplification", "stop"})


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def library_part_is_coding(*, aa: str, category: str = "") -> bool:
    """Is a motif-library entry protein? Yes when it has an amino-acid sequence that is not just stops and is
    not a regulatory, promoter or stop entry. The web builder and the agent use the same rule."""
    return bool(aa) and set(aa) != {"*"} and category not in _NON_CODING_CATEGORIES


@dataclass(frozen=True)
class CargoPart:
    name: str
    start: int  # 0-based, inclusive
    end: int  # 0-based, exclusive
    locked: bool
    why: str  # what the report says about it


@dataclass(frozen=True)
class CargoSpec:
    sequence: str
    parts: Tuple[CargoPart, ...]
    #: Stretches the user excluded, 0-based half-open, merged and sorted. They apply on top of ``parts``.
    excluded: Tuple[Tuple[int, int], ...] = ()
    notes: Tuple[str, ...] = ()

    def locked_mask(self) -> List[bool]:
        mask = [False] * len(self.sequence)
        for p in self.parts:
            if p.locked:
                for i in range(p.start, p.end):
                    mask[i] = True
        for a, b in self.excluded:
            for i in range(a, b):
                mask[i] = True
        return mask

    def part_at(self, i: int) -> CargoPart:
        for p in self.parts:
            if p.start <= i < p.end:
                return p
        raise IndexError(i)


# ---- finding known non-protein elements -----------------------------------------------------------------


def _motif_hits(seq: str) -> List[Tuple[str, int, int]]:
    """``(name, start, end)`` of every known non-protein element, both orientations, sorted, not overlapping.
    A forward hit keeps the run of G's after the core (the transcribed start); a reverse hit keeps the run of
    C's before it (the complement of that run)."""
    core, rc_core = T7_MIN_PROMOTER, _rc(T7_MIN_PROMOTER)
    hits: List[Tuple[str, int, int]] = []
    i = seq.find(core)
    while i >= 0:
        j = i + len(core)
        while j < len(seq) and seq[j] == "G":
            j += 1
        hits.append(("T7 promoter", i, j))
        i = seq.find(core, j)
    i = seq.find(rc_core)
    while i >= 0:
        k = i
        while k > 0 and seq[k - 1] == "C":
            k -= 1
        hits.append(("T7 promoter (reverse)", k, i + len(rc_core)))
        i = seq.find(rc_core, i + len(rc_core))
    hits.sort(key=lambda h: h[1])
    out: List[Tuple[str, int, int]] = []
    for h in hits:
        if out and h[1] < out[-1][2]:
            continue
        out.append(h)
    return out


_MOTIF_WHY = "a known promoter: its sequence is its function, so it is never changed"


# ---- building a spec -------------------------------------------------------------------------------------


def _clean_acgt(label: str, seq: str) -> str:
    s = (seq or "").upper()
    if not s or set(s) - _ACGT:
        raise ValueError(f"{label} must be non-empty ACGT; got {seq!r}")
    return s


def _coding_problem(dna: str) -> Optional[str]:
    """Why ``dna`` cannot be recoded as one in-frame CDS, or ``None`` when it can."""
    if len(dna) % 3:
        return f"{len(dna)} nt is not a multiple of 3"
    aa = [CODE[dna[i : i + 3]] for i in range(0, len(dna), 3)]
    if "*" in aa[:-1]:
        return "it has a stop codon before its last codon"
    return None


def _auto_parts(seq: str) -> Tuple[List[CargoPart], List[str]]:
    notes: List[str] = []
    hits = _motif_hits(seq)
    limit = hits[0][1] if hits else len(seq)  # a coding stretch cannot run into a known element
    whole = limit // 3  # whole codons available before the first known element
    stop = next((k for k in range(whole) if CODE[seq[3 * k : 3 * k + 3]] == "*"), None)
    end = 3 * (stop if stop is not None else whole)
    whole_cargo = not hits and len(seq) % 3 == 0 and (stop is None or stop == len(seq) // 3 - 1)
    orf_codons = end // 3
    coding = orf_codons >= 1 and (whole_cargo or orf_codons >= MIN_ORF_CODONS)
    parts: List[CargoPart] = []
    cursor = 0
    if coding:
        parts.append(
            CargoPart(
                "protein-coding ORF",
                0,
                end,
                False,
                f"auto-detected: frame 0 from the first base to the first stop codon ({orf_codons} codons)",
            )
        )
        cursor = end
    else:
        notes.append(
            "no stretch of the cargo is known to be protein: the first ORF is "
            f"{orf_codons} codon(s) and not the whole cargo, under the {MIN_ORF_CODONS} that would "
            "identify it. Nothing is recoded; mark the coding parts to allow it"
        )

    def gap(a: int, b: int, name: str, why: str) -> None:
        if b > a:
            parts.append(CargoPart(name, a, b, True, why))

    if coding and stop is not None and end == 3 * stop and cursor == end:
        gap(cursor, cursor + 3, "stop codon", "a stop codon is never changed")
        cursor += 3
    rest_why = "outside the protein-coding ORF: kept exactly"
    for name, a, b in hits:  # a stop codon found before the first element cannot reach into it
        gap(cursor, a, "sequence outside the ORF", rest_why)
        parts.append(CargoPart(name, a, b, True, _MOTIF_WHY))
        cursor = b
    gap(cursor, len(seq), "sequence outside the ORF", rest_why)
    return parts, notes


def _parts_from_user(
    seq: str, given: Sequence[Mapping[str, Any]]
) -> Tuple[List[CargoPart], List[str]]:
    parts: List[CargoPart] = []
    notes: List[str] = []
    cursor = 0
    joined: List[str] = []
    for n, item in enumerate(given):
        dna = _clean_acgt(f"cargo_parts[{n}]", item.get("dna") or item.get("sequence") or "")
        joined.append(dna)
        name = str(item.get("name") or f"part {n + 1}")
        end = cursor + len(dna)
        locked = bool(item.get("locked", True))  # unmarked means untouched
        why = (
            "kept exactly as you marked it"
            if locked
            else "marked protein-coding: its codons may change"
        )
        if not locked:
            problem = _coding_problem(dna)
            if problem:
                locked = True
                why = f"marked protein-coding, but {problem}, so it is kept exactly"
                notes.append(f"{name}: {why}")
        parts.append(CargoPart(name, cursor, end, locked, why))
        cursor = end
    if "".join(joined) != seq:
        raise ValueError(
            "cargo_parts do not add up to cargo: their DNA, joined in order, must equal it exactly"
        )
    coding = [p for p in parts if not p.locked]
    if len(coding) > 1:
        first_frame = coding[0].start % 3
        for p in coding[1:]:
            if p.start % 3 != first_frame:
                notes.append(
                    f"{p.name} is read in a different frame from {coding[0].name} (the parts before it "
                    f"are {p.start - coding[0].start} nt apart, not a multiple of 3): intended?"
                )
    return parts, notes


def _excluded_from_user(
    n: int, ranges: Optional[Sequence[Sequence[int]]]
) -> Tuple[Tuple[int, int], ...]:
    """1-based inclusive ``[first, last]`` ranges -> merged 0-based half-open spans."""
    spans: List[Tuple[int, int]] = []
    for r in ranges or ():
        try:
            first, last = int(r[0]), int(r[1])
        except (TypeError, ValueError, IndexError):
            raise ValueError(
                f"cargo_locked entries are [first, last] nucleotide numbers; got {r!r}"
            ) from None
        if not 1 <= first <= last <= n:
            raise ValueError(f"cargo_locked {first}-{last} is outside the cargo (1-{n})")
        spans.append((first - 1, last))
    spans.sort()
    merged: List[Tuple[int, int]] = []
    for a, b in spans:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return tuple(merged)


def build_cargo_spec(
    cargo: str,
    *,
    parts: Optional[Sequence[Mapping[str, Any]]] = None,
    locked: Optional[Sequence[Sequence[int]]] = None,
) -> CargoSpec:
    """The lock model for ``cargo``: from ``parts`` when given, else detected (module docstring)."""
    seq = _clean_acgt("cargo", cargo)
    if parts:
        built, notes = _parts_from_user(seq, parts)
    else:
        built, notes = _auto_parts(seq)
    if not built:
        built = [CargoPart("cargo", 0, len(seq), True, "kept exactly")]
    return CargoSpec(
        sequence=seq,
        parts=tuple(built),
        excluded=_excluded_from_user(len(seq), locked),
        notes=tuple(notes),
    )


# ---- what may change ---------------------------------------------------------------------------------------


def free_codons(spec: CargoSpec) -> List[int]:
    """Start index of every codon the recoder may change: in a coding part, not a stop, touching no locked base."""
    mask = spec.locked_mask()
    out: List[int] = []
    for p in spec.parts:
        if p.locked:
            continue
        for s in range(p.start, p.end - 2, 3):  # whole codons only
            if CODE[spec.sequence[s : s + 3]] == "*" or any(mask[s : s + 3]):
                continue
            out.append(s)
    return out


def cargo_variants(spec: CargoSpec, max_variants: int) -> List[str]:
    """The cargo as given, then up to ``max_variants - 1`` CAI-ranked synonymous variants that change only
    free codons. Every locked base, every stop and every amino acid is identical in every variant.
    """
    out = [spec.sequence]
    free = free_codons(spec)
    if max_variants <= 1 or not free:
        return out
    from pegrna_pipeline.hpc.codon_expander import rank_insert_variants_by_cai  # noqa: PLC0415

    aa = "".join(CODE[spec.sequence[s : s + 3]] for s in free)
    seen = {spec.sequence}
    for cv in rank_insert_variants_by_cai(aa, max_variants=max_variants * 4):
        dna = list(spec.sequence)
        for k, s in enumerate(free):
            dna[s : s + 3] = cv.dna[3 * k : 3 * k + 3]
        new = "".join(dna)
        if new in seen:
            continue
        assert_locks_hold(spec, new)
        seen.add(new)
        out.append(new)
        if len(out) >= max_variants:
            break
    return out


def assert_locks_hold(spec: CargoSpec, new: str) -> None:
    """Raise if ``new`` changes a locked base, a stop codon or an amino acid. A recoder bug must not ship."""
    if len(new) != len(spec.sequence):
        raise AssertionError("a recoded cargo must keep its length")
    mask = spec.locked_mask()
    for i, (a, b) in enumerate(zip(spec.sequence, new)):
        if a != b and mask[i]:
            raise AssertionError(
                f"locked base {i + 1} ({spec.part_at(i).name}) changed: {a} -> {b}"
            )
    free = set(free_codons(spec))
    for i, (a, b) in enumerate(zip(spec.sequence, new)):
        if a == b:
            continue
        p = spec.part_at(i)
        s = p.start + 3 * ((i - p.start) // 3)
        if s not in free:
            raise AssertionError(f"base {i + 1} is not in a codon the recoder may change")
        if CODE[spec.sequence[s : s + 3]] != CODE[new[s : s + 3]]:
            raise AssertionError(f"codon at {s + 1} changed amino acid")


# ---- what changed --------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class BaseChange:
    scope: str  # "cargo" | "genomic arm"
    part: str  # the cargo part, or "replacement arm"
    position: int  # 1-based: in the cargo, or in the design window
    codon: Optional[int]  # 1-based codon number within the part (cargo) or the window's CDS (arm)
    old_codon: str
    new_codon: str
    aa: str
    old: str
    new: str

    def row(self) -> Dict[str, Any]:
        return {
            "scope": self.scope,
            "part": self.part,
            "position": self.position,
            "codon": self.codon,
            "old_codon": self.old_codon,
            "new_codon": self.new_codon,
            "amino_acid": self.aa,
            "old_base": self.old,
            "new_base": self.new,
        }

    def line(self) -> str:
        where = f"{self.part}, nt {self.position}"
        if self.codon:
            return (
                f"{where} (codon {self.codon}, {self.aa}): {self.old_codon} -> {self.new_codon}"
                f"  [{self.old}->{self.new}]"
            )
        return f"{where}: {self.old}->{self.new}"


def diff_cargo(spec: CargoSpec, new: str) -> Tuple[BaseChange, ...]:
    """Every base of ``new`` that differs from the cargo as given, with its part and codon."""
    out: List[BaseChange] = []
    for i, (a, b) in enumerate(zip(spec.sequence, new)):
        if a == b:
            continue
        p = spec.part_at(i)
        s = p.start + 3 * ((i - p.start) // 3)
        old_c, new_c = spec.sequence[s : s + 3], new[s : s + 3]
        out.append(
            BaseChange(
                "cargo",
                p.name,
                i + 1,
                (s - p.start) // 3 + 1,
                old_c,
                new_c,
                CODE.get(old_c, "?"),
                a,
                b,
            )
        )
    return tuple(out)


def arm_changes(
    genome: str, lo: int, hi: int, arm_sense: str, cds_mask: Sequence[bool]
) -> Tuple[BaseChange, ...]:
    """Every base of the design window that the replacement arm changes: window position (1-based), the codon
    it sits in and the amino acid, which the recode keeps. ``genome[lo:hi]`` is what the arm replaces.
    """
    original = genome[lo:hi]
    if len(arm_sense) != len(original):
        raise ValueError("the replacement arm must be as long as the span it replaces")
    window = genome[:lo] + arm_sense + genome[hi:]
    mask = list(cds_mask)
    try:
        cds_start = mask.index(True)
    except ValueError:
        cds_start = None
    out: List[BaseChange] = []
    for i, (a, b) in enumerate(zip(original, arm_sense)):
        if a == b:
            continue
        pos = lo + i
        if cds_start is not None and pos >= cds_start and mask[pos]:
            cs = pos - ((pos - cds_start) % 3)
            old_c, new_c = genome[cs : cs + 3], window[cs : cs + 3]
            out.append(
                BaseChange(
                    "genomic arm",
                    "replacement arm",
                    pos + 1,
                    (cs - cds_start) // 3 + 1,
                    old_c,
                    new_c,
                    CODE.get(old_c, "?"),
                    a,
                    b,
                )
            )
        else:
            out.append(
                BaseChange("genomic arm", "replacement arm", pos + 1, None, "", "", "", a, b)
            )
    return tuple(out)


@dataclass(frozen=True)
class ChangeReport:
    """Everything the pipeline did to what you specified, exactly. Nothing not listed here changed."""

    cargo_changes: Tuple[BaseChange, ...]
    arm_changes: Tuple[BaseChange, ...]
    kept_exactly: Tuple[str, ...]
    cargo_variant: int  # 0 = the cargo as specified
    cargo_length: int
    notes: Tuple[str, ...] = ()

    @property
    def cargo_unchanged(self) -> bool:
        return not self.cargo_changes

    def rows(self) -> List[Dict[str, Any]]:
        return [c.row() for c in (*self.cargo_changes, *self.arm_changes)]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "cargo_unchanged": self.cargo_unchanged,
            "cargo_variant": self.cargo_variant,
            "cargo_length": self.cargo_length,
            "cargo_changes": [c.row() for c in self.cargo_changes],
            "arm_changes": [c.row() for c in self.arm_changes],
            "kept_exactly": list(self.kept_exactly),
            "notes": list(self.notes),
        }

    def lines(self, limit: int = 60) -> List[str]:
        out: List[str] = []
        if self.cargo_unchanged:
            out.append(f"Cargo: used exactly as you gave it ({self.cargo_length} nt, 0 changed).")
        else:
            n = len(self.cargo_changes)
            parts = sorted({c.part for c in self.cargo_changes})
            out.append(
                f"Cargo: RECODED, not as you gave it: {n} nucleotide(s) changed in {', '.join(parts)} "
                f"(variant {self.cargo_variant}; every amino acid is identical)."
            )
            out.extend("  " + c.line() for c in self.cargo_changes[:limit])
            if n > limit:
                out.append(f"  ... {n - limit} more (all in recoded_tjpe_changes.csv)")
        for line in self.kept_exactly:
            out.append(f"  kept exactly: {line}")
        for note in self.notes:
            out.append(f"  note: {note}")
        if self.arm_changes:
            out.append(
                f"Genomic arm: {len(self.arm_changes)} nucleotide(s) of the design window changed by the "
                "silent recode (amino acids identical):"
            )
            out.extend("  " + c.line() for c in self.arm_changes[:limit])
            if len(self.arm_changes) > limit:
                out.append(
                    f"  ... {len(self.arm_changes) - limit} more (all in recoded_tjpe_changes.csv)"
                )
        else:
            out.append("Genomic arm: copied exactly (no nucleotide changed).")
        out.append("Everything not listed above is exactly as you specified.")
        return out


def _spans_text(spec: CargoSpec) -> List[str]:
    lines: List[str] = []
    for p in spec.parts:
        if p.locked:
            lines.append(f"{p.name} (nt {p.start + 1}-{p.end}, {p.end - p.start} nt): {p.why}")
    for a, b in spec.excluded:
        lines.append(f"nt {a + 1}-{b} ({b - a} nt): excluded by you")
    return lines


def build_change_report(
    spec: CargoSpec, used: str, variant: int, arm: Tuple[BaseChange, ...] = ()
) -> ChangeReport:
    """The report for the cargo actually used (``variant`` 0 = as given) and the arm's changes."""
    free = len(free_codons(spec))
    notes = list(spec.notes)
    if not free:
        notes.append(
            "no codon of the cargo may be recoded, so the recoding lever had nothing to change"
        )
    return ChangeReport(
        cargo_changes=diff_cargo(spec, used),
        arm_changes=arm,
        kept_exactly=tuple(_spans_text(spec)),
        cargo_variant=variant,
        cargo_length=len(spec.sequence),
        notes=tuple(notes),
    )
