"""Shared helpers for the nick-convention and Prime Assembly tests.

``make_locus`` builds a small, deterministic locus in which every guide is one we planted: the background
contains no GG and no CC, so it offers no accidental PAM on either strand, and the coding part has no
in-frame stop codon. That makes "which guides exist here" a fact of the fixture and not an accident of a
random sequence. The cut of each planted guide is given in the convention of ``docs/STRUCTURE.md``: the
index of the first base to the right of the cut.

The real CTCF window is used elsewhere, behind ``importorskip("dnachisel")``.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from Bio.Seq import Seq

from pegrna_pipeline.hpc.complex_dock import DockFake
from pegrna_pipeline.hpc.cryptic_splice import MaxEntScanFake
from pegrna_pipeline.hpc.full_structure import FullStructureFake
from pegrna_pipeline.hpc.splice_screen import SpliceFake
from pegrna_pipeline.hpc.structure_importance import ResidueMetrics

STOPS = {"TAA", "TAG", "TGA"}
NON_STOP_CODONS = [
    a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS
]


def rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


@dataclass(frozen=True)
class Locus:
    genome: str
    site_s: int
    cds_mask: List[bool]
    guides: Tuple[Tuple[str, int, str], ...]  # (strand, cut, protospacer), as planted


def _ok_pair(prev: str, base: str, run: List[str]) -> bool:
    if prev == base and base in "GC":  # no GG, no CC
        return False
    return not (len(run) >= 3 and run[-1] == run[-2] == run[-3] == base)  # no run of 4


def _background(rng: random.Random, n: int, coding_end: int) -> List[str]:
    """Bases with no GG, no CC and no run of four, and no in-frame stop codon before ``coding_end``."""
    out: List[str] = []
    while len(out) < n:
        base = rng.choice("ACGT")
        if out and not _ok_pair(out[-1], base, out):
            continue
        if len(out) % 3 == 2 and len(out) < coding_end and "".join(out[-2:]) + base in STOPS:
            continue
        out.append(base)
    return out


def _protospacer(rng: random.Random) -> str:
    while True:
        s = "".join(rng.choice("ACGT") for _ in range(20))
        gc = (s.count("G") + s.count("C")) / 20
        if "GG" in s or "CC" in s or any(c * 4 in s for c in "ACGT") or not 0.40 <= gc <= 0.60:
            continue
        return s


def _isolate(trial: List[str], idx: int) -> None:
    """Set the base beside a planted block to an A or T that does not extend a run of identical bases, and
    that cannot form a PAM (no G or C) with the block's edge."""
    for base in "AT":
        trial[idx] = base
        if not any(c * 4 in "".join(trial[max(0, idx - 3) : idx + 4]) for c in "ACGT"):
            return


def _valid(seq: Sequence[str], coding_end: int, n_top: int, n_bottom: int) -> bool:
    s = "".join(seq)
    if s.count("GG") != n_top or s.count("CC") != n_bottom:
        return False
    if any(c * 4 in s for c in "ACGT"):
        return False
    return not any(s[i : i + 3] in STOPS for i in range(0, coding_end - 2, 3))


def make_locus(
    plants: Sequence[Tuple[str, int]],
    *,
    n: int = 300,
    site_s: int = 150,
    seed: int = 7,
    coding_end: Optional[int] = None,
) -> Locus:
    """A locus of ``n`` bases with S at ``site_s`` and one guide planted per ``(strand, cut)``.

    ``strand`` is "top" or "bottom". A top guide's protospacer is ``[cut-17, cut+3)`` and its PAM
    (``TGG``) ``[cut+3, cut+6)``; a bottom guide shows ``CCA`` at ``[cut-6, cut-3)`` and its protospacer
    reverse-complemented at ``[cut-3, cut+17)``. Coding is everything before ``coding_end``, which is S
    unless given (a PRIDICT-path fixture that walks the design through an arm 3' of S sets it to ``n``).
    """
    coding_end = site_s if coding_end is None else coding_end
    rng = random.Random(seed)
    seq = _background(rng, n, coding_end)
    placed: List[Tuple[str, int, str]] = []
    for strand, cut in plants:
        for _ in range(400):
            proto = _protospacer(rng)
            lo, block = (
                (cut - 17, proto + "TGG") if strand == "top" else (cut - 6, "CCA" + rc(proto))
            )
            trial = list(seq)
            trial[lo : lo + len(block)] = list(block)
            if lo > 0:
                _isolate(trial, lo - 1)
            if lo + len(block) < n:
                _isolate(trial, lo + len(block))
            done = placed + [(strand, cut, proto)]
            n_top = sum(1 for s, _, _ in done if s == "top")
            if _valid(trial, coding_end, n_top, len(done) - n_top):
                seq, placed = trial, done
                break
        else:  # pragma: no cover - a fixture that cannot be built is a test bug
            raise RuntimeError(f"could not plant a {strand} guide at {cut}")
    return Locus(
        genome="".join(seq),
        site_s=site_s,
        cds_mask=[i < coding_end for i in range(n)],
        guides=tuple(placed),
    )


def make_cargo(length: int = 120, seed: int = 11) -> str:
    """An in-frame cargo with no stop codon and no run of four identical bases."""
    rng = random.Random(seed)
    while True:
        cargo = "".join(rng.choice(NON_STOP_CODONS) for _ in range(length // 3))
        if not any(c * 4 in cargo for c in "ACGT"):
            return cargo


# ---- what real PRIDICT2 reports for a nicking guide ----------------------------------------------


def pridict_distance(label: str, true_rel: int, insert_len: int = 24) -> int:
    """The ``Nicking-Position-to-edit`` real PRIDICT2 reported for a guide whose cut is ``true_rel`` from S.

    Measured on the CTCF window with a 24 bp insert (``scripts/analysis/probe_pridict_nick_convention.py``,
    39 guides, all consistent): a forward-strand guide upstream of S reads ``true - 1``, downstream
    ``true + insert_len - 1``; a reverse-strand guide downstream reads ``-true``, upstream
    ``-true + insert_len``. It is not the position of the cut, which is why the designer ignores it.
    """
    if label == "Fw":
        return true_rel - 1 if true_rel < 0 else true_rel + insert_len - 1
    return -true_rel if true_rel > 0 else -true_rel + insert_len


def nicking_row(spacer: str, distance: int, label: str = "Fw", score: float = 0.7) -> dict:
    """A PRIDICT ``*_nicking_guides.csv`` row as ``csv.DictReader`` returns it."""
    return {
        "Nicking-Protospacer": spacer,
        "Nicking-Position-to-edit": str(distance),
        "Nicking-PAMdisrupt": "Nicking_PAM_disrupt",
        "Target_Strand": label,
        "DeepCas9score": str(score),
    }


def guide_row(locus: "Locus", index: int, *, insert_len: int = 24, score: float = 0.7) -> dict:
    """The row real PRIDICT would emit for the ``index``-th guide planted in ``locus``."""
    strand, cut, spacer = locus.guides[index]
    label = "Fw" if strand == "top" else "Rv"
    return nicking_row(
        spacer, pridict_distance(label, cut - locus.site_s, insert_len), label, score
    )


# ---- fake gate scorers -----------------------------------------------------------------------


class AllPassImportance:
    def metrics(self, i):  # noqa: ARG002
        return ResidueMetrics(
            plddt=50.0,
            rsa=0.5,
            contact_number=4,
            is_catalytic=False,
            is_disulfide=False,
            is_positive_phi_gly=False,
            is_structural_pro=False,
        )


def clean_offtarget(protospacers) -> Dict[str, Tuple[float, int, bool]]:
    return {p: (0.8, 0, False) for p in protospacers}


def fake_scorers(**over):
    base = dict(
        splice_scorer=SpliceFake(delta=0.0),
        full_structure_scorer=FullStructureFake(pairs={}),
        importance_scorer=AllPassImportance(),
        dock_scorer=DockFake(apo_iptm=0.85, candidate_iptm=0.82),
        offtarget_scorer=clean_offtarget,
        maxent_scorer=MaxEntScanFake(),
    )
    base.update(over)
    return base
