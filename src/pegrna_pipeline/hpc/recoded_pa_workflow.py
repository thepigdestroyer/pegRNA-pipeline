"""Prime Assembly on the canonical workflow, with the donor restoring the DNA between the nicks (spec 4.2).

Two extended guides on opposite strands write two flaps, and the ends of a separate linear donor match them
(Liu et al., Nature 2026). Here the donor also RESTORES the DNA between the nicks, as the Replacement Arm does in
tj-PE (decision A+): ``donor = arm5 + cargo + arm3``, the product is ``genome[:n1] + donor + genome[n2:]``, and the net
deletion is zero.

The walk is a ranked search over PAM-in guide pairs:

1. enumerate guides, drop those that fail the sequence rules or the off-target screen (one batched call);
2. pair them, and for each pair build the arms, break the PAMs where a synonymous change exists, assemble the donor
   and the product, assert INV-0 on each arm, choose each side's extension and check it with ``simulate_side``;
3. rank the pairs by deterministic risks (:func:`rank_key`), best first;
4. for each pair, in order: the cheap per-pegRNA fold, the product-level splice gates, then the expensive docking,
   every per-pegRNA verdict memoized by the pegRNA's sequence so a pegRNA that appears in many pairs is folded and
   docked once; the first pair that passes everything wins, after ``simulate_pa`` has confirmed the pair.

The cargo is never changed here (the INSERT lever of tj-PE does not apply to a donor); the lock inputs are accepted,
validated and reported. Only the genomic arms are recoded, and only where the bases code for protein.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from pegrna_pipeline.hpc.cargo_spec import (
    BaseChange,
    ChangeReport,
    arm_changes,
    build_cargo_spec,
    build_change_report,
)
from pegrna_pipeline.hpc.complex_dock import DockingServiceError, passes_docking_gate
from pegrna_pipeline.hpc.cryptic_splice import passes_cryptic_splice_gate
from pegrna_pipeline.hpc.full_structure import (
    functional_region_pairing,
    scaffold_extension_pairing,
)
from pegrna_pipeline.hpc.gate_thresholds import THRESHOLDS
from pegrna_pipeline.hpc.guide_sanity import assess_spacer
from pegrna_pipeline.hpc.microhomology import longest_common_run
from pegrna_pipeline.hpc.nick_search import offtarget_rejection
from pegrna_pipeline.hpc.pa_geometry import (
    MIN_NICK_SEPARATION_BP,
    Arm,
    Extension,
    GuidePair,
    PamOutcome,
    assemble_donor,
    assemble_pa_pegrna,
    break_pams,
    build_arm,
    choose_extension,
    find_guides,
    flap_lengths,
    genomic_run_under_flap,
    pa_regions,
    pam_in_pairs,
    product_allele,
)
from pegrna_pipeline.hpc.pa_simulator import PAMismatch, simulate_pa, simulate_side
from pegrna_pipeline.hpc.recoded_tjpe_workflow import (
    NoValidPegRNA,
    _advisory_splice_text,
    _arm_protein_equal,
    _folder_name,
    _frame_at,
    _importance_metrics,
    _insertion_site_footprint,
    _resolve_protein_length,
    resolve_structure_threshold,
)
from pegrna_pipeline.hpc.scarless import assert_inv0
from pegrna_pipeline.hpc.splice_screen import passes_splice_gate
from pegrna_pipeline.hpc.structure_importance import span_blocked

__all__ = [
    "STANDING_CAVEAT",
    "NoValidPADesign",
    "SideCandidate",
    "PairPlan",
    "PairPlanning",
    "PAGuideDesign",
    "RecodedPAResult",
    "rank_key",
    "plan_pairs",
    "run_recoded_pa",
]

STANDING_CAVEAT = (
    "A donor that restores the DNA between the nicks is this project's extension of Prime Assembly (Liu et al., "
    "Nature 2026), which replaced that region with non-genomic sequence. Nothing has validated a donor whose arms "
    "copy the genome under the flaps, so treat this design as untested: check each guide's nick on its own, and "
    "the assembly itself, before relying on it."
)

_PAM_NOTE = {
    "intact-noncoding": "lies in non-coding DNA, which is kept exactly",
    "intact-no-synonymous-change": "lies in coding sequence but no synonymous change breaks it",
}


class NoValidPADesign(NoValidPegRNA):
    """No guide pair survived: the diagnostics say, pair by pair, what rejected it."""


# ---- value types --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class SideCandidate:
    """One pegRNA for one side: a PBS length with the flap length the sequence rules chose for it."""

    extension: Extension
    pegrna: str


@dataclass(frozen=True)
class PairPlan:
    """Everything about a pair that is known before any model runs. ``rank_key`` orders the plans."""

    pair: GuidePair
    arm5: Arm
    arm3: Arm
    pam_outcomes: Tuple[PamOutcome, PamOutcome]
    donor: str
    product: str
    candidates_a: Tuple[SideCandidate, ...]
    candidates_b: Tuple[SideCandidate, ...]
    run_a: int
    run_b: int
    advisories: Tuple[str, ...]
    spacer_advisories: Tuple[Tuple[str, ...], Tuple[str, ...]]
    on_target: Tuple[float, float]
    offtargets: Tuple[int, int]
    rank_key: tuple


@dataclass
class PairPlanning:
    plans: List[PairPlan]  # ranked, best first
    diagnostics: List[str]
    pairs_found: int  # PAM-in pairs left after the guide screens
    pairs_geometric: int = 0  # PAM-in pairs before them


@dataclass
class PAGuideDesign:
    side: str  # "A" (top strand, nick n1) | "B" (bottom strand, nick n2)
    strand: str
    spacer: str
    cut: int
    pam: str
    pam_status: str  # PamOutcome.status
    #: The PAM of the product as this guide reads it: NGG when regenerated, NAG or NGA when only weakened.
    product_pam: str
    flap_len: int
    pbs_len: int
    pbs: str
    rtt: str
    tail: str
    pegRNA_sequence: str
    offtargets: int
    on_target: float
    #: Longest stretch the flap shares with the strand it displaces (a copy of the genome can re-anneal there).
    genomic_run: int


@dataclass
class RecodedPAResult:
    """Outcome of a Prime Assembly design: two pegRNAs, one donor, the predicted product."""

    guide_a: PAGuideDesign
    guide_b: PAGuideDesign
    site_s: int
    span: Tuple[int, int]  # (n1, n2): the nicks; the donor replaces genome[n1:n2)
    cargo: str
    arm5: str
    arm3: str
    donor: str
    product: str
    donor_parts: List[Tuple[str, str]]  # the non-empty parts of arm5 + cargo + arm3, in order
    recode_edits: int
    deletion_size: int  # always 0: INV-0
    candidates_tried: int
    pairs_considered: int
    warnings: List[str] = field(default_factory=list)
    gate_results: List[Tuple[str, str]] = field(default_factory=list)
    #: Genome-window spans (half-open) that were recoded, and that were copied exactly because non-coding.
    recode_spans: List[Tuple[int, int]] = field(default_factory=list)
    locked_spans: List[Tuple[int, int]] = field(default_factory=list)
    #: Every base the donor's arms change from the window, with codon and amino acid (ADR 0029).
    arm_changes: Tuple[BaseChange, ...] = ()
    changes: Optional[ChangeReport] = None
    simulated: bool = False  # simulate_pa confirmed the ordered pegRNAs and donor


# ---- ranking --------------------------------------------------------------------------------------------


def rank_key(
    *,
    advisories: int,
    run_total: int,
    span: int,
    on_target_sum: float,
    a_cut: int,
    b_cut: int,
) -> tuple:
    """Total, deterministic order on pairs, best first (a heuristic: there is no wet-lab data for PA here).

    Fewest applicable advisories; then the smaller combined genomic run under the two flaps; then the smaller span;
    then the HIGHER summed on-target score; then the cut positions, so the order never depends on iteration order.
    It prefers pairs that avoid the two known risks: a flap that is a copy of the genome, and a regenerated PAM.

    The span IS the number of restored base pairs (the donor's arms restore exactly what lies between the nicks),
    so "fewer restored base pairs" and "the smaller span" are one criterion. A shorter span needs fewer recoding
    edits (about 17 for a 91-nt arm), which is why it ranks above the on-target score.
    """
    return (advisories, run_total, span, -round(on_target_sum, 6), a_cut, b_cut)


# ---- argument checks ------------------------------------------------------------------------------------


def _check_args(
    *,
    cargo: str,
    genome: str,
    site_s: int,
    cds_mask: Sequence[bool],
    pbs_lengths: Sequence[int],
    flap_length: int,
    flap_min: int,
    flap_max: int,
    min_span: int,
    max_span: int,
    max_run: int,
    tail: str,
    max_offtargets: int,
    max_pairs: int,
) -> None:
    if not cargo or set(cargo.upper()) - set("ACGT"):
        raise ValueError("cargo must be non-empty DNA (A, C, G, T only)")
    if len(cds_mask) != len(genome):
        raise ValueError(
            f"cds_mask has {len(cds_mask)} entries but the genome window has {len(genome)} bases"
        )
    if not 0 <= site_s <= len(genome):
        raise ValueError(f"site_s must be between 0 and {len(genome)}; got {site_s}")
    if flap_min > flap_max:
        raise ValueError(f"flap_min ({flap_min}) must not exceed flap_max ({flap_max})")
    if not flap_min <= flap_length <= flap_max:
        raise ValueError(
            f"flap_length must be between flap_min and flap_max ({flap_min}-{flap_max}); got {flap_length}"
        )
    if flap_min < 30 or flap_max > 50:
        raise ValueError(
            "flap_min and flap_max must stay within 30-50 nt, the range the paper found robust"
        )
    if not pbs_lengths or any(int(x) != x or not 6 <= x <= 30 for x in pbs_lengths):
        raise ValueError(
            f"pbs_lengths must be a non-empty list of whole numbers from 6 to 30; got {list(pbs_lengths)}"
        )
    if min_span < MIN_NICK_SEPARATION_BP:
        raise ValueError(
            f"min_span must be at least {MIN_NICK_SEPARATION_BP} bp (the project floor); got {min_span}"
        )
    if max_span < min_span:
        raise ValueError(f"max_span ({max_span}) must be at least min_span ({min_span})")
    if int(max_run) != max_run or max_run < 1:
        raise ValueError(f"max_run must be a whole number >= 1; got {max_run!r}")
    if not tail or set(tail.upper()) - {"U"} or not 2 <= len(tail) <= 3:
        raise ValueError(f"tail must be 2-3 uridines (U only); got {tail!r}")
    if int(max_offtargets) != max_offtargets or max_offtargets < 0:
        raise ValueError(f"max_offtargets must be an integer >= 0; got {max_offtargets!r}")
    if int(max_pairs) != max_pairs or max_pairs < 1:
        raise ValueError(f"max_pairs must be a whole number >= 1; got {max_pairs!r}")


def _flap_ceiling(cargo: str, flap_length: int, flap_min: int, flap_max: int) -> int:
    """The longest flap the cargo can support, or a refusal.

    Below ``2 * F`` the two flaps overlap each other, and a flap can reach into the other side's arm. That is
    twinPE, not assembly, and ``design-recoded`` is the tool for a cargo that short.
    """
    if len(cargo) < 2 * flap_length:
        raise NoValidPADesign(
            f"the cargo is {len(cargo)} nt, shorter than two flaps ({2 * flap_length} nt = 2 x {flap_length}): "
            "the two flaps would overlap, which is not Prime Assembly. Use design-recoded (tj-PE) for a cargo "
            "this short, or a longer cargo",
            diagnostics=[f"cargo {len(cargo)} nt < 2 x flap {flap_length} = {2 * flap_length} nt"],
        )
    return min(flap_max, len(cargo) // 2)


# ---- planning ---------------------------------------------------------------------------------------------


def _why_no_extension(
    side: str, genome: str, cut: int, donor: str, pbs_len: int, lengths: List[int]
) -> str:
    from pegrna_pipeline.hpc.pa_geometry import derive_extension

    exts = [derive_extension(side, genome, cut, donor, f, pbs_len) for f in lengths]
    exts = [e for e in exts if e is not None]
    if not exts:
        return "the PBS window leaves the genome window, or the flap is longer than the donor"
    return f"TTTT (a Pol III terminator) in the RTT or PBS at every flap length {min(lengths)}-{max(lengths)}"


def _side_candidates(
    side: str,
    guide,
    *,
    genome: str,
    donor: str,
    pbs_lengths: Sequence[int],
    flap_length: int,
    flap_min: int,
    flap_hi: int,
    tail: str,
    lab: str,
    diagnostics: List[str],
) -> List[SideCandidate]:
    """One checked pegRNA per PBS length: the flap length the sequence rules choose, then ``simulate_side``."""
    out: List[SideCandidate] = []
    seen = set()
    lengths = flap_lengths(flap_length, flap_min, flap_hi)
    for pbs_len in pbs_lengths:
        ext = choose_extension(
            side, genome, guide.cut, donor, pbs_len, default=flap_length, lo=flap_min, hi=flap_hi
        )
        if ext is None:
            continue
        if (ext.flap_len, ext.pbs_len) in seen:
            continue
        seen.add((ext.flap_len, ext.pbs_len))
        peg = assemble_pa_pegrna(spacer=guide.spacer, rtt=ext.rtt, pbs=ext.pbs, tail=tail)
        try:
            simulate_side(
                side=side,
                genome=genome,
                guide=guide,
                donor=donor,
                peg=peg,
                flap_len=ext.flap_len,
                pbs_len=ext.pbs_len,
                tail_len=len(tail),
            )
        except PAMismatch as exc:
            diagnostics.append(f"{lab}: pegRNA {side} (PBS {pbs_len}): PA strand mismatch: {exc}")
            continue
        out.append(SideCandidate(ext, peg))
    if not out and not any(
        "PA strand mismatch" in d and f"pegRNA {side}" in d and lab in d for d in diagnostics
    ):
        diagnostics.append(
            f"{lab}: guide {side} (cut {guide.cut}): no usable extension: "
            + _why_no_extension(side, genome, guide.cut, donor, min(pbs_lengths), lengths)
        )
    return out


def _protein_silent(genome: str, arm: Arm, cds_mask: Sequence[bool]) -> bool:
    if arm.recode_span is None:
        return True
    a, b = arm.recode_span
    frame = _frame_at(cds_mask, a)
    return _arm_protein_equal(genome[a:b], arm.sense[a - arm.lo : b - arm.lo], frame)


def _check_inv0(
    genome: str, cargo: str, site_s: int, arm5: Arm, arm3: Arm, cds_mask: Sequence[bool]
) -> None:
    """INV-0 for each arm, on exactly the shapes tj-PE uses (the other arm put back as the genome has it)."""
    shapes = (
        (arm5, genome[: arm5.lo] + arm5.sense + cargo + genome[site_s:]),  # tj-PE's "nick 5' of S"
        (arm3, genome[:site_s] + cargo + arm3.sense + genome[arm3.hi :]),  # tj-PE's "nick 3' of S"
    )
    for arm, edited in shapes:
        assert_inv0(
            reference=genome,
            edited=edited,
            site_s=site_s,
            cargo=cargo,
            recode_span=arm.recode_span if arm.recode_span is not None else (arm.lo, arm.lo),
            protein_check=lambda _r, _e, _arm=arm: _protein_silent(genome, _arm, cds_mask),
            allowed_silent=None,
            locked_span=arm.locked_span,
            locked_expected=arm.locked_expected,
        )


def _pam_window(guide) -> str:
    lo, hi = guide.pam_span
    return f"nt {lo + 1}-{hi}"


def _advisories(
    pam_outcomes: Sequence[PamOutcome],
    guides: Sequence[Any],
    runs: Sequence[int],
    flaps: Sequence[int],
    spacer_advisories: Sequence[Sequence[str]],
    max_run: int,
) -> List[str]:
    out: List[str] = []
    for side, outcome, guide in zip("AB", pam_outcomes, guides):
        if outcome.status == "weakened":
            out.append(
                f"guide {side}'s PAM ({_pam_window(guide)}) can only be weakened, not broken: it reads "
                f"{outcome.pam} in the product, and SpCas9 still cuts that at a fraction of its NGG rate (about a "
                "fifth for NAG, less for NGA), so the guide may still nick the product"
            )
        note = _PAM_NOTE.get(outcome.status)
        if note:
            out.append(
                f"guide {side}'s PAM ({_pam_window(guide)}) {note}, so the donor restores it and the guide can "
                "cut the product again"
            )
    for side, run, flap in zip("AB", runs, flaps):
        if run > max_run:
            out.append(
                f"flap {side} ({flap} nt) shares a run of {run} identical bases with the strand it displaces "
                f"(more than {max_run}): it can re-anneal where it came from instead of on the donor, and that "
                "stretch is non-coding (or otherwise cannot be recoded)"
            )
    for side, advisories in zip("AB", spacer_advisories):
        out.extend(f"guide {side} spacer advisory: {a}" for a in advisories)
    return out


def plan_pairs(
    *,
    cargo: str,
    genome: str,
    site_s: int,
    cds_mask: Sequence[bool],
    offtarget_scorer: Callable[[List[str]], Dict[str, Tuple[float, int, bool]]],
    pbs_lengths: Sequence[int] = (13, 15, 17),
    flap_length: int = 35,
    flap_min: int = 30,
    flap_max: int = 50,
    min_span: int = 30,
    max_span: int = 100,
    max_run: int = 8,
    tail: str = "UUU",
    max_offtargets: int = 0,
    min_ontarget: float = 0.5,
) -> PairPlanning:
    """Enumerate, screen, pair, build and rank. Everything here is deterministic and runs no model but the
    off-target screen (one batched call). A pair that cannot be built is not an error: it is a diagnostic.
    """
    flap_hi = _flap_ceiling(cargo, flap_length, flap_min, flap_max)
    genome = genome.upper()
    cargo = cargo.upper()
    diagnostics: List[str] = []

    every = find_guides(genome)
    near = [g for g in every if site_s - max_span <= g.cut <= site_s + max_span]
    n_top = sum(1 for g in every if g.strand == "top")
    pairs0 = pam_in_pairs(near, site_s=site_s, min_span=min_span, max_span=max_span)
    if not pairs0:
        diagnostics.append(
            f"0 PAM-in pairs: the window has {n_top} top-strand and {len(every) - n_top} bottom-strand NGG "
            f"guides, but none make a top-strand guide cutting at or 5' of S and a bottom-strand guide cutting "
            f"at or 3' of S that are {min_span}-{max_span} bp apart with their PAMs clear of each other"
        )
        return PairPlanning([], diagnostics, 0, 0)

    members = {}
    for p in pairs0:
        members[(p.a.strand, p.a.cut)] = p.a
        members[(p.b.strand, p.b.cut)] = p.b
    shortest = min(pbs_lengths)
    sane: Dict[Tuple[str, int], Any] = {}
    bad_rules: Dict[str, int] = {}
    spacer_adv: Dict[Tuple[str, int], Tuple[str, ...]] = {}
    for key, g in sorted(members.items()):
        fits = g.cut >= shortest if g.strand == "top" else g.cut + shortest <= len(genome)
        verdict = assess_spacer(g.spacer)
        if not fits:
            bad_rules["the PBS window leaves the genome window"] = (
                bad_rules.get("the PBS window leaves the genome window", 0) + 1
            )
        elif not verdict.ok:
            for failure in verdict.failures:
                bad_rules[failure.split(":")[0]] = bad_rules.get(failure.split(":")[0], 0) + 1
        else:
            sane[key] = g
            spacer_adv[key] = tuple(verdict.advisories)

    scores: Dict[str, Tuple[float, int, bool]] = {}
    usable: Dict[Tuple[str, int], Any] = {}
    screened_out = 0
    if sane:
        spacers = sorted({g.spacer for g in sane.values()})
        scores = offtarget_scorer(spacers)
        for key, g in sane.items():
            score = scores.get(g.spacer, (0.0, 10**9, True))
            if offtarget_rejection(score, max_offtargets=max_offtargets, min_ontarget=min_ontarget):
                screened_out += 1
            else:
                usable[key] = g

    pairs = [
        p for p in pairs0 if (p.a.strand, p.a.cut) in usable and (p.b.strand, p.b.cut) in usable
    ]
    parts = []
    n_bad = len(members) - len(sane)
    if n_bad:
        parts.append(
            f"{n_bad} failed the sequence rules ({'; '.join(sorted(bad_rules)) or 'see above'})"
        )
    if screened_out:
        parts.append(
            f"{screened_out} failed the off-target screen (each needs a unique on-target and at most "
            f"{max_offtargets} off-target site(s) within the gate mismatch tier)"
        )
    diagnostics.append(
        f"{len(pairs0)} PAM-in pair(s) from {len(members)} guide(s)"
        + (f"; {', '.join(parts)}" if parts else "")
        + f"; {len(pairs)} pair(s) remain"
    )
    if not pairs:
        return PairPlanning([], diagnostics, 0, len(pairs0))

    plans: List[PairPlan] = []
    for pair in pairs:
        plan = _plan_one(
            pair,
            cargo=cargo,
            genome=genome,
            site_s=site_s,
            cds_mask=cds_mask,
            pbs_lengths=pbs_lengths,
            flap_length=flap_length,
            flap_min=flap_min,
            flap_hi=flap_hi,
            max_run=max_run,
            tail=tail,
            scores=scores,
            spacer_adv=spacer_adv,
            diagnostics=diagnostics,
        )
        if plan is not None:
            plans.append(plan)
    plans.sort(key=lambda p: p.rank_key)
    return PairPlanning(plans, diagnostics, len(pairs), len(pairs0))


def _plan_one(
    pair: GuidePair,
    *,
    cargo: str,
    genome: str,
    site_s: int,
    cds_mask: Sequence[bool],
    pbs_lengths: Sequence[int],
    flap_length: int,
    flap_min: int,
    flap_hi: int,
    max_run: int,
    tail: str,
    scores: Dict[str, Tuple[float, int, bool]],
    spacer_adv: Dict[Tuple[str, int], Tuple[str, ...]],
    diagnostics: List[str],
) -> Optional[PairPlan]:
    a, b = pair.a, pair.b
    lab = f"pair {a.cut}/{b.cut}"
    arm5 = build_arm(genome, a.cut, site_s, cds_mask, max_run)
    arm3 = build_arm(genome, site_s, b.cut, cds_mask, max_run)
    if arm5 is None or arm3 is None:
        diagnostics.append(
            f"{lab}: arm partition failed (no synonymous CDS solution or non-monotone arm)"
        )
        return None
    arm5, arm3, outcomes = break_pams(genome, arm5, arm3, [("A", a), ("B", b)], cds_mask, max_run)
    donor = assemble_donor(arm5.sense, cargo, arm3.sense)
    product = product_allele(genome, a.cut, b.cut, donor)
    if len(product) != len(genome) + len(cargo):
        diagnostics.append(
            f"{lab}: the product is {len(product)} nt, expected {len(genome) + len(cargo)}"
        )
        return None
    try:
        _check_inv0(genome, cargo, site_s, arm5, arm3, cds_mask)
    except ValueError as exc:
        diagnostics.append(f"{lab}: {exc}")  # assert_inv0's messages already begin "INV-0:"
        return None

    common = dict(
        genome=genome,
        donor=donor,
        pbs_lengths=pbs_lengths,
        flap_length=flap_length,
        flap_min=flap_min,
        flap_hi=flap_hi,
        tail=tail,
        lab=lab,
        diagnostics=diagnostics,
    )
    cands_a = _side_candidates("A", a, **common)
    cands_b = _side_candidates("B", b, **common)
    if not cands_a or not cands_b:
        return None

    fa, fb = cands_a[0].extension.flap_len, cands_b[0].extension.flap_len
    shared = longest_common_run(donor[:fa], donor[-fb:])
    if shared > max_run:
        diagnostics.append(
            f"{lab}: the donor's two ends share a run of {shared} bp (more than {max_run}): the two flaps could "
            "pair with each other instead of with the donor"
        )
        return None

    run_a = genomic_run_under_flap("A", genome, a.cut, b.cut, donor, fa)
    run_b = genomic_run_under_flap("B", genome, a.cut, b.cut, donor, fb)
    key_a, key_b = (a.strand, a.cut), (b.strand, b.cut)
    advisories = _advisories(
        outcomes,
        (a, b),
        (run_a, run_b),
        (fa, fb),
        (spacer_adv.get(key_a, ()), spacer_adv.get(key_b, ())),
        max_run,
    )
    sa, sb = scores.get(a.spacer, (0.0, 0, False)), scores.get(b.spacer, (0.0, 0, False))
    return PairPlan(
        pair=pair,
        arm5=arm5,
        arm3=arm3,
        pam_outcomes=(outcomes[0], outcomes[1]),
        donor=donor,
        product=product,
        candidates_a=tuple(cands_a),
        candidates_b=tuple(cands_b),
        run_a=run_a,
        run_b=run_b,
        advisories=tuple(advisories),
        spacer_advisories=(tuple(spacer_adv.get(key_a, ())), tuple(spacer_adv.get(key_b, ()))),
        on_target=(float(sa[0]), float(sb[0])),
        offtargets=(int(sa[1]), int(sb[1])),
        rank_key=rank_key(
            advisories=len(advisories),
            run_total=run_a + run_b,
            span=pair.span,
            on_target_sum=float(sa[0]) + float(sb[0]),
            a_cut=a.cut,
            b_cut=b.cut,
        ),
    )


# ---- the walk -----------------------------------------------------------------------------------------


def _emitter(on_gate: Optional[Callable[..., None]]) -> Callable[..., None]:
    """The live progress hook, per gate and per pair. It must never break a run, and an old four-argument
    callback (no ``metrics``) still works."""

    def emit(
        cand: int, gate_id: str, info: str, passed: bool, metrics: Optional[dict] = None
    ) -> None:
        if not on_gate:
            return
        try:
            on_gate(cand, gate_id, info, passed, metrics)
        except TypeError:
            try:
                on_gate(cand, gate_id, info, passed)
            except Exception:  # noqa: BLE001
                pass
        except Exception:  # noqa: BLE001
            pass

    return emit


@dataclass(frozen=True)
class _Fold:
    """A folding verdict for one pegRNA, kept by the pegRNA's sequence."""

    error: Optional[str]
    pairs: int = 0
    helix: str = ""
    scaffold_nt: int = 0
    metrics: Optional[dict] = None


@dataclass(frozen=True)
class _Dock:
    """A docking verdict for one pegRNA, kept by the pegRNA's sequence."""

    error: Optional[str]
    ok: bool = False
    why: str = ""
    metrics: Optional[dict] = None


def _tag(side: str, cand: SideCandidate) -> str:
    return f"pegRNA {side} (PBS {cand.extension.pbs_len}, flap {cand.extension.flap_len})"


class _Walker:
    """Runs the gates for one pair at a time and remembers every per-pegRNA verdict.

    A pegRNA depends on its guide, its flap and PBS lengths and the donor end it copies, so the same pegRNA turns
    up in many pairs (every pair that shares guide A and its arm shares pegRNA A). Folding is cheap and docking
    is a ~40 s paid request, so a verdict is computed once per distinct sequence and reused.
    """

    def __init__(
        self,
        *,
        genome: str,
        site_s: int,
        cargo: str,
        cds_mask: Sequence[bool],
        splice_scorer: Any,
        maxent_scorer: Any,
        full_structure_scorer: Any,
        dock_scorer: Any,
        cas9: str,
        max_interface_drop: float,
        structure_threshold: float,
        tail: str,
        max_run: int,
        max_offtargets: int,
        offtarget_scorer: Any,
        emit: Callable[..., None],
    ) -> None:
        self.genome, self.site_s, self.cargo, self.cds_mask = genome, site_s, cargo, cds_mask
        self.splice_scorer, self.maxent_scorer = splice_scorer, maxent_scorer
        self.fold_scorer, self.dock_scorer, self.cas9 = full_structure_scorer, dock_scorer, cas9
        self.max_drop, self.threshold, self.tail, self.max_run = (
            max_interface_drop,
            structure_threshold,
            tail,
            max_run,
        )
        self.max_offtargets, self.offtarget_scorer, self.emit = (
            max_offtargets,
            offtarget_scorer,
            emit,
        )
        self.folds: Dict[str, _Fold] = {}
        self.docks: Dict[str, _Dock] = {}

    # -- per-pegRNA, memoized ------------------------------------------------------------------------

    def fold(self, cand: SideCandidate, guide: Any) -> _Fold:
        peg = cand.pegrna
        if peg in self.folds:
            return self.folds[peg]
        regions = pa_regions(guide.spacer, cand.extension.rtt, cand.extension.pbs, self.tail)
        try:
            # Fold once: the gate and the scaffold advisory read the same probability map.
            pair_map = self.fold_scorer.pair_probs(peg)
            pairing = functional_region_pairing(
                peg,
                regions=regions,
                scorer=self.fold_scorer,
                threshold=self.threshold,
                pairs=pair_map,
            )
            scaffold = scaffold_extension_pairing(
                peg,
                regions=regions,
                scorer=self.fold_scorer,
                threshold=self.threshold,
                pairs=pair_map,
            )
        except (
            Exception
        ) as exc:  # noqa: BLE001 - a scorer crash skips the candidate, with a diagnostic
            verdict = _Fold(error=str(exc))
        else:
            metrics = {
                "metric": "cross-region pairs",
                "value": pairing.cross_region_pairs,
                "cutoff": 0,
                "cmp": "==",
                "extra": (
                    f"base-pair prob ≥ {self.threshold} · extension pairs with {scaffold.paired_nt} scaffold nt"
                    + (f" · {pairing.describe()}" if pairing.cross_region_pairs else "")
                ),
            }
            verdict = _Fold(
                None, pairing.cross_region_pairs, pairing.describe(), scaffold.paired_nt, metrics
            )
        self.folds[peg] = verdict
        return verdict

    def dock(self, peg: str) -> _Dock:
        if peg in self.docks:
            return self.docks[peg]
        metrics: Dict[str, Any] = {}
        try:
            ok, why = passes_docking_gate(
                peg, self.cas9, scorer=self.dock_scorer, max_drop=self.max_drop, metrics=metrics
            )
        except DockingServiceError:
            raise  # an outage is not a verdict on this pegRNA (see its docstring), so it is not cached either
        except Exception as exc:  # noqa: BLE001
            verdict = _Dock(error=str(exc))
        else:
            verdict = _Dock(None, ok, why, metrics)
        self.docks[peg] = verdict
        return verdict

    # -- one pair ---------------------------------------------------------------------------------------

    def _fold_side(
        self,
        side: str,
        guide: Any,
        cands: Sequence[SideCandidate],
        idx: int,
        lab: str,
        diag: List[str],
    ) -> List[Tuple[SideCandidate, _Fold]]:
        """Fold every candidate of one side (cheap); report one verdict for the side."""
        passing: List[Tuple[SideCandidate, _Fold]] = []
        last: Optional[Tuple[SideCandidate, _Fold]] = None
        for cand in cands:
            verdict = self.fold(cand, guide)
            last = (cand, verdict)
            if verdict.error is not None:
                diag.append(f"{lab}: {_tag(side, cand)}: structure scorer error: {verdict.error}")
            elif verdict.pairs > 0:
                diag.append(
                    f"{lab}: {_tag(side, cand)}: cross-region pairing ({_folder_name(self.fold_scorer)}): "
                    f"{verdict.helix or str(verdict.pairs) + ' pair(s)'}"
                )
            else:
                passing.append((cand, verdict))
        gate = f"structure_{side.lower()}"
        shown = passing[0] if passing else last
        if shown is not None:
            cand, v = shown
            if v.error is not None:
                self.emit(idx, gate, f"structure scorer error: {v.error}", False, None)
            else:
                info = f"{v.pairs} cross-region pairs (p>={self.threshold})" + (
                    f": {v.helix}" if v.pairs else ""
                )
                self.emit(idx, gate, info, v.pairs == 0, v.metrics)
        return passing

    def _dock_side(
        self,
        side: str,
        passing: Sequence[Tuple[SideCandidate, _Fold]],
        idx: int,
        lab: str,
        diag: List[str],
    ) -> Optional[Tuple[SideCandidate, _Fold, _Dock]]:
        gate = f"docking_{side.lower()}"
        last: Optional[_Dock] = None
        for cand, fold in passing:
            verdict = self.dock(cand.pegrna)
            last = verdict
            if verdict.error is not None:
                diag.append(f"{lab}: {_tag(side, cand)}: docking scorer error: {verdict.error}")
                continue
            if verdict.ok:
                self.emit(idx, gate, verdict.why, True, verdict.metrics)
                return cand, fold, verdict
            diag.append(f"{lab}: {_tag(side, cand)}: docking: {verdict.why}")
        if last is not None:
            if last.error is not None:
                self.emit(idx, gate, f"docking scorer error: {last.error}", False, None)
            else:
                self.emit(idx, gate, last.why, False, last.metrics)
        return None

    def try_pair(self, plan: PairPlan, idx: int, diag: List[str]) -> Optional[Dict[str, Any]]:
        """Run every gate for one pair. Returns what the result needs, or ``None`` with the reasons in ``diag``."""
        a, b = plan.pair.a, plan.pair.b
        lab = f"pair {a.cut}/{b.cut}"

        # 1. folding, both sides, the cheapest model gate
        ok_a = self._fold_side("A", a, plan.candidates_a, idx, lab, diag)
        ok_b = self._fold_side("B", b, plan.candidates_b, idx, lab, diag)
        if not ok_a or not ok_b:
            return None

        # 2. the product-level splice gates (they do not depend on PBS or flap lengths)
        info: Dict[str, str] = {}
        splice_advisory: Optional[str] = None
        if self.maxent_scorer is not None:
            ms: Dict[str, Any] = {}
            try:
                ms_ok, ms_why = passes_cryptic_splice_gate(
                    self.genome, plan.product, scorer=self.maxent_scorer, metrics=ms
                )
            except Exception as exc:  # noqa: BLE001
                diag.append(f"{lab}: cryptic-splice scorer error: {exc}")
                return None
            self.emit(idx, "cryptic", ms_why, ms_ok, ms)
            if not ms_ok:
                diag.append(f"{lab}: cryptic splice: {ms_why}")
                return None
            info["cryptic"] = ms_why

        sp: Dict[str, Any] = {}
        try:
            sp_ok, sp_why = passes_splice_gate(
                self.genome, plan.product, scorer=self.splice_scorer, metrics=sp
            )
        except ValueError:
            info["splice"] = "advisory skip (needs >=10kb context; MaxEnt is primary)"
            self.emit(
                idx,
                "splice",
                "advisory skip (<10kb context)",
                True,
                {"metric": "SpliceAI Δ", "value": None, "note": "advisory skip (<10kb context)"},
            )
        except Exception as exc:  # noqa: BLE001
            diag.append(f"{lab}: splice cross-check error: {exc}")
            return None
        else:
            binding = bool(getattr(self.splice_scorer, "binding", True))
            where = str(getattr(self.splice_scorer, "context_note", "") or "")
            if binding:
                self.emit(idx, "splice", sp_why, sp_ok, sp)
                if not sp_ok:
                    diag.append(f"{lab}: splice (SpliceAI cross-check): {sp_why}")
                    return None
                info["splice"] = f"{sp_why} ({where})" if where else sp_why
            else:
                advisory = _advisory_splice_text(sp_ok, sp, where, primary_passed="cryptic" in info)
                info["splice"] = advisory
                if not sp_ok:
                    splice_advisory = advisory
                self.emit(
                    idx,
                    "splice",
                    advisory,
                    True,
                    {
                        "metric": "SpliceAI Δ",
                        "value": sp.get("value"),
                        "note": f"SpliceAI Δ {sp.get('value')}: advisory, not binding ({where})",
                    },
                )

        # 3. docking, the expensive gate, last; A before B so a failure on A spares B
        won_a = self._dock_side("A", ok_a, idx, lab, diag)
        if won_a is None:
            return None
        won_b = self._dock_side("B", ok_b, idx, lab, diag)
        if won_b is None:
            return None
        cand_a, fold_a, dock_a = won_a
        cand_b, fold_b, dock_b = won_b

        # 4. the two chosen pegRNAs together: the ordered sequences must reproduce the product through the strand
        #    simulator. (That the donor's two ends cannot pair with each other was settled when the pair was planned:
        #    every candidate of a side has the same flap length, because the TTTT rule on the PBS does not depend on
        #    the flap and the RTT rules do not depend on the PBS.)
        fa, fb = cand_a.extension.flap_len, cand_b.extension.flap_len
        try:
            simulate_pa(
                genome=self.genome,
                a=a,
                b=b,
                donor=plan.donor,
                peg_a=cand_a.pegrna,
                peg_b=cand_b.pegrna,
                flap_a=fa,
                pbs_a=cand_a.extension.pbs_len,
                flap_b=fb,
                pbs_b=cand_b.extension.pbs_len,
                tail_len=len(self.tail),
                expected_product=plan.product,
            )
        except PAMismatch as exc:
            diag.append(f"{lab}: PA strand mismatch: {exc}")
            return None
        return {
            "cand_a": cand_a,
            "cand_b": cand_b,
            "fold_a": fold_a,
            "fold_b": fold_b,
            "dock_a": dock_a,
            "dock_b": dock_b,
            "info": info,
            "splice_advisory": splice_advisory,
        }


# ---- assembling the answer --------------------------------------------------------------------------


def _tier_text(offtarget_scorer: Any) -> str:
    tier = getattr(offtarget_scorer, "count_mismatch_threshold", None)
    return f"within {tier} mismatches" if tier is not None else "in the mismatch tier"


def _arm_change_rows(
    genome: str, plan: PairPlan, cds_mask: Sequence[bool]
) -> Tuple[BaseChange, ...]:
    """Every base the donor's arms change from the window, each tagged with the arm it is in, and with the guide
    whose PAM it breaks (or only weakens) where that is why it was made."""
    verbs = {"broken": "breaks", "weakened": "weakens"}
    pam_edits = {
        o.position: (o.side, verbs[o.status])
        for o in plan.pam_outcomes
        if o.status in verbs and o.position is not None
    }
    rows: List[BaseChange] = []
    for arm, label in ((plan.arm5, "donor arm 5'"), (plan.arm3, "donor arm 3'")):
        if not arm.sense or not arm.n_edits:
            continue
        for c in arm_changes(genome, arm.lo, arm.hi, arm.sense, cds_mask, part=label):
            edit = pam_edits.get(c.position - 1)
            if edit:
                side, verb = edit
                c = replace(c, part=f"{label}, {verb} the PAM of guide {side}")
            rows.append(c)
    return tuple(rows)


def _kept_arm_lines(plan: PairPlan) -> List[str]:
    lines: List[str] = []
    for arm, label in ((plan.arm5, "donor arm 5'"), (plan.arm3, "donor arm 3'")):
        if arm.locked_span is not None:
            a, b = arm.locked_span
            lines.append(
                f"{label} (nt {a + 1}-{b} of the window, {b - a} nt): non-coding DNA copied exactly "
                "(Decision 15: only protein-coding bases are ever recoded)"
            )
    return lines


def _recode_info(plan: PairPlan) -> str:
    total = plan.arm5.n_edits + plan.arm3.n_edits
    bits = []
    for side, outcome in zip("AB", plan.pam_outcomes):
        if outcome.status == "broken" and outcome.position is not None:
            bits.append(
                f"guide {side}'s PAM broken at nt {outcome.position + 1} by one synonymous change"
            )
        elif outcome.status == "weakened" and outcome.position is not None:
            bits.append(
                f"guide {side}'s PAM weakened to {outcome.pam} at nt {outcome.position + 1} (the only "
                "synonymous change there; a weak PAM, see the warnings)"
            )
        elif outcome.status == "weakened":
            bits.append(
                f"guide {side}'s PAM left weak ({outcome.pam}) by the recode (see the warnings)"
            )
        elif outcome.status == "broken-by-recode":
            bits.append(f"guide {side}'s PAM broken by the recode")
        else:
            bits.append(f"guide {side}'s PAM regenerated (see the warnings)")
    return (
        f"{total} synonymous nucleotide edit(s) in the donor arms (5' arm {plan.arm5.n_edits}, "
        f"3' arm {plan.arm3.n_edits}); " + "; ".join(bits)
    )


def _build_result(
    plan: PairPlan,
    won: Dict[str, Any],
    *,
    spec: Any,
    genome: str,
    site_s: int,
    cds_mask: Sequence[bool],
    cargo: str,
    tail: str,
    max_run: int,
    max_offtargets: int,
    offtarget_scorer: Any,
    structure_threshold: float,
    importance_info: str,
    candidates_tried: int,
    pairs_considered: int,
) -> RecodedPAResult:
    a, b = plan.pair.a, plan.pair.b
    cand_a, cand_b = won["cand_a"], won["cand_b"]
    fa, fb = cand_a.extension.flap_len, cand_b.extension.flap_len
    run_a = genomic_run_under_flap("A", genome, a.cut, b.cut, plan.donor, fa)
    run_b = genomic_run_under_flap("B", genome, a.cut, b.cut, plan.donor, fb)

    warnings: List[str] = [STANDING_CAVEAT]
    warnings.extend(
        _advisories(
            plan.pam_outcomes, (a, b), (run_a, run_b), (fa, fb), plan.spacer_advisories, max_run
        )
    )
    if won["splice_advisory"]:
        warnings.append("splice cross-check: " + won["splice_advisory"])
    for side, cand, fold in (("A", cand_a, won["fold_a"]), ("B", cand_b, won["fold_b"])):
        if fold.scaffold_nt >= THRESHOLDS.scaffold_pairing_advisory_nt:
            warnings.append(
                f"pegRNA {side}'s extension pairs with {fold.scaffold_nt} scaffold nucleotides "
                f"(p >= {structure_threshold}; flagged from {THRESHOLDS.scaffold_pairing_advisory_nt}). Ordinary "
                "designs pair with far fewer; designs built to tie up the scaffold pair with 18 or more. Neither "
                "the folding gate nor the docking gate reliably catches this, so check the predicted structure "
                "before ordering"
            )
    for side, dock in (("A", won["dock_a"]), ("B", won["dock_b"])):
        metrics = dock.metrics or {}
        if metrics.get("borderline"):
            warnings.append(f"docking gate (pegRNA {side}): {metrics.get('note', dock.why)}")
        if metrics.get("advisory"):
            warnings.append(f"docking gate is advisory-only for pegRNA {side}: {dock.why}")
    tier = _tier_text(offtarget_scorer)
    for side, n in zip("AB", plan.offtargets):
        if n > 0:
            warnings.append(
                f"guide {side} has {n} off-target site(s) {tier} (max_offtargets={max_offtargets}). Each is a "
                "single-strand nick, not a double-strand break, but check them before ordering"
            )

    rows = _arm_change_rows(genome, plan, cds_mask)
    report = build_change_report(spec, cargo, 0, rows)
    report = replace(report, arm_kept_exactly=tuple(_kept_arm_lines(plan)))

    info = won["info"]
    gate_info: Dict[str, str] = {
        "offtarget": (
            f"guide A (cut {a.cut}): {plan.offtargets[0]} off-target site(s) {tier}; "
            f"guide B (cut {b.cut}): {plan.offtargets[1]} off-target site(s) {tier}"
        ),
        "recode": _recode_info(plan),
        "inv0": "0 bp endogenous deletion; recode protein-silent (INV-0 held)",
        "strand": "both pegRNAs reproduce the donor ends and the product (simulate_pa)",
        **info,
        "importance": importance_info,
    }
    for side, fold in (("a", won["fold_a"]), ("b", won["fold_b"])):
        gate_info[f"structure_{side}"] = (
            f"{fold.pairs} cross-region pairs (threshold p>={structure_threshold}); extension pairs with "
            f"{fold.scaffold_nt} scaffold nt"
        )
    gate_info["docking_a"] = won["dock_a"].why
    gate_info["docking_b"] = won["dock_b"].why
    order = [
        "offtarget",
        "recode",
        "inv0",
        "strand",
        "cryptic",
        "splice",
        "importance",
        "structure_a",
        "structure_b",
        "docking_a",
        "docking_b",
    ]

    def guide(side, g, cand, outcome, offs, on, run) -> PAGuideDesign:
        return PAGuideDesign(
            side=side,
            strand=g.strand,
            spacer=g.spacer,
            cut=g.cut,
            pam=g.pam,
            pam_status=outcome.status,
            product_pam=outcome.pam,
            flap_len=cand.extension.flap_len,
            pbs_len=cand.extension.pbs_len,
            pbs=cand.extension.pbs,
            rtt=cand.extension.rtt,
            tail=tail,
            pegRNA_sequence=cand.pegrna,
            offtargets=offs,
            on_target=on,
            genomic_run=run,
        )

    parts = [("arm5", plan.arm5.sense), ("cargo", cargo), ("arm3", plan.arm3.sense)]
    return RecodedPAResult(
        guide_a=guide(
            "A", a, cand_a, plan.pam_outcomes[0], plan.offtargets[0], plan.on_target[0], run_a
        ),
        guide_b=guide(
            "B", b, cand_b, plan.pam_outcomes[1], plan.offtargets[1], plan.on_target[1], run_b
        ),
        site_s=site_s,
        span=(a.cut, b.cut),
        cargo=cargo,
        arm5=plan.arm5.sense,
        arm3=plan.arm3.sense,
        donor=plan.donor,
        product=plan.product,
        donor_parts=[(name, seq) for name, seq in parts if seq],
        recode_edits=plan.arm5.n_edits + plan.arm3.n_edits,
        deletion_size=0,
        candidates_tried=candidates_tried,
        pairs_considered=pairs_considered,
        warnings=warnings,
        gate_results=[(k, gate_info[k]) for k in order if k in gate_info],
        recode_spans=[s for s in (plan.arm5.recode_span, plan.arm3.recode_span) if s is not None],
        locked_spans=[s for s in (plan.arm5.locked_span, plan.arm3.locked_span) if s is not None],
        arm_changes=rows,
        changes=report,
        simulated=True,
    )


def run_recoded_pa(
    *,
    cargo: str,
    genome: str,
    site_s: int,
    cds_mask: Sequence[bool],
    # injected scorers (fakes in CI, real adapters with the models)
    splice_scorer: Any,
    full_structure_scorer: Any,
    importance_scorer: Any,
    dock_scorer: Any,
    offtarget_scorer: Callable[[List[str]], Dict[str, Tuple[float, int, bool]]],
    maxent_scorer: Any = None,
    pbs_lengths: Sequence[int] = (13, 15, 17),
    flap_length: int = 35,
    flap_min: int = 30,
    flap_max: int = 50,
    min_span: int = 30,
    max_span: int = 100,
    max_run: int = 8,
    tail: str = "UUU",
    cas9: str = "nCas9",
    structure_threshold: Optional[float] = None,
    max_interface_drop: float = 0.15,
    protein_length: Optional[int] = None,
    insertion_importance_context: int = 1,
    max_offtargets: int = 0,
    min_ontarget: float = 0.5,
    max_pairs: int = 30,
    cargo_locked: Optional[Sequence[Sequence[int]]] = None,
    cargo_parts: Optional[Sequence[Dict[str, Any]]] = None,
    on_gate: Optional[Callable[..., None]] = None,
) -> RecodedPAResult:
    """Design a Prime Assembly insertion: two pegRNAs, one donor, the predicted product.

    Raises :class:`NoValidPADesign` (fail loud, with per-pair diagnostics) when no pair survives every gate or the
    cargo is too short for assembly, and :class:`ValueError` for arguments that make no sense. The cargo is used
    exactly as given; ``cargo_locked`` / ``cargo_parts`` are validated and reported (ADR 0029) though nothing in a
    donor ever recodes the cargo. Only the genomic arms are recoded, and only where the bases code for protein.
    """
    cargo = (cargo or "").upper()
    genome = genome.upper()
    _check_args(
        cargo=cargo,
        genome=genome,
        site_s=site_s,
        cds_mask=cds_mask,
        pbs_lengths=pbs_lengths,
        flap_length=flap_length,
        flap_min=flap_min,
        flap_max=flap_max,
        min_span=min_span,
        max_span=max_span,
        max_run=max_run,
        tail=tail,
        max_offtargets=max_offtargets,
        max_pairs=max_pairs,
    )
    spec = build_cargo_spec(
        cargo, parts=cargo_parts, locked=cargo_locked
    )  # ValueError on bad lock inputs
    _flap_ceiling(
        cargo, flap_length, flap_min, flap_max
    )  # the cargo rule: refuse before any model runs
    emit = _emitter(on_gate)

    # The insertion site does not depend on which guides are used, so its structural tolerance is judged once, first:
    # if S is a buried or special residue no pair can help, and there is no point screening guides.
    footprint = _insertion_site_footprint(
        site_s,
        cds_mask,
        context=insertion_importance_context,
        protein_length=_resolve_protein_length(protein_length, importance_scorer),
    )
    if footprint:
        try:
            blocked = span_blocked(importance_scorer, footprint)
        except Exception as exc:  # noqa: BLE001
            msg = f"importance scorer error: {exc}"
            raise NoValidPADesign(msg, diagnostics=[msg]) from exc
        emit(
            1,
            "importance",
            f"residues {footprint} "
            + ("blocked (buried/special)" if blocked else "tolerated (surface/flexible)"),
            not blocked,
            _importance_metrics(importance_scorer, footprint),
        )
        if blocked:
            msg = (
                f"insertion site S residue(s) {footprint} structurally important (buried/packed/special) -- "
                "inserting a tag/linker here risks disrupting the protein; choose a surface/loop/terminal "
                "insertion site (this does not depend on the guide pair, so no pair was tried)"
            )
            raise NoValidPADesign(f"no design: {msg}", diagnostics=[msg])
        importance_info = f"residues {footprint} tolerated (surface/flexible)"
    else:
        importance_info = "no flanking CDS at insertion site (n/a)"
        emit(1, "importance", "no flanking CDS (n/a)", True)

    # Warm the apo baseline for the legacy absolute-floor gate; a calibrated scorer never uses it. Fail here,
    # before the minutes-long guide screen, if the docking service cannot run at all.
    try:
        if not hasattr(dock_scorer, "interface_iptm_samples"):
            dock_scorer.apo_interface_iptm(cas9)
    except DockingServiceError:
        raise
    except Exception as exc:  # noqa: BLE001
        raise NoValidPADesign(f"docking apo baseline failed to compute: {exc}") from exc

    planning = plan_pairs(
        cargo=cargo,
        genome=genome,
        site_s=site_s,
        cds_mask=cds_mask,
        offtarget_scorer=offtarget_scorer,
        pbs_lengths=pbs_lengths,
        flap_length=flap_length,
        flap_min=flap_min,
        flap_max=flap_max,
        min_span=min_span,
        max_span=max_span,
        max_run=max_run,
        tail=tail,
        max_offtargets=max_offtargets,
        min_ontarget=min_ontarget,
    )
    diagnostics = list(planning.diagnostics)
    if not planning.plans:
        raise NoValidPADesign(
            f"no pair survived all gates after 0 tried ({planning.pairs_found} PAM-in pair(s) after the guide "
            f"screens); diagnostics: {'; '.join(diagnostics[:10])}",
            diagnostics=diagnostics,
        )

    threshold = resolve_structure_threshold(structure_threshold, full_structure_scorer)
    walker = _Walker(
        genome=genome,
        site_s=site_s,
        cargo=cargo,
        cds_mask=cds_mask,
        splice_scorer=splice_scorer,
        maxent_scorer=maxent_scorer,
        full_structure_scorer=full_structure_scorer,
        dock_scorer=dock_scorer,
        cas9=cas9,
        max_interface_drop=max_interface_drop,
        structure_threshold=threshold,
        tail=tail,
        max_run=max_run,
        max_offtargets=max_offtargets,
        offtarget_scorer=offtarget_scorer,
        emit=emit,
    )
    tried = 0
    for plan in planning.plans[:max_pairs]:
        tried += 1
        won = walker.try_pair(plan, tried, diagnostics)
        if won is not None:
            return _build_result(
                plan,
                won,
                spec=spec,
                genome=genome,
                site_s=site_s,
                cds_mask=cds_mask,
                cargo=cargo,
                tail=tail,
                max_run=max_run,
                max_offtargets=max_offtargets,
                offtarget_scorer=offtarget_scorer,
                structure_threshold=threshold,
                importance_info=importance_info,
                candidates_tried=tried,
                pairs_considered=planning.pairs_found,
            )
    raise NoValidPADesign(
        f"no pair survived all gates after {tried} tried ({planning.pairs_found} PAM-in pair(s) after the guide "
        f"screens); diagnostics: {'; '.join(diagnostics[:10])}",
        diagnostics=diagnostics,
    )
