"""Codon expansion and optimisation.

Core building blocks used across the pipeline:
- :data:`CODONS_BY_AA`  -- canonical codon map, cytosine-first per aa.
- :data:`HUMAN_CODON_FREQUENCIES` -- flat dict of {codon: fraction} from the
  ``h_sapiens_9606`` table shipped by ``python-codon-tables``.
- :func:`codons_for_aa` / :func:`codon_frequency` -- convenience accessors.
- :func:`optimise_bulk` -- pick one codon per residue using a weighted score.
- (added in later tasks) :func:`enumerate_variants`,
  :func:`enumerate_with_metadata`, quality filters, etc.
"""

from __future__ import annotations

import heapq
import itertools
import math
import re
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from typing import Dict, Iterator, List, Tuple

import python_codon_tables as _pct
import yaml

# ---------------------------------------------------------------------- tables


def _load_human_codon_table() -> Dict[str, Dict[str, float]]:
    return _pct.get_codons_table("h_sapiens_9606")


def _build_codons_by_aa() -> Dict[str, List[str]]:
    """Canonical codon lists per amino acid, sorted cytosine-first then lex.

    Cytosine-first ordering matters for downstream enumeration where the
    pipeline prefers C-rich codons (GC-content stability in the C-terminal
    edit window).
    """
    table = _load_human_codon_table()
    out: Dict[str, List[str]] = {}
    for aa, codon_map in table.items():
        codons = list(codon_map.keys())
        # Stable, cytosine-first, then lexicographic.
        codons.sort(key=lambda c: (0 if c.startswith("C") else 1, c))
        out[aa] = codons
    return out


CODONS_BY_AA: Dict[str, List[str]] = _build_codons_by_aa()


def _build_flat_frequency_table() -> Dict[str, float]:
    table = _load_human_codon_table()
    flat: Dict[str, float] = {}
    for aa, codon_map in table.items():
        for codon, freq in codon_map.items():
            flat[codon] = float(freq)
    return flat


HUMAN_CODON_FREQUENCIES: Dict[str, float] = _build_flat_frequency_table()


# ------------------------------------------------------------- CAI (ADR 0009)

# ``w = f(codon) / max_{c in family} f(c)`` -- the codon's relative adaptiveness
# against the human genome-wide table. Single-codon families (Met, Trp) and stop
# codons carry no discriminating signal and are excluded from the CAI geometric
# mean (they are identical across every variant of one aa sequence).
_CAI_ZERO_FLOOR = 0.01  # rounded table can read 0.00; floor keeps log finite


def _build_cai_tables() -> Tuple[Dict[str, float], frozenset]:
    table = _load_human_codon_table()
    rel: Dict[str, float] = {}
    excluded: set = set()
    for aa, codon_map in table.items():
        fmax = max(codon_map.values()) or 1.0
        single = len(codon_map) == 1
        for codon, freq in codon_map.items():
            cu = codon.upper()
            rel[cu] = (float(freq) / fmax) if fmax else 0.0
            if single or aa == "*":
                excluded.add(cu)
    return rel, frozenset(excluded)


_CAI_REL_ADAPT, _CAI_EXCLUDED = _build_cai_tables()


def human_cai(dna: str, frame: int = 0) -> float:
    """Codon Adaptation Index of ``dna`` against the human genome-wide table.

    Geometric mean of per-codon relative adaptiveness ``w = f / f_max_in_family``
    over the in-frame codons, excluding single-codon families (Met ``ATG``, Trp
    ``TGG``) and stop codons. ``w`` is floored at :data:`_CAI_ZERO_FLOOR` so a
    codon that the 2-decimal table rounds to ``0.00`` degrades smoothly instead of
    raising ``math domain error``.

    Empty input returns ``0.0``; a non-empty sequence of only excluded codons
    returns ``1.0`` (trivially maximal, no discriminating signal). Case-insensitive.
    Raises :class:`ValueError` on a ragged (non-multiple-of-3) length at the given
    frame, or on an unknown codon. The caller is responsible for frame-correctness.
    """
    if not dna:
        return 0.0
    seq = dna.upper()
    if frame:
        seq = seq[frame:]
    if len(seq) % 3 != 0:
        raise ValueError(f"human_cai: length {len(seq)} not a multiple of 3 at frame {frame}")
    logs: List[float] = []
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        if codon not in _CAI_REL_ADAPT:
            raise ValueError(f"human_cai: unknown codon {codon!r}")
        if codon in _CAI_EXCLUDED:
            continue
        logs.append(math.log(max(_CAI_REL_ADAPT[codon], _CAI_ZERO_FLOOR)))
    if not logs:
        return 1.0
    return math.exp(sum(logs) / len(logs))


# ---------------------------------------------------------------- CodonVariant


@dataclass(frozen=True)
class CodonVariant:
    """A candidate DNA variant for a given aa sequence.

    Attributes
    ----------
    dna:
        Uppercase ACGT string encoding the aa sequence.
    c_count:
        Number of cytosine residues in the DNA (higher -> preferred).
    gc_content:
        GC fraction in [0, 1].
    rank:
        Position in a sorted variant list (0 = best).
    """

    dna: str
    c_count: int
    gc_content: float
    rank: int


# --------------------------------------------------------------------- lookups


def codons_for_aa(aa: str) -> List[str]:
    """Return the cytosine-first list of codons encoding ``aa``.

    Case-insensitive. Raises :class:`KeyError` for unknown aa letters.
    """
    key = aa.upper()
    if key not in CODONS_BY_AA:
        raise KeyError(aa)
    return list(CODONS_BY_AA[key])


def codon_frequency(aa: str, codon: str) -> float:
    """Return the human frequency of ``codon`` within the synonymous group for ``aa``."""
    key = aa.upper()
    table = _load_human_codon_table()
    try:
        return float(table[key][codon.upper()])
    except KeyError as exc:
        raise KeyError(f"codon {codon!r} not listed for aa {aa!r}") from exc


# --------------------------------------------------------- variant enumeration


def _surviving_variants(aa_sequence: str) -> Iterator[Tuple[str, int, float]]:
    """Yield (dna, c_count, gc) for every synonymous DNA combination that
    passes all quality filters. Intentionally a generator so the caller can
    stream into heapq for top-K selection.
    """
    codon_choices = [CODONS_BY_AA[aa.upper()] for aa in aa_sequence]
    for combo in itertools.product(*codon_choices):
        dna = "".join(combo)
        if _has_internal_stop(dna):
            continue
        if _has_cpg(dna):
            continue
        if _has_restriction_site(dna):
            continue
        if _has_homopolymer(dna):
            continue
        yield dna, dna.count("C"), _gc_content(dna)


_RANKINGS = ("gc_c", "cai")


def _check_ranking(ranking: str) -> None:
    if ranking not in _RANKINGS:
        raise ValueError(f"ranking must be one of {_RANKINGS}; got {ranking!r}")


def _nlargest_key(ranking: str):
    """Key for :func:`heapq.nlargest` (larger = better) over ``(dna, c_count, gc)``."""
    if ranking == "cai":
        # CAI primary, then the gc_c order as a deterministic tiebreak.
        return lambda t: (human_cai(t[0]), t[1], -abs(t[2] - 0.55))
    return lambda t: (t[1], -abs(t[2] - 0.55))


def _final_sort_key(ranking: str):
    """Ascending sort key mirroring :func:`_nlargest_key`, with ``dna`` for stability."""
    if ranking == "cai":
        return lambda t: (-human_cai(t[0]), -t[1], abs(t[2] - 0.55), t[0])
    return lambda t: (-t[1], abs(t[2] - 0.55), t[0])


def enumerate_variants(
    aa_sequence: str, max_variants: int = 100, ranking: str = "gc_c"
) -> List[CodonVariant]:
    """Return top ``max_variants`` surviving variants, globally ranked.

    ``ranking="gc_c"`` (default): ``(c_count DESC, abs(gc_content - 0.55) ASC)`` --
    the synthesis-stability order. ``ranking="cai"``: human Codon Adaptation Index
    DESC (see :func:`human_cai`), with the gc_c order as a deterministic tiebreak,
    so the most human-adapted synonymous variant ranks first. Uses
    :func:`heapq.nlargest` so the top-K is globally correct even when only a subset
    of combinations is kept in memory.
    """
    _check_ranking(ranking)
    if not aa_sequence:
        return []

    top = heapq.nlargest(
        max_variants,
        _surviving_variants(aa_sequence),
        key=_nlargest_key(ranking),
    )

    # nlargest produces items in descending key order; ties may land arbitrarily,
    # so re-sort by the canonical ordering before assigning ranks (deterministic).
    top.sort(key=_final_sort_key(ranking))

    return [
        CodonVariant(dna=dna, c_count=c, gc_content=gc, rank=i)
        for i, (dna, c, gc) in enumerate(top)
    ]


@dataclass(frozen=True)
class EnumerationResult:
    """Return type of :func:`enumerate_with_metadata`.

    ``theoretical_max`` is the number of synonymous DNA combinations that
    exist before filters are applied (i.e. the pure product of per-aa codon
    counts). ``was_capped`` is True when at least one valid variant was
    discarded because the ``max_variants`` cap was reached.
    """

    variants: list
    theoretical_max: int
    was_capped: bool


def enumerate_with_metadata(
    aa_sequence: str, max_variants: int = 100, ranking: str = "gc_c"
) -> EnumerationResult:
    """Like :func:`enumerate_variants` but also returns counting metadata."""
    _check_ranking(ranking)
    if not aa_sequence:
        return EnumerationResult(variants=[], theoretical_max=0, was_capped=False)

    theoretical_max = 1
    for aa in aa_sequence:
        theoretical_max *= len(CODONS_BY_AA[aa.upper()])

    # We need to know whether more valid (filter-passing) variants existed
    # than the cap allowed, so pull one extra item past max_variants and
    # compare.
    peek = max_variants + 1
    picked = heapq.nlargest(
        peek,
        _surviving_variants(aa_sequence),
        key=_nlargest_key(ranking),
    )
    was_capped = len(picked) > max_variants
    picked = picked[:max_variants]
    picked.sort(key=_final_sort_key(ranking))
    variants = [
        CodonVariant(dna=dna, c_count=c, gc_content=gc, rank=i)
        for i, (dna, c, gc) in enumerate(picked)
    ]
    return EnumerationResult(
        variants=variants,
        theoretical_max=theoretical_max,
        was_capped=was_capped,
    )


# ---------------------------------------------------------- bulk optimisation


# ---------------------------------------------------------------- filters


_STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})


@lru_cache(maxsize=1)
def load_restriction_sites() -> List[str]:
    """Return the sorted list of recognition sequences bundled with the package."""
    raw = resources.files("pegrna_pipeline.data").joinpath("restriction_sites.yaml").read_text()
    data = yaml.safe_load(raw) or {}
    sites = sorted({e["recognition"].upper() for e in data.get("enzymes", [])})
    return sites


def _has_restriction_site(dna: str) -> bool:
    """Return True if ``dna`` contains any bundled restriction site as a substring."""
    upper = dna.upper()
    return any(site in upper for site in load_restriction_sites())


def _gc_content(dna: str) -> float:
    """Return the GC fraction of ``dna`` in [0.0, 1.0]. Empty input -> 0.0."""
    if not dna:
        return 0.0
    upper = dna.upper()
    return (upper.count("G") + upper.count("C")) / len(upper)


def _has_homopolymer(dna: str, min_run: int = 5) -> bool:
    """Return True if ``dna`` contains a homopolymer run of length >= min_run."""
    if min_run <= 0:
        return bool(dna)
    pattern = re.compile(rf"(.)\1{{{min_run - 1},}}")
    return pattern.search(dna) is not None


def _has_cpg(dna: str) -> bool:
    """Return True if ``dna`` contains a CpG dinucleotide (literal 'CG')."""
    return "CG" in dna


def _has_internal_stop(dna: str) -> bool:
    """Return True if ``dna`` contains any in-frame stop codon.

    Walks indices ``range(0, len(dna) - 2, 3)`` so every full codon is
    inspected -- including the last one (off-by-one fix).
    """
    for i in range(0, len(dna) - 2, 3):
        if dna[i : i + 3] in _STOP_CODONS:
            return True
    return False


# Penalty constants used by :func:`optimise_bulk`. Exposed as module-level
# constants so external callers (and tests) can reference or monkey-patch them.
_CPG_JUNCTION_PENALTY = 10.0
_RESTRICTION_SITE_PENALTY = 10.0
_HOMOPOLYMER_PENALTY = 10.0
_C_COUNT_BONUS = 0.3
_JUNCTION_WINDOW = 10  # characters to inspect around the junction for sites


def _window_has_site(prev: str, codon: str) -> bool:
    """Return True if the 10-nt window spanning the junction contains a
    bundled restriction site.
    """
    half = _JUNCTION_WINDOW - len(codon)
    half = max(half, 0)
    left = prev[-half:] if half else ""
    window = left + codon
    return _has_restriction_site(window)


def optimise_bulk(aa_sequence: str) -> str:
    """Weighted greedy codon optimiser.

    For each residue, each candidate codon is scored:

        score = codon_frequency(aa, codon)
              + 0.3 * codon.count("C")
              - 10.0 if prev DNA's final C concatenates with codon's leading G
                      (CpG at the exact junction)
              - 10.0 if the 10-nt window around the junction contains any
                      bundled restriction site
              - 10.0 if placing this codon creates a homopolymer run of
                      length >= 5

    The highest-scoring codon is picked; ties are broken by the
    cytosine-first ordering of :data:`CODONS_BY_AA`. A warning is logged if
    the best codon still violates one or more filters (i.e. no clean option
    exists at that position).

    An empty string returns ``""``. Unknown amino acids raise
    :class:`ValueError`.
    """
    if not aa_sequence:
        return ""

    import logging  # deferred to keep import cost low for simple users

    logger = logging.getLogger(__name__)

    dna_parts: List[str] = []
    for aa in aa_sequence:
        key = aa.upper()
        if key not in CODONS_BY_AA:
            raise ValueError(f"unknown amino acid {aa!r}")
        candidates = CODONS_BY_AA[key]
        prev = "".join(dna_parts)

        best_codon = None
        best_score = float("-inf")
        best_index = -1
        best_violations: List[str] = []

        for idx, codon in enumerate(candidates):
            score = codon_frequency(key, codon)
            score += _C_COUNT_BONUS * codon.count("C")

            violations: List[str] = []

            # CpG at junction.
            if prev and prev[-1:] == "C" and codon[:1] == "G":
                score -= _CPG_JUNCTION_PENALTY
                violations.append("cpg_junction")

            # Restriction site within junction window.
            if _window_has_site(prev, codon):
                score -= _RESTRICTION_SITE_PENALTY
                violations.append("restriction_site")

            # Homopolymer run after concatenation.
            if _has_homopolymer(prev + codon, min_run=5):
                score -= _HOMOPOLYMER_PENALTY
                violations.append("homopolymer")

            # Prefer higher score; tie-break favours earlier (cytosine-first) entry.
            if (score > best_score) or (score == best_score and idx < best_index):
                best_score = score
                best_codon = codon
                best_index = idx
                best_violations = violations

        assert best_codon is not None  # at least one candidate always exists
        if best_violations:
            logger.warning(
                "optimise_bulk: best codon %s for %s still violates %s",
                best_codon,
                key,
                best_violations,
            )
        dna_parts.append(best_codon)

    return "".join(dna_parts)


# ------------------------------------------ Stage B: scalable top-N by CAI (ADR 0015)
#
# Ranking synonymous variants by CAI is separable: CAI = exp((1/K) * sum_i log w_i)
# with K (count of non-excluded codons) fixed by the aa sequence, so CAI-descending
# order == descending sum of per-codon log-w. The top-N is therefore a k-best-over-
# independent-choices problem (best-first / Lawler-style k-best on the codon-index
# lattice), solved lazily so a 40 bp end (~10^5-10^7 combos) or a 500 bp middle never
# needs full enumeration. Quality filters (stop/CpG/restriction/homopolymer) are
# non-separable, so they gate YIELDING (not exploration) as a lazy post-filter; since
# the stream is CAI-descending, filtering preserves order.


def _cai_weight(codon: str) -> float:
    """Floored relative adaptiveness used for CAI ordering. Excluded codons
    (single-codon families, stops) are CAI-neutral (weight 1.0 -> log 0), matching
    :func:`human_cai`, which drops them from the geometric mean."""
    cu = codon.upper()
    if cu in _CAI_EXCLUDED:
        return 1.0
    return max(_CAI_REL_ADAPT.get(cu, 0.0), _CAI_ZERO_FLOOR)


def _cai_choices(aa_sequence: str) -> Tuple[List[List[str]], List[List[float]]]:
    """Per-position codon lists sorted by CAI weight desc (codon asc tiebreak), with
    the matching per-choice log-weights. Raises :class:`KeyError` on an unknown aa."""
    choices: List[List[str]] = []
    logws: List[List[float]] = []
    for aa in aa_sequence:
        key = aa.upper()
        if key not in CODONS_BY_AA:
            raise KeyError(aa)
        ordered = sorted(CODONS_BY_AA[key], key=lambda c: (-_cai_weight(c), c))
        choices.append(ordered)
        logws.append([math.log(_cai_weight(c)) for c in ordered])
    return choices, logws


def _cai_kbest(
    aa_sequence: str, *, max_examined: int, max_states: int, stats: "dict | None" = None
) -> Iterator[Tuple[str, int, float]]:
    """Lazily yield ``(dna, c_count, gc)`` for filter-passing synonymous variants of
    ``aa_sequence`` in CAI-descending order, via best-first (Lawler-style) k-best over
    the codon-index lattice.

    The heap key is ``(-score, state)``; ``score`` (sum of per-codon log-w) is carried
    incrementally (child = parent - old_logw + new_logw, clamped non-increasing so pop
    order stays monotone under float round-off), and the DNA is materialised only when a
    state is popped (for the filter check and yield) -- so a pushed successor costs no
    full rebuild. ``state`` is unique per variant, giving a total, deterministic order.
    Exploration is bounded by ``max_examined`` (time) and ``max_states`` (memory, = size
    of the visited set); when either cap is hit, ``stats['truncated']`` is set (if
    ``stats`` is given) and the generator stops, so a filter-hostile locus terminates
    instead of walking the whole lattice.
    """
    if not aa_sequence:
        return
    choices, logws = _cai_choices(aa_sequence)
    length = len(choices)
    start = (0,) * length
    start_score = math.fsum(logws[i][0] for i in range(length))
    heap: List[Tuple[float, Tuple[int, ...]]] = [(-start_score, start)]
    visited = {start}
    examined = 0
    while heap:
        neg_score, state = heapq.heappop(heap)
        examined += 1
        parent_score = -neg_score
        # Explore successors regardless of this state's filter outcome (a filtered state
        # may sit between two passing ones on the lattice). Score is incremental.
        for i in range(length):
            nxt = state[i] + 1
            if nxt < len(choices[i]):
                child = state[:i] + (nxt,) + state[i + 1 :]
                if child not in visited:
                    visited.add(child)
                    # Non-increasing by construction (choices sorted by weight desc);
                    # clamp guards float round-up so the heap pops in monotone CAI order.
                    child_score = min(
                        parent_score, parent_score - logws[i][state[i]] + logws[i][nxt]
                    )
                    heapq.heappush(heap, (-child_score, child))
        dna = "".join(choices[i][state[i]] for i in range(length))
        if not (
            _has_internal_stop(dna)
            or _has_cpg(dna)
            or _has_restriction_site(dna)
            or _has_homopolymer(dna)
        ):
            yield dna, dna.count("C"), _gc_content(dna)
        if examined >= max_examined or len(visited) >= max_states:
            if stats is not None:
                stats["truncated"] = 1
            return


def rank_insert_variants_by_cai(
    aa_sequence: str,
    *,
    max_variants: int = 2000,
    max_examined: int | None = None,
    max_states: int | None = None,
) -> List[CodonVariant]:
    """Return up to ``max_variants`` synonymous variants of ``aa_sequence`` in
    CAI-descending order (most human-adapted **filter-passing** variant first) via the
    lazy k-best search (no full enumeration). Within an exact CAI tie the order is a
    deterministic function of the codon indices (not a c-count/GC preference).

    ``max_examined`` (default ``max(200_000, 200*max_variants)``, time) and ``max_states``
    (default scales inversely with insert length to bound memory) cap exploration; under
    very restrictive quality filters the result may contain fewer than ``max_variants``
    and a warning is logged. Empty input -> ``[]``. Raises :class:`KeyError` on an unknown
    amino-acid letter (validated up front, independent of ``max_variants``).
    """
    if not aa_sequence:
        return []
    _cai_choices(aa_sequence)  # validate the aa alphabet up front (raises KeyError)
    if max_variants <= 0:
        return []
    length = len(aa_sequence)
    if max_examined is None:
        max_examined = max(200_000, 200 * max_variants)
    if max_states is None:
        max_states = max(100_000, 20_000_000 // max(length, 1))
    stats = {"truncated": 0}
    out: List[CodonVariant] = []
    for dna, c_count, gc in _cai_kbest(
        aa_sequence, max_examined=max_examined, max_states=max_states, stats=stats
    ):
        out.append(CodonVariant(dna=dna, c_count=c_count, gc_content=gc, rank=len(out)))
        if len(out) >= max_variants:
            break
    if len(out) < max_variants and stats["truncated"]:
        import logging  # noqa: PLC0415

        logging.getLogger(__name__).warning(
            "rank_insert_variants_by_cai: returned %d < %d requested for a %d-codon "
            "insert; exploration cap hit under restrictive filters. Raise max_examined/"
            "max_states or relax the insert if more variants are needed.",
            len(out),
            max_variants,
            length,
        )
    return out


def select_top_n_for_pridict(
    ranked: List[CodonVariant], n: int = 2000, *, lo: int = 1000, hi: int = 5000
) -> List[CodonVariant]:
    """Select the top ``n`` CAI-ranked variants to hand to PRIDICT2, with ``n`` clamped
    to ``[lo, hi]`` (default 1000-5000, 2000). Never invents variants: returns at most
    ``len(ranked)``. Assumes ``ranked`` is already CAI-ordered (from
    :func:`rank_insert_variants_by_cai`)."""
    n = max(lo, min(hi, int(n)))
    return list(ranked[:n])


def split_insert(insert_dna: str, *, end_len: int = 40) -> Tuple[str, str, str]:
    """Split an insert into ``(first_end, middle, last_end)`` such that
    ``first_end + middle + last_end == insert_dna`` exactly.

    The two ends (up to ``end_len`` nt, PRIDICT2's practical reach) are the parts a
    PRIDICT design touches; the middle is CAI-ranked independently (PRIDICT never sees
    it). For an insert shorter than ``2 * end_len`` the middle is empty and the ends do
    not overlap; for one ``<= end_len`` the whole insert is the first end. This is a
    plain byte split. The middle's reading frame within the insert is
    ``len(first_end) % 3`` -- 0 only when ``end_len`` is a multiple of 3; the default
    ``end_len=40`` is NOT, so the middle is generally frame-shifted. A caller CAI-ranking
    the middle MUST translate it at that frame (handled in Stage C) or it will silently
    score the wrong reading frame."""
    n = len(insert_dna)
    if n <= end_len:
        return insert_dna, "", ""
    if n <= 2 * end_len:
        return insert_dna[:end_len], "", insert_dna[end_len:]
    return insert_dna[:end_len], insert_dna[end_len : n - end_len], insert_dna[n - end_len :]
