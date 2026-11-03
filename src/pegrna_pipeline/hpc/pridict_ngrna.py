"""PRIDICT-driven ngRNA selection for the recoded TJ-PE auto path (wiring spec 2026-07-10).

The ngRNA (pegRNA#2 spacer) is sourced from PRIDICT2's ``*_nicking_guides.csv`` rows
-- returned per sequence under the ``"nicking"`` key of
:meth:`pegrna_pipeline.hpc.pridict_local.PridictLocalRunner.run_batch` -- while the
PIPELINE keeps ownership of *where* it nicks:

* :func:`nicking_rows_to_ngrnas` maps the raw nicking dicts to
  :class:`pegrna_pipeline.hpc.data_types.NgRNADesign`.
* :func:`build_pridict_ngrnas` selects among them by DISTANCE WINDOW + opposite strand
  + DeepCas9 score (reusing :func:`pegrna_pipeline.hpc.tjpe_workflow.select_tjpe_ngrna`),
  gates each choice on the injected off-target scorer, and on a gate failure advances to
  the NEXT nicking candidate (the NICK lever, per the unified plan's gate->lever map).
  rcPBS2 stays genome-derived via
  :func:`pegrna_pipeline.hpc.tjpe_workflow.derive_rc_pbs2` (``rcPBS2 = RC(TS)`` over the
  window 5' of the nick; plan D4). The surviving, off-target-clean candidates are
  returned best-first as :class:`InjectedNgRNA` records for
  :func:`pegrna_pipeline.hpc.recoded_tjpe_workflow.run_recoded_tjpe`.
* :func:`resolve_pridict_context` picks the sequence window handed to PRIDICT for the
  pegRNA#1 (and therefore the nicking) design: the provided window when it already has
  ``>= min_flank`` bp each side of the insert, else a wider GRCh38 window auto-fetched
  via :func:`pegrna_pipeline.hpc.genome_ingest.fetch_genomic_context` with the provided
  window localised inside it by exact substring match. It is fail-loud when neither
  yields adequate flanks.

This module imports no heavy dependency at load: the default genome fetcher is imported
lazily so CI (fake runner + fake fetcher) never touches the network.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from pegrna_pipeline.hpc.data_types import NgRNADesign
from pegrna_pipeline.hpc.nick_search import nick_usable_for_arm, offtarget_rejection
from pegrna_pipeline.hpc.tjpe_workflow import derive_rc_pbs2, select_tjpe_ngrna

__all__ = [
    "InjectedNgRNA",
    "nicking_rows_to_ngrnas",
    "build_pridict_ngrnas",
    "resolve_pridict_context",
]

#: The four columns :func:`pegrna_pipeline.hpc.pridict_parser.parse_ngrna_csv` requires;
#: a nicking row missing any of them is a hard error (never silently dropped).
_NICKING_REQUIRED = (
    "Nicking-Protospacer",
    "Nicking-Position-to-edit",
    "Target_Strand",
    "DeepCas9score",
)
#: PAM-disrupt column name (PRIDICT 2.0 vs pre-2.0), matched in priority order.
_PAM_DISRUPT_COLS = ("Nicking-PAMdisrupt", "PAM_disrupted")


@dataclass(frozen=True)
class InjectedNgRNA:
    """A pipeline-selected, off-target-clean PRIDICT ngRNA ready for injection.

    ``nick_position`` is the 0-indexed genome-window coordinate of the second nick
    (``site_s + distance_to_edit``); ``protospacer`` is the PRIDICT nicking protospacer
    (which becomes the ngRNA spacer); ``rcpbs2`` is the genome-derived
    ``derive_rc_pbs2`` window 5' of the nick.
    """

    nick_position: int
    protospacer: str
    rcpbs2: str
    distance_to_edit: int
    target_strand: str
    deepcas9_score: float
    #: Off-target sites the screen counted for this guide (0 unless a tolerance was granted).
    offtargets: int = 0


def _coerce_flag(val: object) -> bool:
    """Map a PRIDICT boolean-ish flag string to ``bool`` (leading ``no_`` -> False)."""
    s = str(val).strip().lower()
    if not s or s.startswith("no_") or s in {"false", "0", "no", "n", "na", "nan"}:
        return False
    return True


def nicking_rows_to_ngrnas(rows: Sequence[dict]) -> List[NgRNADesign]:
    """Map PRIDICT ``*_nicking_guides.csv`` rows (``csv.DictReader`` dicts) to
    :class:`NgRNADesign`.

    Mirrors :func:`pegrna_pipeline.hpc.pridict_parser.parse_ngrna_csv`'s field
    contract; a row missing any required key raises :class:`ValueError` (fail-loud).
    """
    out: List[NgRNADesign] = []
    for i, row in enumerate(rows):
        missing = [k for k in _NICKING_REQUIRED if k not in row]
        if missing:
            raise ValueError(f"nicking row {i} missing required key(s): {', '.join(missing)}")
        pam_col = next((c for c in _PAM_DISRUPT_COLS if c in row), None)
        out.append(
            NgRNADesign(
                spacer=str(row["Nicking-Protospacer"]),
                distance_to_edit=int(row["Nicking-Position-to-edit"]),
                target_strand=str(row["Target_Strand"]),
                deepcas9_score=float(row["DeepCas9score"]),
                pe3b=_coerce_flag(row["PE3b"]) if "PE3b" in row else False,
                pam_disrupt=_coerce_flag(row[pam_col]) if pam_col else False,
            )
        )
    return out


def _rank_admissible(
    ngrnas: Sequence[NgRNADesign],
    *,
    min_distance: int,
    max_distance: int,
    pegrna_strand: Optional[str],
) -> List[NgRNADesign]:
    """Admissible ngRNAs (distance window + opposite strand + spacer sanity), best-first.

    Reuses :func:`select_tjpe_ngrna` (which picks the single highest-DeepCas9 survivor)
    repeatedly, removing the winner each round, so the returned order is exactly the
    NICK-lever order: try the strongest guide first, then the next, and so on.
    """
    pool = list(ngrnas)
    ordered: List[NgRNADesign] = []
    while True:
        best = select_tjpe_ngrna(
            pool,
            min_distance=min_distance,
            max_distance=max_distance,
            pegrna_strand=pegrna_strand,
        )
        if best is None:
            break
        ordered.append(best)
        pool.remove(best)
    return ordered


def build_pridict_ngrnas(
    nicking_rows: Sequence[dict],
    *,
    genome: str,
    site_s: int,
    pegrna_strand: str,
    ngrna_distance_min: int,
    ngrna_distance_max: int,
    offtarget_scorer: Callable[[List[str]], Dict[str, Tuple[float, int, bool]]],
    pbs2_length: int = 17,
    min_ontarget: float = 0.5,
    max_offtargets: int = 0,
    arm_side: str = "auto",
) -> List[InjectedNgRNA]:
    """Select + off-target-gate PRIDICT nicking guides into injectable ngRNAs (best-first).

    Selection order (mirrors :func:`pegrna_pipeline.hpc.nick_search.screen_nick_candidates`
    off-target semantics):

    1. Admissible = distance window + opposite strand + spacer sanity (via
       :func:`_rank_admissible` / :func:`select_tjpe_ngrna`), ordered by DeepCas9 desc.
    2. Arm side: a nick on the wrong side of S for ``arm_side`` is dropped HERE, before
       the (slow) off-target scan -- see :func:`nick_usable_for_arm` for why this must
       happen when the pool is built and not only when it is walked.
    3. Off-target gate on each protospacer (one batched scorer call), by
       :func:`pegrna_pipeline.hpc.nick_search.offtarget_rejection`: a candidate with more than
       ``max_offtargets`` in-tier hits, or on-target below ``min_ontarget``, is skipped -- the
       NICK lever advances to the next candidate. ``max_offtargets`` defaults to 0 (ADR 0010).
    4. rcPBS2 = :func:`derive_rc_pbs2` over the window 5' of ``site_s + distance_to_edit``;
       a candidate whose PBS2 window falls off the genome window is skipped.

    Returns ``[]`` when no admissible, off-target-clean, rcPBS2-derivable ngRNA exists,
    so the caller can fall back to the genome-enumerated path. ``arm_side="auto"`` (the
    default) applies no side filter, preserving the behaviour of existing callers.
    """
    ngrnas = nicking_rows_to_ngrnas(nicking_rows)
    ordered = _rank_admissible(
        ngrnas,
        min_distance=ngrna_distance_min,
        max_distance=ngrna_distance_max,
        pegrna_strand=pegrna_strand,
    )
    ordered = [
        ng for ng in ordered if nick_usable_for_arm(site_s + ng.distance_to_edit, site_s, arm_side)
    ]
    if not ordered:
        return []

    scores = offtarget_scorer([ng.spacer for ng in ordered])
    out: List[InjectedNgRNA] = []
    for ng in ordered:
        score = scores.get(ng.spacer, (0.0, 10**9, True))
        if offtarget_rejection(score, max_offtargets=max_offtargets, min_ontarget=min_ontarget):
            continue  # NICK lever: skip this nick, try the next admissible candidate
        off = int(score[1])
        nick_pos = site_s + ng.distance_to_edit
        try:
            rcpbs2 = derive_rc_pbs2(genome, ngrna_nick_position=nick_pos, pbs2_length=pbs2_length)
        except ValueError:
            continue  # PBS2 window off the genome window -> not a usable nick
        out.append(
            InjectedNgRNA(
                nick_position=nick_pos,
                protospacer=ng.spacer,
                rcpbs2=rcpbs2,
                distance_to_edit=ng.distance_to_edit,
                target_strand=ng.target_strand,
                deepcas9_score=ng.deepcas9_score,
                offtargets=off,
            )
        )
    return out


def _default_context_fetcher(transcript_id: str) -> str:
    """Fetch the reference (hap1) GRCh38 region sequence for a transcript.

    Imported lazily so the module (and CI) never pulls the Ensembl/requests path unless
    a real auto-fetch is actually needed.
    """
    from pegrna_pipeline.hpc.genome_ingest import fetch_genomic_context  # noqa: PLC0415

    ctx = fetch_genomic_context(transcript_id)
    return ctx.haplotypes[0].genomic_sequence


def resolve_pridict_context(
    *,
    genome: str,
    site_s: int,
    min_flank: int,
    transcript: Optional[str] = None,
    fetcher: Optional[Callable[[str], str]] = None,
    target_flank: int = 175,
) -> Tuple[str, int]:
    """Return ``(context_sequence, context_site_s)`` with ``>= min_flank`` bp each side.

    * If the provided ``genome`` window already has ``>= min_flank`` bp on both sides of
      ``site_s``, it is used verbatim (identity return).
    * Otherwise a wider GRCh38 window is fetched (``fetcher`` or the Ensembl-backed
      default) and the provided window is located inside it by exact substring match to
      recover the offset (mapping ``site_s`` into the fetched sequence). The fetched
      context is trimmed to ~``target_flank`` bp each side (never below ``min_flank``).

    Raises :class:`ValueError` (fail-loud, actionable) when the window is too short and
    there is no transcript to fetch from, when the provided window cannot be localised
    inside the fetched context (e.g. a strand flip or an applied variant), or when the
    fetched context still lacks ``min_flank`` bp on a side.
    """
    left = site_s
    right = len(genome) - site_s
    if left >= min_flank and right >= min_flank:
        return genome, site_s

    if not transcript:
        raise ValueError(
            f"genome window too short for PRIDICT context: need >= {min_flank} bp each "
            f"side of site_s={site_s} (have left={left} right={right}) and no transcript "
            f"was given to auto-fetch a wider GRCh38 window"
        )

    fetch = fetcher or _default_context_fetcher
    wide = fetch(transcript)
    idx = wide.find(genome)
    if idx < 0:
        raise ValueError(
            f"could not localise the provided genome window inside the fetched GRCh38 "
            f"context for transcript {transcript!r} (exact substring match failed; a "
            f"strand flip or an applied variant in the window would cause this)"
        )
    ctx_site = idx + site_s
    left2 = ctx_site
    right2 = len(wide) - ctx_site
    if left2 < min_flank or right2 < min_flank:
        raise ValueError(
            f"fetched GRCh38 context for transcript {transcript!r} is still too short: "
            f"need >= {min_flank} bp each side of the insert (have left={left2} "
            f"right={right2})"
        )

    want = max(min_flank, target_flank)
    lo = max(0, ctx_site - want)
    hi = min(len(wide), ctx_site + want)
    return wide[lo:hi], ctx_site - lo
