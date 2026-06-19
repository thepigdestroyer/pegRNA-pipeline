"""Scarless tjPE: restored-segment construction + INV-0 verification.

This module implements the zero-endogenous-deletion ("scarless") tjPE design
(see docs/superpowers/plans/2026-06-09-scarless-tjpe.md and the vault notes
305/306). It is deliberately free of cluster / Biopython-heavy imports beyond
``Bio.Seq`` so it can be unit-tested in isolation.

Key pieces:

* :func:`build_restored_segment` — the exact A→B genomic span (revcomp for Rv).
* :func:`predict_edited_allele`  — reconstruct the predicted *sense* allele.
* :func:`assert_inv0`            — HARD GATE: zero endogenous deletion / ≤1 silent base.
* :func:`silent_pam_break`       — a single synonymous base change that abolishes a PAM.

Risk note (adversarial review, accepted under the A-only policy): the restored
segment is 100% homologous to the strand the flap must displace — the exact
microhomology condition Zheng 2023 (Fah) engineered away. Callers should surface
``microhomology_flag`` accordingly. See spec §10 risk register.
"""

from __future__ import annotations

from typing import Callable, List, Optional, Tuple, Union

from Bio.Seq import Seq

__all__ = [
    "build_restored_segment",
    "predict_edited_allele",
    "assert_inv0",
    "silent_pam_break",
]


def build_restored_segment(
    genome: str,
    *,
    a_nick: int,
    b_nick: int,
    pegrna_strand: str,
) -> str:
    """Return the exact A→B genomic span, reverse-complemented for Rv pegRNAs.

    ``a_nick`` / ``b_nick`` are 0-indexed positions on the sense ``genome``; the
    restored span is ``genome[min:max]``. A forward pegRNA templates the span in
    sense orientation; a reverse-strand pegRNA templates its reverse complement
    so that the *sense* strand reads the intended sequence after synthesis.
    """
    lo, hi = sorted((a_nick, b_nick))
    if lo < 0 or hi > len(genome):
        raise ValueError(f"A→B span [{lo},{hi}) out of bounds for genome length {len(genome)}")
    span = genome[lo:hi]
    if pegrna_strand == "Fw":
        return span
    if pegrna_strand == "Rv":
        return str(Seq(span).reverse_complement())
    raise ValueError(f"pegrna_strand must be 'Fw' or 'Rv'; got {pegrna_strand!r}")


def predict_edited_allele(
    genome: str,
    *,
    site_s: int,
    cargo: str,
    a_nick: int,  # noqa: ARG001 — part of the design contract / documentation
    b_nick: int,  # noqa: ARG001
    pegrna_strand: str,  # noqa: ARG001
    silent_pam_edit: Optional[dict],
) -> str:
    """Reconstruct the predicted post-edit **sense** allele.

    Scarless ⇒ reference with ``cargo`` inserted at ``site_s`` and the A→B span
    restored (so nothing endogenous is deleted). The sense allele is identical
    regardless of which strand was synthesised — the reverse-complement handling
    lives only in :func:`build_restored_segment` / cargo templating. If
    ``silent_pam_edit`` is given it applies exactly one base change at ``pos``
    (sense coordinates, pre-insertion); the coordinate is shifted by the cargo
    length when it lies at/after ``site_s``.
    """
    edited = genome[:site_s] + cargo + genome[site_s:]
    if silent_pam_edit is not None:
        pos = silent_pam_edit["pos"]
        shifted = pos + len(cargo) if pos >= site_s else pos
        edited = edited[:shifted] + silent_pam_edit["alt"] + edited[shifted + 1 :]
    return edited


def assert_inv0(
    *,
    reference: str,
    edited: str,
    site_s: int,
    cargo: str,
    allowed_silent: Union[dict, List[dict], None],
    recode_span: Optional[Tuple[int, int]] = None,
    protein_check: Optional[Callable[[str, str], bool]] = None,
    locked_span: Optional[Tuple[int, int]] = None,
    locked_expected: Optional[str] = None,
) -> None:
    """HARD GATE — zero endogenous deletion + bounded, justified changes.

    ``edited`` must equal ``reference`` with ``cargo`` inserted at ``site_s``,
    except for changes that are explicitly permitted:

    * Inside ``recode_span`` (CDS): any number of changes, provided they are
      protein-silent (verified by the injected ``protein_check(ref, edited)``).
    * Inside ``locked_span`` (non-coding; decision 15): **zero** changes — the
      slice must equal ``locked_expected`` (the targeted allele's bytes) byte for
      byte. No silent edits allowed in non-coding.
    * Outside both spans: at most ``len(allowed_silent)`` single-base silent PAM
      edits (B2-ii: ≤2). A 3rd outside change ⇒ fail.

    No endogenous insertion/deletion is tolerated (length must match exactly).

    **Backward compatibility:** with ``recode_span=None`` and ``allowed_silent``
    either ``None`` or a single ``dict``, the legacy ≤1-silent-base behaviour is
    preserved byte-for-byte. ``allowed_silent`` may be a ``dict`` (legacy) or a
    ``list[dict]`` (B2-ii).
    """
    expected = reference[:site_s] + cargo + reference[site_s:]
    if len(expected) != len(edited):
        raise ValueError(
            "INV-0: length mismatch — endogenous insertion/deletion detected "
            f"(expected {len(expected)} nt, got {len(edited)} nt)"
        )

    # Normalise allowed_silent to a list and count the budget.
    if allowed_silent is None:
        silent_budget = 0
    elif isinstance(allowed_silent, dict):
        silent_budget = 1
    else:
        silent_budget = len(allowed_silent)
    if silent_budget > 2:
        raise ValueError(f"INV-0: {silent_budget} silent PAM edits declared; hard cap is 2 (B2-ii)")

    # ---- legacy path (no recode span): preserve exact historical semantics ----
    if recode_span is None:
        if edited == expected:
            return
        diffs = [i for i, (a, b) in enumerate(zip(expected, edited)) if a != b]
        if silent_budget == 0:
            raise ValueError(f"INV-0: {len(diffs)} endogenous change(s) but no silent edit allowed")
        if len(diffs) != silent_budget:
            raise ValueError(
                f"INV-0: {len(diffs)} endogenous changes; at most {silent_budget} "
                "silent PAM base(s) allowed"
            )
        return

    # ---- span-aware path ----
    # A span coordinate is in *reference* space; cargo is inserted at site_s.
    # The span's start shifts only if strictly past S; the (exclusive) end shifts
    # only if strictly past S, so a span ending exactly at S keeps the cargo OUT.
    def _shift_lo(pos: int) -> int:
        return pos + len(cargo) if pos >= site_s else pos

    def _shift_hi(pos: int) -> int:
        return pos + len(cargo) if pos > site_s else pos

    rs = (_shift_lo(recode_span[0]), _shift_hi(recode_span[1]))
    ls = (_shift_lo(locked_span[0]), _shift_hi(locked_span[1])) if locked_span else None

    # Locked (non-coding) span must be byte-exact to the targeted allele.
    if ls is not None:
        if locked_expected is None:
            raise ValueError("INV-0: locked_span given without locked_expected")
        observed = edited[ls[0] : ls[1]]
        if observed != locked_expected:
            raise ValueError(
                "INV-0: non-coding locked span altered — must be byte-exact to the "
                f"targeted allele (expected {locked_expected!r}, got {observed!r})"
            )

    diffs = [i for i, (a, b) in enumerate(zip(expected, edited)) if a != b]
    outside = []
    for i in diffs:
        if rs[0] <= i < rs[1]:
            continue  # inside recode span — allowed if protein-silent (checked below)
        if ls is not None and ls[0] <= i < ls[1]:
            continue  # handled by the byte-exact check above
        outside.append(i)

    if len(outside) > silent_budget:
        raise ValueError(
            f"INV-0: {len(outside)} change(s) outside recode_span/cargo/locked_span "
            f"(max {silent_budget} silent PAM base(s) allowed)"
        )

    # Recode-span changes must be protein-silent.
    if any(rs[0] <= i < rs[1] for i in diffs):
        if protein_check is not None and not protein_check(reference, edited):
            raise ValueError("INV-0: recode_span changes are not protein-silent")
    return


def silent_pam_break(
    *,
    codon_context: str,
    pam_start: int,
    pam: str,  # noqa: ARG001 — documented contract; the two G positions are derived from pam_start
    frame: int,
) -> Optional[dict]:
    """Find a single base change inside the PAM that abolishes an NGG/NAG motif
    **and** is translationally silent.

    Only the two positions after the leading ``N`` of an ``NGG`` PAM are
    considered (mutating ``N`` does not reliably abolish recognition). Returns a
    dict describing the change, or ``None`` if no synonymous PAM-breaking base
    exists.
    """
    ctx = list(codon_context.upper())
    for offset in (1, 2):  # the two PAM positions after N
        pos = pam_start + offset
        if pos >= len(ctx):
            continue
        ref = ctx[pos]
        for alt in "ACGT":
            if alt == ref:
                continue
            trial = ctx.copy()
            trial[pos] = alt
            cstart = ((pos - frame) // 3) * 3 + frame
            before = "".join(ctx[cstart : cstart + 3])
            after = "".join(trial[cstart : cstart + 3])
            if len(before) < 3 or len(after) < 3:
                continue
            aa_before = str(Seq(before).translate())
            aa_after = str(Seq(after).translate())
            if aa_before == aa_after:
                return {
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "codon_before": before,
                    "codon_after": after,
                    "aa_before": aa_before,
                    "aa_after": aa_after,
                }
    return None
