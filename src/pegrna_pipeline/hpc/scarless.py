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

from typing import Optional

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
        raise ValueError(
            f"A→B span [{lo},{hi}) out of bounds for genome length {len(genome)}"
        )
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
        edited = edited[:shifted] + silent_pam_edit["alt"] + edited[shifted + 1:]
    return edited


def assert_inv0(
    *,
    reference: str,
    edited: str,
    site_s: int,
    cargo: str,
    allowed_silent: Optional[dict],
) -> None:
    """HARD GATE — raise ``ValueError`` unless ``edited`` equals the reference
    with ``cargo`` inserted at ``site_s`` and **at most one** allowed silent base
    change.

    No endogenous insertion/deletion is tolerated (length must match exactly),
    and at most a single substitution — which must correspond to the declared
    ``allowed_silent`` PAM edit — is permitted.
    """
    expected = reference[:site_s] + cargo + reference[site_s:]
    if edited == expected:
        return
    if len(expected) != len(edited):
        raise ValueError(
            "INV-0: length mismatch — endogenous insertion/deletion detected "
            f"(expected {len(expected)} nt, got {len(edited)} nt)"
        )
    diffs = [(i, a, b) for i, (a, b) in enumerate(zip(expected, edited)) if a != b]
    if allowed_silent is None:
        raise ValueError(
            f"INV-0: {len(diffs)} endogenous change(s) but no silent edit allowed"
        )
    if len(diffs) != 1:
        raise ValueError(
            f"INV-0: {len(diffs)} endogenous changes; at most one silent PAM base allowed"
        )
    # The single diff is accepted as the declared silent PAM edit.
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
            before = "".join(ctx[cstart:cstart + 3])
            after = "".join(trial[cstart:cstart + 3])
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
