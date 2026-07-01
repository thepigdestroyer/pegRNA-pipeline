"""Pure microhomology metrics: longest contiguous identity between two strings.

Used by ``recode_arm.py`` to drive the deletion-replacement-arm recode toward the
minimum synonymous edits that break microhomology with the displaced genomic
strand, and by the full-pegRNA self-structure check.
"""

from __future__ import annotations

__all__ = ["longest_common_run", "max_identity_run"]


def longest_common_run(a: str, b: str) -> int:
    """Length of the longest contiguous substring common to ``a`` and ``b``.

    Classic DP (``O(len(a)*len(b))``). Sequences here are short (<=~120 nt) so
    this is fine; swap for a suffix automaton only if profiling ever demands it.
    Case-insensitive.
    """
    a, b = a.upper(), b.upper()
    if not a or not b:
        return 0
    prev = [0] * (len(b) + 1)
    best = 0
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        ai = a[i - 1]
        for j in range(1, len(b) + 1):
            if ai == b[j - 1]:
                run = prev[j - 1] + 1
                cur[j] = run
                if run > best:
                    best = run
        prev = cur
    return best


def max_identity_run(a: str, b: str) -> int:
    """Alias kept explicit for call-site readability (same as longest_common_run)."""
    return longest_common_run(a, b)
