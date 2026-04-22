"""Motif library: bundled + user override, case-insensitive lookup, fuzzy suggestions."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Iterable

import yaml


class MotifNotFoundError(KeyError):
    """Raised when a motif name is not present in the library."""


def _levenshtein(a: str, b: str) -> int:
    """Classic O(len(a)*len(b)) edit distance with a rolling row."""
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        curr = [i] + [0] * len(b)
        for j, cb in enumerate(b, start=1):
            cost = 0 if ca == cb else 1
            curr[j] = min(
                curr[j - 1] + 1,  # insertion
                prev[j] + 1,  # deletion
                prev[j - 1] + cost,  # substitution
            )
        prev = curr
    return prev[-1]


def fuzzy_suggestions(query: str, library: "MotifLibrary", max_distance: int = 2) -> list[str]:
    """Return motif names within Levenshtein distance <= max_distance of query.

    Comparison is case-insensitive; returned names preserve library casing and
    are ordered by (distance, name). Duplicate lowercased names are collapsed
    so aliases such as ZeocinR / ZEOCINr aren't double-counted.
    """
    q = query.lower()
    scored: list[tuple[int, str]] = []
    seen_lower: set[str] = set()
    for motif in library:
        lower = motif.name.lower()
        if lower in seen_lower:
            continue
        seen_lower.add(lower)
        d = _levenshtein(q, lower)
        if d <= max_distance:
            scored.append((d, motif.name))
    scored.sort(key=lambda t: (t[0], t[1].lower()))
    return [name for _, name in scored]


@dataclass(frozen=True)
class Motif:
    """A single library entry."""

    name: str
    category: str
    dna: str
    aa: str
    description: str = ""

    @classmethod
    def from_mapping(cls, data: dict) -> "Motif":
        return cls(
            name=data["name"],
            category=data.get("category", "unknown"),
            dna=data["dna"].upper(),
            aa=data.get("aa", ""),
            description=data.get("description", ""),
        )


class MotifLibrary:
    """Collection of motifs with case-insensitive lookup.

    Motif names that differ only in casing are treated as distinct entries in
    ``list_names`` (so library size counts them separately), but all of them
    are reachable case-insensitively via ``lookup``. The last registration
    wins for collisions, so user overrides replace bundled motifs.
    """

    def __init__(self, motifs: Iterable[Motif]) -> None:
        self._entries: list[Motif] = []
        self._by_lower: dict[str, Motif] = {}
        for m in motifs:
            # Preserve the original distinct entries (ZeocinR / ZEOCINr count as 2).
            self._entries.append(m)
            # Later entries with the same lowercased name override earlier ones
            # for lookup purposes so user overrides can replace bundled motifs.
            self._by_lower[m.name.lower()] = m

    # ------------------------------------------------------------------ I/O

    @classmethod
    def load_bundled(cls) -> "MotifLibrary":
        """Load the packaged motif_library.yaml."""
        data = yaml.safe_load(
            resources.files("pegrna_pipeline.data").joinpath("motif_library.yaml").read_text()
        )
        motifs = [Motif.from_mapping(m) for m in data.get("motifs", [])]
        return cls(motifs)

    @classmethod
    def from_path(cls, path: Path) -> "MotifLibrary":
        data = yaml.safe_load(path.read_text())
        motifs = [Motif.from_mapping(m) for m in data.get("motifs", [])]
        return cls(motifs)

    @classmethod
    def _user_override_path(cls) -> Path:
        return Path.home() / ".pegrna" / "motifs.yaml"

    @classmethod
    def load_user_default(cls) -> "MotifLibrary":
        """Load ~/.pegrna/motifs.yaml, or return an empty library if missing/empty."""
        path = cls._user_override_path()
        if not path.is_file():
            return cls([])
        raw = path.read_text().strip()
        if not raw:
            return cls([])
        data = yaml.safe_load(raw) or {}
        motifs = [Motif.from_mapping(m) for m in data.get("motifs", [])]
        return cls(motifs)

    @classmethod
    def load_full(cls) -> "MotifLibrary":
        """Bundled library plus user overrides (user wins on name collisions)."""
        bundled = cls.load_bundled()
        user = cls.load_user_default()
        # User entries come last so they override bundled entries with the same
        # lowercased name (see _by_lower construction).
        return cls([*bundled, *user])

    # -------------------------------------------------------------- Queries

    def lookup(self, name: str) -> Motif:
        try:
            return self._by_lower[name.lower()]
        except KeyError as exc:
            raise MotifNotFoundError(name) from exc

    def __contains__(self, name: str) -> bool:
        return name.lower() in self._by_lower

    def list_names(self) -> list[str]:
        return sorted((m.name for m in self._entries), key=str.lower)

    def __iter__(self):
        return iter(self._entries)

    def __len__(self) -> int:
        return len(self._entries)
