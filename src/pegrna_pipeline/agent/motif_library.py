"""Motif library: bundled + user override, case-insensitive lookup, fuzzy suggestions."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Iterable

import yaml


class MotifNotFoundError(KeyError):
    """Raised when a motif name is not present in the library."""


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
            resources.files("pegrna_pipeline.data")
            .joinpath("motif_library.yaml")
            .read_text()
        )
        motifs = [Motif.from_mapping(m) for m in data.get("motifs", [])]
        return cls(motifs)

    @classmethod
    def from_path(cls, path: Path) -> "MotifLibrary":
        data = yaml.safe_load(path.read_text())
        motifs = [Motif.from_mapping(m) for m in data.get("motifs", [])]
        return cls(motifs)

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
