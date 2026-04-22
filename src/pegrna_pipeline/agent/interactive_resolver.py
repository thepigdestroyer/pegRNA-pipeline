"""Build an interactive resolver callable for unknown motif names.

Presents a numbered menu containing fuzzy suggestions followed by three
"always-available" choices: raw aa input, raw dna input, and abort.
Returns a :class:`RawComponent` or raises :class:`AbortedByUser`.
"""

from __future__ import annotations

from typing import Callable

from pegrna_pipeline.agent.construct_parser import (
    RawComponent,
    RawComponentKind,
)
from pegrna_pipeline.agent.motif_library import (
    MotifLibrary,
    fuzzy_suggestions,
)


class AbortedByUser(RuntimeError):
    """Raised when the user selects ``abort`` from the interactive menu."""


PromptFn = Callable[[str], str]


def _default_prompt(msg: str) -> str:  # pragma: no cover - thin wrapper
    return input(msg)


def build_interactive_resolver(
    library: MotifLibrary,
    *,
    prompt: PromptFn = _default_prompt,
    max_suggestions: int = 5,
) -> Callable[[str], RawComponent]:
    """Return a resolver(name) -> RawComponent callable bound to ``library``.

    The returned callable is called for each unknown motif encountered during
    construct resolution. It prints a numbered menu to stdout and reads the
    user's choice via ``prompt``.
    """

    def resolve(unknown: str) -> RawComponent:
        suggestions = fuzzy_suggestions(unknown, library)[:max_suggestions]
        while True:
            _render_menu(unknown, suggestions)
            raw = prompt("Choose: ").strip()

            # Keyword shortcuts.
            lowered = raw.lower()
            if lowered in {"abort", "a"}:
                raise AbortedByUser(f"user aborted while resolving {unknown!r}")
            if lowered in {"aa", "aa:"}:
                payload = prompt("Enter amino acid sequence (single-letter): ").strip().upper()
                return RawComponent(RawComponentKind.RAW_AA, payload)
            if lowered in {"dna", "dna:"}:
                payload = prompt("Enter DNA sequence (ACGT): ").strip().upper()
                return RawComponent(RawComponentKind.RAW_DNA, payload)

            # Numeric choice maps to suggestion index.
            if raw.isdigit():
                idx = int(raw)
                if 1 <= idx <= len(suggestions):
                    return RawComponent(
                        RawComponentKind.MOTIF, suggestions[idx - 1]
                    )
                # Indices beyond suggestions correspond to aa / dna / abort.
                extra = idx - len(suggestions)
                if extra == 1:
                    payload = prompt("Enter amino acid sequence (single-letter): ").strip().upper()
                    return RawComponent(RawComponentKind.RAW_AA, payload)
                if extra == 2:
                    payload = prompt("Enter DNA sequence (ACGT): ").strip().upper()
                    return RawComponent(RawComponentKind.RAW_DNA, payload)
                if extra == 3:
                    raise AbortedByUser(
                        f"user aborted while resolving {unknown!r}"
                    )

            print(f"Unrecognised choice {raw!r}; please try again.")

    return resolve


def _render_menu(unknown: str, suggestions: list[str]) -> None:
    print(f"Unknown motif {unknown!r}. Did you mean:")
    next_index = 1
    for suggestion in suggestions:
        print(f"  [{next_index}] {suggestion}")
        next_index += 1
    print(f"  [{next_index}] aa: enter raw amino acid sequence")
    print(f"  [{next_index + 1}] dna: enter raw DNA sequence")
    print(f"  [{next_index + 2}] abort")
