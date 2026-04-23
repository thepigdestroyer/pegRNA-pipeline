"""Chunking helpers for splitting work across PBS jobs."""

from __future__ import annotations

from typing import Iterable, Iterator, List, TypeVar

T = TypeVar("T")


def partition(items: Iterable[T], *, max_per_chunk: int) -> Iterator[List[T]]:
    """Yield successive slices of ``items`` of size up to ``max_per_chunk``.

    The last chunk may be shorter if ``len(items)`` is not a multiple of
    the chunk size. ``max_per_chunk`` must be a positive integer.
    """
    if max_per_chunk <= 0:
        raise ValueError(f"max_per_chunk must be > 0, got {max_per_chunk}")
    buf: List[T] = []
    for item in items:
        buf.append(item)
        if len(buf) >= max_per_chunk:
            yield buf
            buf = []
    if buf:
        yield buf
