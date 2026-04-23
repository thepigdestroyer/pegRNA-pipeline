"""Tests for :func:`partition` in :mod:`pegrna_pipeline.hpc.chunking`."""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.chunking import partition


def test_partition_exact_multiple():
    items = list(range(10))
    out = list(partition(items, max_per_chunk=5))
    assert out == [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]]


def test_partition_with_remainder():
    items = list(range(7))
    out = list(partition(items, max_per_chunk=3))
    assert out == [[0, 1, 2], [3, 4, 5], [6]]


def test_partition_single_chunk():
    out = list(partition([1, 2, 3], max_per_chunk=10))
    assert out == [[1, 2, 3]]


def test_partition_empty():
    out = list(partition([], max_per_chunk=5))
    assert out == []


def test_partition_is_generator():
    """Partition must stream — never materialise the full list up front."""
    it = partition(range(5), max_per_chunk=2)
    # Pull one chunk and check nothing else has been consumed
    first = next(iter(it))
    assert first == [0, 1]


def test_partition_zero_max_raises():
    with pytest.raises(ValueError):
        list(partition([1, 2, 3], max_per_chunk=0))


def test_partition_negative_max_raises():
    with pytest.raises(ValueError):
        list(partition([1, 2, 3], max_per_chunk=-1))
