"""Tests for the ViennaRNA structure wrapper. Skipped if ViennaRNA absent."""

from __future__ import annotations

import pytest

RNA = pytest.importorskip("RNA")  # ViennaRNA python bindings

from pegrna_pipeline.hpc.rnafold import fold_metrics, passes_structure  # noqa: E402


def test_fold_metrics_returns_mfe_and_accessibility():
    m = fold_metrics("GGGGAAAACCCC", regions={"pbs": (8, 12)})
    assert "mfe" in m and isinstance(m["mfe"], float)
    assert "norm_mfe" in m
    assert 0.0 <= m["accessibility"]["pbs"] <= 1.0


def test_passes_structure_returns_bool():
    assert isinstance(passes_structure("AUAUAUAUAUAU", regions={}, max_norm_mfe=0.5), bool)


def test_passes_structure_rejects_strong_hairpin():
    # A perfect self-complementary stem should be flagged as too structured.
    strong = "GGGGGGGGGGCCCCCCCCCC"
    assert passes_structure(strong, regions={}, max_norm_mfe=-0.01) is False
