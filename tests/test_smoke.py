"""Smoke tests to verify package installs and imports."""

import importlib


def test_top_level_import():
    module = importlib.import_module("pegrna_pipeline")
    assert module is not None


def test_subpackage_imports():
    assert importlib.import_module("pegrna_pipeline.agent") is not None
    assert importlib.import_module("pegrna_pipeline.hpc") is not None
    assert importlib.import_module("pegrna_pipeline.data") is not None


def test_version_attribute():
    module = importlib.import_module("pegrna_pipeline")
    assert hasattr(module, "__version__")
    assert module.__version__ == "1.0.0"
