"""Tests for :func:`retry_with_doubled_resources` (Plan 3 T3)."""

from __future__ import annotations

import pytest

from pegrna_pipeline.hpc.pbs import (
    MEM_CAP_GB,
    PBSJobSpec,
    WALLTIME_CAP_HOURS,
    retry_with_doubled_resources,
)


def _spec(**over):
    base = dict(
        name="x",
        ncpus=2,
        mem_gb=8,
        walltime_hours=4,
        log_path="/x.log",
        body="true\n",
    )
    base.update(over)
    return PBSJobSpec(**base)


def test_walltime_doubled():
    out = retry_with_doubled_resources(_spec(walltime_hours=4), reason="walltime")
    assert out.walltime_hours == 8
    # Mem unchanged
    assert out.mem_gb == 8


def test_walltime_capped():
    out = retry_with_doubled_resources(_spec(walltime_hours=20), reason="walltime")
    assert out.walltime_hours == WALLTIME_CAP_HOURS
    assert out.walltime_hours == 24


def test_oom_doubled():
    out = retry_with_doubled_resources(_spec(mem_gb=16), reason="oom")
    assert out.mem_gb == 32
    assert out.walltime_hours == 4


def test_oom_capped():
    out = retry_with_doubled_resources(_spec(mem_gb=48), reason="oom")
    assert out.mem_gb == MEM_CAP_GB
    assert out.mem_gb == 64


def test_unknown_reason_raises():
    with pytest.raises(ValueError):
        retry_with_doubled_resources(_spec(), reason="disk")  # type: ignore[arg-type]


def test_spec_is_immutable():
    """retry returns a new spec — original must be untouched."""
    original = _spec(mem_gb=8)
    out = retry_with_doubled_resources(original, reason="oom")
    assert original.mem_gb == 8
    assert out.mem_gb == 16
