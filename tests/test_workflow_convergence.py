"""Stage H: one canonical tj-PE design workflow; the legacy orchestrator is deprecated.

``recoded_tjpe_workflow.run_recoded_tjpe`` is the single canonical design workflow (the
CLI + web use it via ``recoded_driver.design_recoded``). ``tjpe_workflow.run_tjpe_workflow``
is the legacy variant-ranking orchestrator (the HPC ``driver`` path), retained only for the
shared primitives + backward compatibility and now deprecated. These tests pin that
convergence and guard against the two paths drifting apart on pegRNA construction (the
divergence that produced the historical rcPBS2/PBS2 bug the unified plan set out to kill).
"""

from __future__ import annotations

import inspect

import pytest

from pegrna_pipeline.hpc import recoded_driver, recoded_tjpe_workflow, tjpe_workflow


def test_run_tjpe_workflow_is_deprecated_legacy_orchestrator():
    # The deprecation fires before the input-source contract check, so a no-source call
    # both warns (legacy) and then raises the contract ValueError.
    with pytest.warns(DeprecationWarning, match="canonical"), pytest.raises(ValueError):
        tjpe_workflow.run_tjpe_workflow(gene=None, insert_full_dna="ATG")


def test_recoded_is_the_canonical_design_entrypoint():
    # design_recoded (the CLI/web entrypoint) drives the canonical run_recoded_tjpe.
    assert recoded_driver.run_recoded_tjpe is recoded_tjpe_workflow.run_recoded_tjpe
    assert "run_recoded_tjpe" in inspect.getsource(recoded_driver.design_recoded)


def test_both_workflows_share_the_authoritative_assembler_no_drift():
    # Both the canonical and the legacy path construct the pegRNA through the SAME
    # authoritative assembler (assemble_tj_pegrna_v2, Spacer.Scaffold.rcPBS2.RTT.PBS1.PolyU),
    # so they cannot silently diverge on construction again (ADR 0014/0020/0021).
    legacy_src = inspect.getsource(tjpe_workflow.run_tjpe_workflow)
    # run_recoded_tjpe is the thin INSERT-lever wrapper; the assembly lives in the core.
    canonical_src = inspect.getsource(recoded_tjpe_workflow._run_recoded_tjpe_single)
    assert "assemble_tj_pegrna_v2" in legacy_src
    assert "assemble_tj_pegrna_v2" in canonical_src
