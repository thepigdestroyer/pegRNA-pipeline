"""A design report must never pass a stand-in gate verdict off as a real one.

`api` mode with no NVIDIA key swaps the Boltz-2 docking gate for `DockFake`, and a
run with no genome FASTA swaps cas-offinder for a clean-pass stub. Before this was
disclosed, the resulting report was byte-identical to one where both gates really
ran -- the "confident wrong answer" failure mode that gates exist to prevent. DNA
is ordered off these reports, so the disclosure is a safety contract, not cosmetics.
"""

import json
from pathlib import Path

from pegrna_pipeline.hpc.report_writer import render_recoded_report
from pegrna_pipeline.hpc.scorer_factory import make_scorers


def test_api_bundle_without_key_declares_docking_stub(monkeypatch):
    monkeypatch.delenv("NVIDIA_API_KEY", raising=False)
    scorers = make_scorers(mode="api", uniprot_acc="P49711")
    assert "docking" in scorers.stub_gates
    # No genome FASTA supplied -> off-target is a clean-pass stand-in too.
    assert "offtarget" in scorers.stub_gates


def test_api_bundle_with_key_does_not_declare_docking_stub(monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "test-key-not-used-offline")
    scorers = make_scorers(mode="api", uniprot_acc="P49711")
    assert "docking" not in scorers.stub_gates


def test_fake_bundle_declares_every_gate_as_stand_in():
    scorers = make_scorers(mode="fake")
    for gate in ("splice", "structure", "importance", "docking", "offtarget", "cryptic"):
        assert gate in scorers.stub_gates


def _ctcf_result():
    """A genuine RecodedTJPEResult from the shipped CTCF example (fake mode)."""
    from pegrna_pipeline.hpc.recoded_driver import design_recoded

    cfg = json.loads(Path("examples/ctcf_cterm/config.json").read_text())
    return design_recoded(
        spacer=cfg["spacer"],
        rtt_aa=cfg["rtt_aa"],
        cargo=cfg["cargo"],
        genome=cfg["genome"],
        site_s=cfg["site_s"],
        cds_mask=cfg["cds_mask"],
        pbs_lengths=cfg["pbs_lengths"],
        ngrna_distance_min=cfg["ngrna_distance_min"],
        ngrna_distance_max=cfg["ngrna_distance_max"],
        pegrna_strand=cfg["pegrna_strand"],
        arm_side=cfg["arm_side"],
        max_run=cfg["max_run"],
        mode="fake",
    ).result


def test_report_banner_names_the_stand_in_gates():
    """The banner sits at the top of the report and names the gates."""
    text = render_recoded_report(
        _ctcf_result(),
        manifest={"scorer_mode": "api", "stand_in_gates": "docking, offtarget"},
    )
    head = text.splitlines()[:6]
    assert any("STAND-IN GATES IN USE" in line for line in head)
    assert any("docking, offtarget" in line for line in head)


def test_report_has_no_banner_when_every_gate_is_real():
    text = render_recoded_report(_ctcf_result(), manifest={"scorer_mode": "real"})
    assert "STAND-IN" not in text


def test_gates_json_carries_stand_in_gates(tmp_path):
    """The web UI reads gates.json, so the disclosure has to survive into it."""
    from pegrna_pipeline.hpc.recoded_driver import design_recoded

    cfg = json.loads(Path("examples/ctcf_cterm/config.json").read_text())
    design_recoded(
        spacer=cfg["spacer"],
        rtt_aa=cfg["rtt_aa"],
        cargo=cfg["cargo"],
        genome=cfg["genome"],
        site_s=cfg["site_s"],
        cds_mask=cfg["cds_mask"],
        pbs_lengths=cfg["pbs_lengths"],
        ngrna_distance_min=cfg["ngrna_distance_min"],
        ngrna_distance_max=cfg["ngrna_distance_max"],
        pegrna_strand=cfg["pegrna_strand"],
        arm_side=cfg["arm_side"],
        max_run=cfg["max_run"],
        mode="fake",
        output_dir=str(tmp_path),
    )
    gates = json.loads((tmp_path / "recoded_tjpe_gates.json").read_text())
    assert gates["scorer_mode"] == "fake"
    assert "docking" in gates["stand_in_gates"]
