"""Unit tests for driver helpers (Plan 3 T10): LocalShellClient,
_local_compute_ssh factory, _build_variant_scores_from_batch,
_extract_winning_guides.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from pegrna_pipeline.hpc.data_types import (
    GeneRecord,
    NgRNADesign,
    PegRNADesign,
    VariantScore,
)
from pegrna_pipeline.hpc.driver import (
    LocalShellClient,
    _build_variant_scores_from_batch,
    _extract_winning_guides,
    _local_compute_ssh,
)
from pegrna_pipeline.hpc.pa_workflow import PAResult
from pegrna_pipeline.hpc.pridict_submitter import BatchResult
from pegrna_pipeline.hpc.tjpe_workflow import TJPEResult


# ------------------------------------------------------ LocalShellClient
def test_local_shell_client_exec_returns_stdout():
    with LocalShellClient() as ssh:
        out, err, rc = ssh.exec("echo hello")
        assert out.strip() == "hello"
        assert rc == 0


def test_local_shell_client_upload_copies_file(tmp_path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("hello")
    with LocalShellClient() as ssh:
        ssh.upload(str(src), str(dst))
    assert dst.read_text() == "hello"


def test_local_shell_client_download_copies_file(tmp_path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("hello")
    with LocalShellClient() as ssh:
        ssh.download(str(src), str(dst))
    assert dst.read_text() == "hello"


def test_local_shell_client_context_manager_noop():
    """No network to connect or disconnect."""
    ssh = LocalShellClient()
    ssh.connect()
    ssh.disconnect()


# ------------------------------------------------------ _local_compute_ssh
def test_local_compute_ssh_local_mode():
    config = {"compute_mode": "local"}
    ssh = _local_compute_ssh(config)
    assert isinstance(ssh, LocalShellClient)


def test_local_compute_ssh_back_ssh_mode(monkeypatch):
    """When compute_mode=back_ssh, return a real SSHClient (unconnected)."""
    from pegrna_pipeline.hpc.ssh_client import SSHClient

    config = {
        "compute_mode": "back_ssh",
        "hostname": "cx-login",
        "username": "mj",
        "password": "pw",
    }
    ssh = _local_compute_ssh(config)
    assert isinstance(ssh, SSHClient)
    # Should not have connected (no real host reachable in tests)
    assert not ssh.is_connected


# ------------------------------------------------------ batch -> variant scores
def test_build_variant_scores_from_batch_matches_outputs(tmp_path):
    # Create two fake output dirs with csv pairs.
    out0 = tmp_path / "out_0"
    out1 = tmp_path / "out_1"
    out0.mkdir()
    out1.mkdir()
    # Minimal PRIDICT-style CSVs — parse_pegrna_csv needs the headers
    # used elsewhere in the codebase. Let's inspect pridict_parser to
    # confirm, but for now we skip those details and craft rows that
    # parse_pegrna_csv accepts (covered by its own existing tests).

    # Instead of going through parse_pegrna_csv here, the helper can be
    # monkey-patched in tests. We'll use a stub that returns canned design
    # lists.

    designs_by_prefix = {}

    from pegrna_pipeline.hpc import driver as driver_mod

    def fake_parse_peg(path):
        return [
            PegRNADesign(
                sequence_name=Path(path).stem.replace("_pegRNA_Pridict_full", ""),
                spacer="GATCGATCGATCGATCGATC",
                pbs="AAAAAAAAAA",
                rt="TTTTTTTTTT",
                pbs_length=10,
                rt_length=10,
                target_strand="Fw",
                editing_position=5,
                pe2_hek293t=0.8,
                pe2_k562=0.7,
                pegrna_full="ACGT",
            )
        ]

    def fake_parse_ng(path):
        return [
            NgRNADesign(
                spacer="CCCCCCCCCCCCCCCCCCCC",
                distance_to_edit=60,
                target_strand="Rv",
                deepcas9_score=0.9,
                pe3b=False,
                pam_disrupt=False,
            )
        ]

    # Create files matching the parser discovery scheme.
    (out0 / "v_0_pegRNA_Pridict_full.csv").write_text("stub")
    (out0 / "v_0_nicking_guides.csv").write_text("stub")
    (out1 / "v_1_pegRNA_Pridict_full.csv").write_text("stub")
    (out1 / "v_1_nicking_guides.csv").write_text("stub")

    import pegrna_pipeline.hpc.pridict_parser as pp

    orig_peg = pp.parse_pegrna_csv
    orig_ng = pp.parse_ngrna_csv
    pp.parse_pegrna_csv = fake_parse_peg
    pp.parse_ngrna_csv = fake_parse_ng
    # The driver imports these names at the module level of driver — patch there too:
    driver_mod.parse_pegrna_csv = fake_parse_peg
    driver_mod.parse_ngrna_csv = fake_parse_ng
    try:
        batch = BatchResult(
            submitted=["j0", "j1"],
            failed=[],
            output_dirs=[str(out0), str(out1)],
        )
        variants = [
            {"sequence_name": "v_0", "editseq": "x"},
            {"sequence_name": "v_1", "editseq": "y"},
        ]
        scores = _build_variant_scores_from_batch(batch, variants)
        assert len(scores) == 2
        assert all(isinstance(s, VariantScore) for s in scores)
        assert {s.sequence_name for s in scores} == {"v_0", "v_1"}
    finally:
        pp.parse_pegrna_csv = orig_peg
        pp.parse_ngrna_csv = orig_ng


# ------------------------------------------------------ _extract_winning_guides
def test_extract_winning_guides_from_tjpe_result():
    vs = VariantScore(
        sequence_name="v", codon_variant_dna="A",
        pegRNA_designs=[
            PegRNADesign(
                sequence_name="v", spacer="GATCGATCGATCGATCGATC",
                pbs="A", rt="T", pbs_length=1, rt_length=1,
                target_strand="Fw", editing_position=0,
                pe2_hek293t=0.5, pe2_k562=0.4, pegrna_full="ACGT",
            )
        ],
    )
    res = TJPEResult(
        winning_variant=vs,
        full_insert_dna="ACGT",
        tjpe_pegRNA_sequence="AAA",
        spacer="GATCGATCGATCGATCGATC",
        pbs1="AAAAAAAAAA",
        rc_pbs2="TTTTT",
        rc_pbs2_length=5,
        ngrna_spacer="CCCCCCCCCCCCCCCCCCCC",
        ngrna_distance_to_edit=60,
        ngrna_deepcas9_score=0.9,
        pegrna_pe2_hek293t=0.5,
        edit_position_in_mrna=100,
        ngrna_nick_in_mrna=160,
    )
    guides = _extract_winning_guides(res)
    assert len(guides) == 2
    ids = {g["guide_id"] for g in guides}
    assert ids == {"pegrna_spacer", "ngrna_spacer"}


def test_extract_winning_guides_from_pa_result():
    from pegrna_pipeline.hpc.pair_scorer import PegRNAPair

    up = PegRNADesign(
        sequence_name="u", spacer="AAAAAAAAAAAAAAAAAAAA",
        pbs="A", rt="T", pbs_length=1, rt_length=1,
        target_strand="Fw", editing_position=0,
        pe2_hek293t=0.5, pe2_k562=0.4, pegrna_full="ACGT",
    )
    down = PegRNADesign(
        sequence_name="d", spacer="TTTTTTTTTTTTTTTTTTTT",
        pbs="A", rt="T", pbs_length=1, rt_length=1,
        target_strand="Rv", editing_position=0,
        pe2_hek293t=0.6, pe2_k562=0.5, pegrna_full="TGCA",
    )
    pair = PegRNAPair(upstream=up, downstream=down, joint_score=1.0, nick_distance=40)
    result = PAResult(
        winning_pair=pair,
        top_pairs=[pair],
        linear_donor_dna="ACGT",
        upstream_pegRNA_sequence="U",
        downstream_pegRNA_sequence="D",
        upstream_flap_dna="A",
        downstream_flap_dna="G",
        bulk_dna="C",
        full_integrated_product="ACGT",
        upstream_pegRNA_off_target="",
        downstream_pegRNA_off_target="",
        total_pairs_evaluated=2,
        pairs_after_compat_filter=1,
    )
    guides = _extract_winning_guides(result)
    ids = {g["guide_id"] for g in guides}
    assert ids == {"upstream_pegrna_spacer", "downstream_pegrna_spacer"}


def test_extract_winning_guides_empty_results():
    # No winning guides.
    res = TJPEResult(
        winning_variant=None,
        full_insert_dna="ACGT",
        tjpe_pegRNA_sequence="",
        spacer="",
        pbs1="",
        rc_pbs2="",
        rc_pbs2_length=0,
        ngrna_spacer="",
        ngrna_distance_to_edit=0,
        ngrna_deepcas9_score=0.0,
        pegrna_pe2_hek293t=0.0,
        edit_position_in_mrna=0,
        ngrna_nick_in_mrna=0,
    )
    assert _extract_winning_guides(res) == []
