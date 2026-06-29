"""Literature-derived gate thresholds: values + that the gates consume them."""

from pegrna_pipeline.hpc.gate_thresholds import THRESHOLDS


def test_splice_delta_is_high_recall_cutoff():
    # Jaganathan 2019 Cell: 0.2 high-recall (safety gate -> maximise sensitivity).
    assert THRESHOLDS.splice_max_delta == 0.2


def test_alphafold_plddt_confident_band():
    # Jumper 2021 / AFDB: 70-90 "confident".
    assert THRESHOLDS.plddt_min == 70.0


def test_alphafold_rsa_buried_levy2010():
    # Levy 2010 J Mol Biol: interior/core rASA < 0.25.
    assert THRESHOLDS.rsa_buried == 0.25


def test_microhomology_run_at_meps():
    # RecA MEPS ~8 bp / MMEJ efficient onset ~6-8 bp.
    assert THRESHOLDS.microhomology_max_run == 8


def test_docking_bands_from_af_multimer():
    # Evans 2022 / Abramson 2024: ipTM>0.8 high, 0.6-0.8 medium, <0.6 unreliable.
    assert THRESHOLDS.docking_min_iptm == 0.6
    assert THRESHOLDS.docking_max_drop == 0.2


def test_rinalmo_pairing_decision_boundary():
    assert THRESHOLDS.rinalmo_pair_prob == 0.5


def test_each_threshold_has_a_citation():
    for key, cite in THRESHOLDS.citations.items():
        assert isinstance(cite, str) and len(cite) > 10, key


def test_gates_consume_thresholds():
    # The gate functions default to the literature values.
    import inspect

    from pegrna_pipeline.hpc.complex_dock import passes_docking_gate
    from pegrna_pipeline.hpc.splice_screen import passes_splice_gate
    from pegrna_pipeline.hpc.structure_importance import residue_blocked

    assert inspect.signature(passes_splice_gate).parameters["max_delta"].default == 0.2
    assert inspect.signature(passes_docking_gate).parameters["max_drop"].default == 0.2
    assert inspect.signature(residue_blocked).parameters["rsa_buried"].default == 0.25
