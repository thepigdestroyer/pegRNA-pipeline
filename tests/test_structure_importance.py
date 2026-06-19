"""Task 9: AlphaFold structural-importance hard gate (rules pure; AFDB/DSSP deferred)."""

from pegrna_pipeline.hpc.structure_importance import (
    ResidueMetrics,
    indel_blocked,
    residue_blocked,
    span_blocked,
)


def _m(**kw):
    base = dict(
        plddt=50.0, rsa=0.5, contact_number=4, is_catalytic=False,
        is_disulfide=False, is_positive_phi_gly=False, is_structural_pro=False,
    )
    base.update(kw)
    return ResidueMetrics(**base)


def test_buried_ordered_residue_blocks():
    assert residue_blocked(_m(plddt=85, rsa=0.10, contact_number=20)) is True


def test_exposed_loop_residue_passes():
    assert residue_blocked(_m(plddt=55, rsa=0.45, contact_number=4)) is False


def test_catalytic_always_blocks():
    assert residue_blocked(_m(plddt=40, rsa=0.6, is_catalytic=True)) is True


def test_disulfide_and_structural_residues_block():
    assert residue_blocked(_m(is_disulfide=True)) is True
    assert residue_blocked(_m(is_positive_phi_gly=True)) is True
    assert residue_blocked(_m(is_structural_pro=True)) is True


def test_high_plddt_alone_does_not_block():
    # pLDDT is a gate-ENABLER, not a standalone block (must pair with burial)
    assert residue_blocked(_m(plddt=95, rsa=0.5, contact_number=4)) is False


def test_high_contact_number_blocks():
    assert residue_blocked(_m(contact_number=20)) is True


def test_span_blocked_is_or_over_residues():
    scorer = {0: _m(rsa=0.5), 1: _m(plddt=90, rsa=0.05), 2: _m(rsa=0.5)}

    class S:
        def metrics(self, i):
            return scorer[i]

    assert span_blocked(S(), [0, 2]) is False
    assert span_blocked(S(), [0, 1, 2]) is True  # residue 1 is buried+ordered


def test_indel_in_secondary_structure_blocks():
    # decision/AlphaFold-gate: in-frame indel landing in helix/strand is CRITICAL
    assert indel_blocked("H") is True   # alpha helix
    assert indel_blocked("E") is True   # beta strand
    assert indel_blocked("G") is True   # 3-10 helix
    assert indel_blocked("C") is False  # coil/loop tolerated
    assert indel_blocked("T") is False  # turn tolerated
