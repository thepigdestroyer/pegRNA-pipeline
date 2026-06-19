"""Task 17: het-SNP advisory + allele-specific pairing + allele-hunt nick movement."""

from pegrna_pipeline.hpc.allele_specificity import (
    classify_snp_position,
    discrimination_strength,
    enumerate_nick_movements,
    find_allele_specific,
)


def test_position_class_pam_seed_distal():
    assert classify_snp_position(element="spacer", offset=19, length=20) == "seed"
    assert classify_snp_position(element="spacer", offset=1, length=20) == "distal"
    assert classify_snp_position(element="pegrna_pam", offset=1, length=3) == "pam"


def test_discrimination_strength_hierarchy():
    assert discrimination_strength("pam", "pegrna") > discrimination_strength("seed", "pegrna")
    assert discrimination_strength("seed", "pegrna") > discrimination_strength("pbs1", "pegrna")
    assert discrimination_strength("distal", "pegrna") == 0.0
    # ngRNA-side is valid (nonzero) but modestly down-weighted
    assert discrimination_strength("pam", "ngrna") > 0.0
    assert discrimination_strength("pam", "ngrna") < discrimination_strength("pam", "pegrna")


def test_find_allele_specific_returns_option_when_discriminating_snp_present():
    het = [{"element": "pegrna_pam", "offset": 1, "haplotype": 1, "phased": True}]
    opts = find_allele_specific(het)
    assert len(opts) >= 1
    assert opts[0].allele_specific is True
    assert opts[0].side in ("pegrna", "ngrna")


def test_find_allele_specific_empty_when_only_distal():
    het = [{"element": "spacer", "offset": 2, "haplotype": 1, "phased": True, "length": 20}]
    assert find_allele_specific(het) == []


def test_find_allele_specific_skips_unphased():
    het = [{"element": "pegrna_pam", "offset": 1, "haplotype": 1, "phased": False}]
    assert find_allele_specific(het) == []


def test_ngrna_side_carries_byproduct_note_not_scar():
    het = [{"element": "ngrna_pam", "offset": 1, "haplotype": 1, "phased": True}]
    opts = find_allele_specific(het)
    assert opts and opts[0].side == "ngrna"
    assert "byproduct" in opts[0].note.lower()
    assert "scar" not in opts[0].note.lower()


def test_nick_movement_partitions_cds_vs_noncoding_and_locks():
    genome = "ACGT" * 50  # 200 nt
    s_site = 120
    cds_mask = [i < 120 for i in range(200)]  # CDS upstream of S; 3'UTR downstream
    het = [{"pos": 130, "element": "spacer", "offset": 1, "haplotype": 1, "phased": True}]
    moves = enumerate_nick_movements(
        s_site=s_site, genome=genome, cds_mask=cds_mask, max_shift=20, het_snps=het,
    )
    assert moves, "expected at least one movement candidate"
    three = [m for m in moves if m.direction == "3prime"]
    five = [m for m in moves if m.direction == "5prime"]
    # 3' move into 3'UTR => locked (non-coding), microhomology advisory
    assert three and three[0].locked_span is not None
    assert three[0].microhomology_advisory is True
    # 5' move into CDS => recodable
    assert five and five[0].recode_span is not None
    # ranked by smallest shift first
    assert [m.shift_abs for m in moves] == sorted(m.shift_abs for m in moves)
