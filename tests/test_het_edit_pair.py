"""ADR 0003: design_recoded_with_het_edit emits the antisense primary pair and,
when a discriminating phased het exists, a sense-window het-edit pair (byte-exact
3'UTR replacement). Gated on the het: no discriminating het -> no second pair."""

from pegrna_pipeline.hpc.recoded_driver import design_recoded_with_het_edit

_GENOME = "ACGATGCATGGCATGACTGA" * 16  # 320 nt; known-passing fixture
_SPACER = "ACGTACGTACGTACGTACGT"

_PRIMARY = dict(  # antisense recode-into-CDS: arm all-CDS
    spacer=_SPACER,
    rtt_aa="MKP",
    cargo="ATGAAACCCGGG",
    genome=_GENOME,
    site_s=150,
    cds_mask=[True] * 320,
    pbs_lengths=(15,),
    pegrna_strand="Rv",
    max_run=8,
)
_HET_WINDOW = dict(  # sense window: CDS 5' of S, 3'UTR (non-coding) 3' of S
    spacer=_SPACER,
    rtt_aa="MKP",
    cargo="ATGAAACCCGGG",
    genome=_GENOME,
    site_s=150,
    cds_mask=[i < 150 for i in range(320)],
    pbs_lengths=(15,),
    pegrna_strand="Rv",
    max_run=8,
)


def test_emits_both_pairs_when_discriminating_het_present():
    het = [{"element": "pegrna_pam", "offset": 1, "haplotype": 1, "phased": True}]
    out = design_recoded_with_het_edit(
        primary=_PRIMARY, het_edit_window=_HET_WINDOW, het_snps=het, mode="fake"
    )
    assert out.primary.result.pegRNA_sequence
    assert out.primary.result.ngrna_protospacer  # full pair (pegRNA + ngRNA)
    assert out.het_edit is not None
    assert out.het_edit.result.pegRNA_sequence
    assert out.het_edit.result.ngrna_protospacer
    assert out.het_edit.result.recode_edits == 0  # 3'UTR byte-exact (exact replacement)
    assert "het-edit" in out.het_edit_note.lower()


def test_no_het_edit_pair_without_discriminating_het():
    out = design_recoded_with_het_edit(
        primary=_PRIMARY, het_edit_window=_HET_WINDOW, het_snps=[], mode="fake"
    )
    assert out.primary.result.pegRNA_sequence
    assert out.het_edit is None
    assert "no discriminating" in out.het_edit_note.lower()


def test_primary_only_when_no_het_window_given():
    out = design_recoded_with_het_edit(primary=_PRIMARY, mode="fake")
    assert out.primary.result.pegRNA_sequence
    assert out.het_edit is None
