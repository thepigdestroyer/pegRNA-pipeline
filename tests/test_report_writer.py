"""Tests for the Plan 4 annotator + report writers."""

from __future__ import annotations

from pathlib import Path

from pegrna_pipeline.hpc.annotator import (
    ComponentSource,
    annotate_pa_result,
    annotate_tjpe_result,
)
from pegrna_pipeline.hpc.data_types import GeneRecord, PegRNADesign, VariantScore
from pegrna_pipeline.hpc.offtarget_scorer import OffTargetReport
from pegrna_pipeline.hpc.pa_workflow import PAResult
from pegrna_pipeline.hpc.pair_scorer import PegRNAPair
from pegrna_pipeline.hpc.report_writer import (
    ALL_VARIANTS_CSV_COLUMNS,
    OFFTARGET_CSV_COLUMNS,
    PA_CSV_COLUMNS,
    TJPE_CSV_COLUMNS,
    write_reports,
)
from pegrna_pipeline.hpc.tjpe_workflow import TJPEResult


def _fake_gene() -> GeneRecord:
    return GeneRecord(
        gene_name="CTCF",
        refseq_accession="NM_006565.4",
        mrna_sequence="A" * 200 + "TAA" + "G" * 200,
        cds_start=0,
        cds_end=203,
        stop_codon_start=200,
        stop_codon="TAA",
        upstream_flank="A" * 100,
        downstream_flank="TAA" + "G" * 100,
        flanking_bp=100,
    )


def _fake_pegrna(name="v1", pe2=42.3, strand="Fw"):
    return PegRNADesign(
        sequence_name=name,
        spacer="GCAGCAGCATGCATGCAGCT",
        pbs="TCACCAGCAGCAGC",
        rt="GACAAGGACGACGAT",
        pbs_length=14,
        rt_length=15,
        target_strand=strand,
        editing_position=100,
        pe2_hek293t=pe2,
        pe2_k562=35.0,
        pegrna_full="PEG_FULL",
    )


def _fake_tjpe_result() -> TJPEResult:
    variant = VariantScore(
        sequence_name="v000001",
        codon_variant_dna="GGCGGCGGCGGCAGCGACTACAAG",
    )
    variant.pegRNA_designs = [_fake_pegrna()]
    return TJPEResult(
        winning_variant=variant,
        full_insert_dna="GGCGGCGGCGGCAGCGACTACAAG" + "A" * 444,
        tjpe_pegRNA_sequence="SPACER_SCAFFOLD_INSERT_PBS1_RCPBS2",
        spacer="GCAGCAGCATGCATGCAGCT",
        pbs1="TCACCAGCAGCAGC",
        rc_pbs2="CACGTTGCATCGATC",
        rc_pbs2_length=17,
        ngrna_spacer="GCTAGCATGCAGCATGCATC",
        ngrna_distance_to_edit=-50,
        ngrna_deepcas9_score=0.89,
        pegrna_pe2_hek293t=42.3,
        edit_position_in_mrna=200,
        ngrna_nick_in_mrna=150,
        all_evaluated_variants=[variant],
    )


def _fake_pa_result() -> PAResult:
    up = _fake_pegrna(name="u1", pe2=38.0, strand="Fw")
    # downstream editing_position differs to give non-zero distance
    dn_obj = PegRNADesign(
        sequence_name="d1",
        spacer="AAACCCGGGTTTAAACCCGG",
        pbs="CCAGCAGCAGCAGC",
        rt="ATCGTCGTCCTTGT",
        pbs_length=14,
        rt_length=14,
        target_strand="Rv",
        editing_position=160,
        pe2_hek293t=40.0,
        pe2_k562=30.0,
        pegrna_full="PEG_FULL_DN",
    )
    pair = PegRNAPair(
        upstream=up,
        downstream=dn_obj,
        nick_distance=60,
        joint_score=(38.0 * 40.0) ** 0.5,
    )
    return PAResult(
        winning_pair=pair,
        top_pairs=[pair],
        linear_donor_dna="GGCGGCGGCGGCAGCGACTACAAG" + "A" * 1000 + "GTGTGCGTGCACTTCGTGGCCGAG",
        upstream_pegRNA_sequence="UP_PEGRNA",
        downstream_pegRNA_sequence="DOWN_PEGRNA",
        upstream_flap_dna="GGCGGCGGCGGCAGCGACTACAAG",
        downstream_flap_dna="GTGTGCGTGCACTTCGTGGCCGAG",
        bulk_dna="A" * 1000,
        full_integrated_product="<product>",
        upstream_pegRNA_off_target="",
        downstream_pegRNA_off_target="",
        total_pairs_evaluated=4,
        pairs_after_compat_filter=1,
    )


def _fake_offtarget():
    return [
        OffTargetReport(
            guide_id="tjpe_pegRNA",
            spacer="GCAGCAGCATGCATGCAGCT",
            deepspcas9=0.85,
            deephf=0.78,
            cfd_specificity=0.89,
            mm0_count=1,
            mm1_count=0,
            mm2_count=2,
            mm3_count=37,
        ),
    ]


# ------------------------------------------------------------- annotator


def test_annotate_tjpe_yields_expected_components():
    annotations = annotate_tjpe_result(_fake_tjpe_result(), _fake_gene())
    names = [a.name for a in annotations]
    assert names == [
        "spacer",
        "sgRNA_scaffold",
        "insert",
        "PBS1",
        "RC-PBS2",
        "ngRNA_spacer",
        "ngRNA_scaffold",
    ]
    spacer = next(a for a in annotations if a.name == "spacer")
    assert spacer.source == ComponentSource.PRIDICT2
    rc = next(a for a in annotations if a.name == "RC-PBS2")
    assert rc.source == ComponentSource.GENOMIC_CONTEXT


def test_annotate_pa_includes_donor():
    annotations = annotate_pa_result(_fake_pa_result(), _fake_gene())
    names = [a.name for a in annotations]
    assert "linear_donor_DNA" in names
    assert "upstream_flap_DNA" in names
    assert "downstream_flap_DNA" in names


# ----------------------------------------------------- report writers


def test_write_reports_tjpe_produces_expected_files(tmp_path: Path):
    paths = write_reports(
        result=_fake_tjpe_result(),
        gene=_fake_gene(),
        construct_str="GGGGS-FLAG",
        offtarget_reports=_fake_offtarget(),
        all_evaluated_variants=[_fake_tjpe_result().winning_variant],
        reports_dir=tmp_path,
        base_name="smoke",
    )
    assert paths["txt"].exists()
    txt = paths["txt"].read_text()
    assert "TJ-PE mode" in txt
    assert "COMPONENT BREAKDOWN" in txt
    assert "CAVEATS" in txt
    assert paths["csv"].exists()
    assert paths["offtarget"].exists()
    assert paths["all_variants"].suffix == ".gz"


def test_write_reports_pa_produces_expected_files(tmp_path: Path):
    paths = write_reports(
        result=_fake_pa_result(),
        gene=_fake_gene(),
        construct_str="GGGGS-FLAG-P2A-ZEOCINr",
        offtarget_reports=_fake_offtarget(),
        all_evaluated_variants=[],
        reports_dir=tmp_path,
        base_name="pa_smoke",
    )
    txt = paths["txt"].read_text()
    assert "Prime Assembly mode" in txt
    assert "LINEAR DONOR DNA" in txt
    assert "COMPONENT BREAKDOWN" in txt


def test_column_schemas_exposed():
    # Sanity: column constants are non-empty so downstream tooling can rely on them
    assert TJPE_CSV_COLUMNS[0] == "rank"
    assert PA_CSV_COLUMNS[0] == "rank"
    assert OFFTARGET_CSV_COLUMNS[0] == "guide_id"
    assert ALL_VARIANTS_CSV_COLUMNS[0] == "variant_rank"
