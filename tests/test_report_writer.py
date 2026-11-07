"""Tests for the Plan 4 annotator + report writers."""

from __future__ import annotations

from pathlib import Path

import pytest

from pegrna_pipeline.hpc.annotator import (
    ComponentSource,
    annotate_pe7_result,
    annotate_tjpe_result,
)
from pegrna_pipeline.hpc.data_types import GeneRecord, PegRNADesign, VariantScore
from pegrna_pipeline.hpc.offtarget_scorer import OffTargetReport
from pegrna_pipeline.hpc.pe7_workflow import PE7DesignPair, PE7Result
from pegrna_pipeline.hpc.report_writer import (
    ALL_VARIANTS_CSV_COLUMNS,
    OFFTARGET_CSV_COLUMNS,
    PE7_CSV_COLUMNS,
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
        pegrna_pe2_score=42.3,
        pegrna_pe2_cell_line="k562",
        edit_position_in_mrna=200,
        ngrna_nick_in_mrna=150,
        all_evaluated_variants=[variant],
    )


def _fake_pe7_result() -> PE7Result:
    variant = VariantScore(
        sequence_name="pe7_var_0",
        codon_variant_dna="GACTACAAGGACGACGATGAC",
    )
    variant.pegRNA_designs = [_fake_pegrna(name="pe7_top", pe2=55.0, strand="Fw")]
    return PE7Result(
        winning_variant=variant,
        full_insert_dna="GACTACAAGGACGACGATGAC",  # 21 bp FLAG-tag DNA
        pegRNA_full="PE7_PEG_FULL",
        spacer="GCAGCAGCATGCATGCAGCT",
        pbs="TCACCAGCAGCAGC",
        rt="GACAAGGACGACGAT",
        pegrna_strand="Fw",
        pegrna_pe2_score=55.0,
        pegrna_pe2_cell_line="k562",
        ngrna_spacer="GCTAGCATGCAGCATGCATC",
        ngrna_distance_to_edit=-50,
        ngrna_strand="Rv",
        ngrna_deepcas9_score=0.91,
        edit_position_in_mrna=200,
        ngrna_nick_in_mrna=150,
        all_evaluated_variants=[variant],
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


def test_annotate_pe7_yields_expected_components():
    annotations = annotate_pe7_result(_fake_pe7_result(), _fake_gene())
    names = [a.name for a in annotations]
    assert names == [
        "spacer",
        "sgRNA_scaffold",
        "RT_template",
        "PBS",
        "ngRNA_spacer",
        "ngRNA_scaffold",
    ]
    rt = next(a for a in annotations if a.name == "RT_template")
    assert rt.source == ComponentSource.PRIDICT2


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


def test_write_reports_refuses_a_result_it_has_no_writer_for(tmp_path: Path):
    """Prime Assembly writes its own files (``pa_report``). Its result handed to this writer would be a mistake, and
    it must fail loudly instead of writing a report of the wrong kind."""
    with pytest.raises(TypeError, match="unknown result type"):
        write_reports(
            result=object(),
            gene=_fake_gene(),
            construct_str="FLAG",
            offtarget_reports=[],
            all_evaluated_variants=[],
            reports_dir=tmp_path,
            base_name="nope",
        )
    assert not list(tmp_path.glob("nope*")), "nothing is written for a result it does not know"


@pytest.mark.parametrize("make_result", [_fake_tjpe_result, _fake_pe7_result])
def test_the_shared_caveats_point_at_design_pa_and_not_at_a_removed_mode(
    tmp_path: Path, make_result
):
    """The caveats are printed under both the TJ-PE and the PE7 report. They used to send a long insert to a
    "Prime Assembly mode" of ``run`` and to quote that prototype's own settings; that mode is gone.
    """
    paths = write_reports(
        result=make_result(),
        gene=_fake_gene(),
        construct_str="FLAG",
        offtarget_reports=[],
        all_evaluated_variants=[],
        reports_dir=tmp_path,
        base_name="caveats",
    )
    txt = paths["txt"].read_text()
    caveats = txt.split("CAVEATS")[-1]
    assert "pegrna design-pa" in caveats
    for gone in ("Prime Assembly mode", "PA 30 bp", "house rule", "30-80 bp"):
        assert gone not in txt, gone


def test_write_reports_pe7_produces_expected_files(tmp_path: Path):
    paths = write_reports(
        result=_fake_pe7_result(),
        gene=_fake_gene(),
        construct_str="FLAG",
        offtarget_reports=_fake_offtarget(),
        all_evaluated_variants=[_fake_pe7_result().winning_variant],
        reports_dir=tmp_path,
        base_name="pe7_smoke",
    )
    assert paths["txt"].exists()
    txt = paths["txt"].read_text()
    assert "PE7" in txt
    assert "COMPONENT BREAKDOWN" in txt
    assert "CAVEATS" in txt
    assert paths["csv"].exists()
    csv_text = paths["csv"].read_text()
    assert "pegrna_spacer" in csv_text
    assert "GCAGCAGCATGCATGCAGCT" in csv_text
    assert paths["offtarget"].exists()
    assert paths["all_variants"].suffix == ".gz"


def test_write_reports_pe7_with_top_pairs_shows_orderable_sequences(tmp_path: Path):
    """top_pairs table + ORDERABLE SEQUENCES block are rendered correctly."""
    pair = PE7DesignPair(
        rank=1,
        variant_sequence_name="pe7_var_0",
        pegRNA_full="PE7_PEG_FULL_OLIGO",
        spacer="GCAGCAGCATGCATGCAGCT",
        pbs="TCACCAGCAGCAGC",
        rt="GACAAGGACGACGAT",
        pegrna_strand="Fw",
        pegrna_pe2_score=55.0,
        pegrna_pe2_cell_line="k562",
        ngrna_spacer="GCTAGCATGCAGCATGCATC",
        ngrna_distance_to_edit=-50,
        ngrna_strand="Rv",
        ngrna_deepcas9_score=0.91,
        ngrna_nick_in_mrna=150,
        joint_pe3_score=0.50,
    )
    result = _fake_pe7_result()
    result.top_pairs.append(pair)
    paths = write_reports(
        result=result,
        gene=_fake_gene(),
        construct_str="FLAG",
        offtarget_reports=[],
        all_evaluated_variants=[],
        reports_dir=tmp_path,
        base_name="pe7_top_pairs",
    )
    txt = paths["txt"].read_text()
    assert "ORDERABLE SEQUENCES" in txt
    assert "PE7_PEG_FULL_OLIGO" in txt
    # ngRNA oligo must be spacer + scaffold (not just the scaffold alone)
    from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD

    expected_ngrna = "GCTAGCATGCAGCATGCATC" + STANDARD_SGRNA_SCAFFOLD
    assert expected_ngrna in txt


def test_write_reports_pe7_empty_result_no_misleading_ngrna(tmp_path: Path):
    """An empty PE7 result must not print the scaffold alone as a valid ngRNA oligo."""
    from pegrna_pipeline.hpc.pe7_workflow import _empty_result

    gene = _fake_gene()
    result = _empty_result(
        gene=gene,
        insert_full_dna="ACGT",
        variants=[],
        score_cell_line="k562",
    )
    paths = write_reports(
        result=result,
        gene=gene,
        construct_str="ACGT",
        offtarget_reports=[],
        all_evaluated_variants=[],
        reports_dir=tmp_path,
        base_name="pe7_empty",
    )
    txt = paths["txt"].read_text()
    # The ngRNA oligo line must show "<not assembled>", not the 76-nt scaffold alone.
    assert "FULL ngRNA OLIGO (0 nt):" in txt
    assert "<not assembled>" in txt
    # The label must never be followed by just the bare scaffold as the "oligo"
    assert "FULL ngRNA OLIGO (76 nt):" not in txt


def test_column_schemas_exposed():
    # Sanity: column constants are non-empty so downstream tooling can rely on them
    assert TJPE_CSV_COLUMNS[0] == "rank"
    assert PE7_CSV_COLUMNS[0] == "rank"
    assert OFFTARGET_CSV_COLUMNS[0] == "guide_id"
    assert ALL_VARIANTS_CSV_COLUMNS[0] == "variant_rank"
