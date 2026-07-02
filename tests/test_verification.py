"""Verification report model (manual cross-check DOCX). The pure annotation model
is where verification correctness lives, so it is tested exhaustively here; the
python-docx renderer is tested only if python-docx is installed (CX3/report-only).
"""

import pytest

from pegrna_pipeline.hpc.recoded_tjpe_workflow import RecodedTJPEResult
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD
from pegrna_pipeline.hpc.verification import build_verification_report

PREFIX = "TTTTTTTTTT"  # 10
ARM_ORIG = "GAAGAAGAAG"  # 10, genome[10:20]
SUFFIX = "CCCCCCCCCC"  # 10
GENOME = PREFIX + ARM_ORIG + SUFFIX  # 30 nt; S=10, nick=20
ARM_SENSE = "GAGGAAGAAG"  # recoded arm: single change vs ARM_ORIG at index 2 (A->G)
CARGO = "ATGCATGC"  # 8
EDITED = PREFIX + CARGO + ARM_SENSE + SUFFIX  # 38 nt (nick>S: prefix+cargo+arm+suffix)


def _result():
    spacer, pbs2, pbs1, polyu = "UUUUUUUUUU", "TTTT", "AAAA", "UUU"
    rtt_components = [("insert", CARGO), ("replacement_gdna", ARM_SENSE)]
    rtt = CARGO + ARM_SENSE
    peg = spacer + STANDARD_SGRNA_SCAFFOLD + pbs2 + rtt + pbs1 + polyu
    return RecodedTJPEResult(
        pegRNA_sequence=peg,
        spacer=spacer,
        pbs2=pbs2,
        rtt=rtt,
        pbs1=pbs1,
        polyu=polyu,
        ngrna_protospacer="CCCCCCCCCC",  # locates in the suffix (20..30)
        site_s=10,
        nick_position=20,
        deletion_size=0,
        recode_edits=1,
        candidates_tried=1,
        rtt_components=rtt_components,
        insert_seq=CARGO,
        replacement_gdna=ARM_SENSE,
        edited_allele=EDITED,
        recode_span=(10, 20),
        locked_span=None,
    )


def _spans(section):
    return {(a.category, a.start, a.end) for a in section.annotations}


def test_genomic_context_annotations():
    rep = build_verification_report(_result(), genome=GENOME, cargo=CARGO)
    spans = _spans(rep.genomic_context)
    assert ("insert_site", 10, 10) in spans
    assert ("arm", 10, 20) in spans
    assert ("recode", 10, 20) in spans
    assert ("pbs1", 6, 10) in spans  # genome[S-len(pbs1):S]
    assert ("nick", 20, 20) in spans
    assert ("ngrna", 20, 30) in spans
    assert ("spacer", 0, 10) in spans  # UUUUUUUUUU -> TTTTTTTTTT at prefix


def test_pegrna_component_offsets_tile_the_construct():
    rep = build_verification_report(_result(), genome=GENOME, cargo=CARGO)
    seq = rep.pegrna.sequence
    cats = {a.category: (a.start, a.end) for a in rep.pegrna.annotations}
    assert cats["spacer"] == (0, 10)
    assert cats["recode"][0] == cats["insert"][1]  # replacement gDNA follows the insert in the RTT
    assert max(a.end for a in rep.pegrna.annotations) == len(seq)
    assert not any("offsets sum" in n for n in rep.pegrna.notes)


def test_edited_allele_cargo_and_recode_edits():
    rep = build_verification_report(_result(), genome=GENOME, cargo=CARGO)
    spans = _spans(rep.edited_allele)
    assert ("insert", 10, 18) in spans  # cargo block
    assert ("recode", 18, 28) in spans  # recoded arm block
    edits = [a for a in rep.edited_allele.annotations if a.category == "recode_edit"]
    assert len(edits) == 1
    assert (edits[0].start, edits[0].end) == (20, 21)  # the single changed base


def test_components_complete():
    rep = build_verification_report(_result(), genome=GENOME, cargo=CARGO)
    names = [c.name for c in rep.components]
    for expected in [
        "spacer",
        "scaffold",
        "PBS2",
        "insert",
        "replacement_gDNA",
        "PBS1",
        "PolyU",
        "full pegRNA",
        "ngRNA protospacer",
        "cargo (sense)",
    ]:
        assert expected in names
    fp = next(c for c in rep.components if c.name == "full pegRNA")
    assert len(fp.sequence) == len(rep.pegrna.sequence)


def test_nick_5prime_of_s_geometry():
    # nick < S: edited = genome[:nick] + arm_sense + cargo + genome[S:]
    r = _result()
    r.site_s = 20
    r.nick_position = 10
    r.recode_span = (10, 20)
    r.edited_allele = PREFIX + ARM_SENSE + CARGO + SUFFIX
    rep = build_verification_report(r, genome=GENOME, cargo=CARGO)
    spans = _spans(rep.edited_allele)
    assert ("recode", 10, 20) in spans  # arm at [nick, nick+arm_len)
    assert ("insert", 20, 28) in spans  # cargo after the arm


def test_docx_renders_when_available(tmp_path):
    docx = pytest.importorskip("docx")  # python-docx is CX3/report-only
    from pegrna_pipeline.hpc.verification import render_verification_docx

    rep = build_verification_report(_result(), genome=GENOME, cargo=CARGO)
    out = render_verification_docx(rep, tmp_path / "verify.docx")
    assert (tmp_path / "verify.docx").exists()
    text = "\n".join(p.text for p in docx.Document(out).paragraphs).lower()
    assert "verification" in text
