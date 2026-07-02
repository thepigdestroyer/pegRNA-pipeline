"""Manual cross-check verification report for a winning recoded TJ-PE design.

Produces a fully-annotated, colour-coded DOCX so a user can verify by eye that the
pipeline assembled a functional pegRNA with nothing missed. Four sections:

  1. The genomic context around the edit, annotated with every region the
     pegRNA / tjPE-pegRNA binds or uses and every region that will be changed.
  2. The full pegRNA sequence, broken into every component (Spacer, Scaffold,
     PBS2, RTT = insert + replacement gDNA, PBS1, PolyU).
  3. The predicted fully-edited genomic DNA, annotated with the inserted cargo
     and every recoded (changed) base.
  4. A table of every component sequence the pipeline used to assemble the insert
     and the winning pegRNA.

Design: a PURE, testable annotation model (:func:`build_verification_report` and
the ``AnnotatedSequence`` / ``Annotation`` dataclasses) is separated from the
python-docx renderer (:func:`render_verification_docx`). The model is where
verification correctness lives and is unit-tested without python-docx; the
renderer is a thin colour/highlight layer (python-docx is CX3/report-only).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD

__all__ = [
    "Annotation",
    "AnnotatedSequence",
    "ComponentRow",
    "VerificationReport",
    "build_verification_report",
    "render_verification_docx",
    "CATEGORY_STYLE",
]


@dataclass(frozen=True)
class Annotation:
    """A labelled region ``[start, end)`` (0-based, half-open) of a sequence.

    ``start == end`` denotes a zero-width point/insertion marker between bases.
    ``category`` keys into :data:`CATEGORY_STYLE` for colour + legend.
    """

    start: int
    end: int
    label: str
    category: str


@dataclass
class AnnotatedSequence:
    title: str
    sequence: str
    annotations: List[Annotation] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class ComponentRow:
    name: str
    sequence: str
    role: str
    part_of: str


@dataclass
class VerificationReport:
    title: str
    genomic_context: AnnotatedSequence
    pegrna: AnnotatedSequence
    edited_allele: AnnotatedSequence
    components: List[ComponentRow]
    metadata: Dict[str, object] = field(default_factory=dict)


# category -> (legend label, docx highlight WD_COLOR_INDEX name, render priority).
# Higher priority wins when annotations overlap on a base. Markers (zero-width)
# are listed for the legend but drawn as carets, not highlights.
CATEGORY_STYLE: Dict[str, Tuple[str, str, int]] = {
    "recode_edit": ("Recoded base changed (synonymous, protein-identical)", "PINK", 100),
    "insert": ("Inserted cargo / payload", "BRIGHT_GREEN", 90),
    "recode": ("Recoded replacement arm (synonymous region)", "YELLOW", 80),
    "locked": ("Non-coding, re-templated byte-exact (Decision 15)", "GRAY_25", 75),
    "spacer": ("Spacer (targets + nicks at S)", "TURQUOISE", 60),
    "pbs1": ("PBS1 (primes off the pegRNA nick 3'-OH at S)", "TEAL", 55),
    "pbs2": ("PBS2 (second-strand capture at the ngRNA nick)", "VIOLET", 55),
    "rtt": ("RTT (reverse-transcription template)", "YELLOW", 40),
    "scaffold": ("sgRNA scaffold (binds nCas9)", "GRAY_25", 30),
    "polyu": ("3' PolyU tail", "DARK_YELLOW", 30),
    "ngrna": ("ngRNA protospacer (second nick guide)", "RED", 50),
    "arm": ("Replacement / recode arm (re-templated span)", "GRAY_25", 20),
    "insert_site": ("Insert site S (cargo inserted here)", "RED", 0),
    "nick": ("ngRNA nick position", "RED", 0),
}


def _rc(seq: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "U": "A"}
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


def _rna_to_dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def _find_span(haystack: str, needle: str) -> Optional[Tuple[int, int]]:
    """Locate ``needle`` (or its reverse complement) in ``haystack``; span or None.

    Returns the first match. Tries the sense orientation first, then the reverse
    complement (the guide/PBS may be defined on the opposite strand).
    """
    if not needle:
        return None
    h = haystack.upper()
    n = _rna_to_dna(needle)
    i = h.find(n)
    if i >= 0:
        return (i, i + len(n))
    r = _rc(n)
    j = h.find(r)
    if j >= 0:
        return (j, j + len(r))
    return None


def _clamp_span(span: Optional[Tuple[int, int]], n: int) -> Optional[Tuple[int, int]]:
    if span is None:
        return None
    lo, hi = span
    lo, hi = max(0, min(lo, n)), max(0, min(hi, n))
    return (lo, hi) if hi >= lo else None


def _genomic_context_section(result, genome: str) -> AnnotatedSequence:
    genome = genome.upper()
    n = len(genome)
    anns: List[Annotation] = []
    notes: List[str] = []
    s = result.site_s
    nick = result.nick_position

    # Insert site S (cargo inserted here; the pegRNA nick is also at S).
    if 0 <= s <= n:
        anns.append(
            Annotation(s, s, "Insert site S (cargo inserted here; pegRNA nick at S)", "insert_site")
        )
    # Replacement / recode arm span between S and the ngRNA nick.
    lo, hi = sorted((s, nick))
    arm = _clamp_span((lo, hi), n)
    if arm and arm[1] > arm[0]:
        anns.append(Annotation(arm[0], arm[1], "Replacement / recode arm (re-templated)", "arm"))
    # Recoded (synonymous) and locked (non-coding byte-exact) spans.
    rec = _clamp_span(result.recode_span, n)
    if rec and rec[1] > rec[0]:
        anns.append(
            Annotation(rec[0], rec[1], "Recoded synonymous CDS (protein-identical)", "recode")
        )
    lock = _clamp_span(result.locked_span, n)
    if lock and lock[1] > lock[0]:
        anns.append(Annotation(lock[0], lock[1], "Non-coding, re-templated byte-exact", "locked"))
    # PBS1 primes off the L genomic bases ending at S (the pegRNA nick).
    pbs1_len = len(result.pbs1)
    if pbs1_len and s > 0:
        w = _clamp_span((s - pbs1_len, s), n)
        if w and w[1] > w[0]:
            anns.append(Annotation(w[0], w[1], "PBS1 binds here (pegRNA nick 3'-OH at S)", "pbs1"))
    # ngRNA nick + protospacer.
    if 0 <= nick <= n:
        anns.append(Annotation(nick, nick, "ngRNA nick position", "nick"))
    ng = _find_span(genome, result.ngrna_protospacer)
    if ng:
        anns.append(Annotation(ng[0], ng[1], "ngRNA protospacer (second nick guide)", "ngrna"))
    else:
        notes.append(
            "ngRNA protospacer not located in the local genomic window (may lie outside it)."
        )
    # pegRNA spacer / protospacer (nicks at S).
    sp = _find_span(genome, result.spacer)
    if sp:
        anns.append(Annotation(sp[0], sp[1], "pegRNA spacer / protospacer (nicks at S)", "spacer"))
    else:
        notes.append("pegRNA spacer not located as an exact substring of the local window.")
    return AnnotatedSequence("1. Genomic context around the edit", genome, anns, notes)


def _pegrna_section(result) -> AnnotatedSequence:
    seq = result.pegRNA_sequence
    anns: List[Annotation] = []
    notes: List[str] = []
    scaffold = STANDARD_SGRNA_SCAFFOLD
    # Authoritative 5'->3' order: Spacer . Scaffold . PBS2 . RTT . PBS1 . PolyU.
    i = 0
    ordered = [
        ("spacer", result.spacer, "Spacer"),
        ("scaffold", scaffold, "Scaffold"),
        ("pbs2", result.pbs2, "PBS2"),
    ]
    for cat, sub, label in ordered:
        if sub:
            anns.append(Annotation(i, i + len(sub), label, cat))
            i += len(sub)
    # RTT, split into its labelled sub-parts (insert / replacement gDNA).
    rtt_start = i
    if result.rtt_components:
        for label, sub in result.rtt_components:
            cat = "insert" if label == "insert" else "recode"
            pretty = "Insert (cargo)" if label == "insert" else "Replacement gDNA (recoded arm)"
            anns.append(Annotation(i, i + len(sub), pretty, cat))
            i += len(sub)
    else:
        anns.append(Annotation(i, i + len(result.rtt), "RTT", "rtt"))
        i += len(result.rtt)
    rtt_end = i
    for cat, sub, label in [("pbs1", result.pbs1, "PBS1"), ("polyu", result.polyu, "PolyU")]:
        if sub:
            anns.append(Annotation(i, i + len(sub), label, cat))
            i += len(sub)
    if i != len(seq):
        notes.append(
            f"component offsets sum to {i} nt but the assembled pegRNA is {len(seq)} nt "
            "-- verify the assembly order."
        )
    notes.append(
        f"tjPE template-jumping portion = PBS2.RTT.PBS1 (positions {rtt_start - len(result.pbs2)}..{rtt_end + len(result.pbs1)})."
    )
    return AnnotatedSequence("2. Full pegRNA, component breakdown", seq, anns, notes)


def _edited_allele_section(result, genome: str, cargo: str) -> AnnotatedSequence:
    genome = genome.upper()
    edited = (result.edited_allele or "").upper()
    anns: List[Annotation] = []
    notes: List[str] = []
    if not edited:
        notes.append("edited allele unavailable (result.edited_allele empty).")
        return AnnotatedSequence("3. Fully-edited genomic DNA (predicted)", edited, anns, notes)

    # The edited allele is built with a KNOWN geometry (deletion_size == 0, so the
    # recode is length-preserving), so the component spans are computed directly
    # rather than diff-inferred (a common-prefix/suffix diff would absorb the
    # matching tail of the arm and mis-size it):
    #   nick > S:  genome[:S]   + cargo + arm_sense + genome[nick:]
    #   nick < S:  genome[:nick]+ arm_sense + cargo + genome[S:]
    s, nick = result.site_s, result.nick_position
    cargo_dna = _rna_to_dna(cargo)
    arm_len = abs(nick - s)
    if nick > s:
        cargo_span = (s, s + len(cargo_dna))
        arm_e = (s + len(cargo_dna), s + len(cargo_dna) + arm_len)
    else:
        arm_e = (nick, nick + arm_len)
        cargo_span = (nick + arm_len, nick + arm_len + len(cargo_dna))

    cargo_span = _clamp_span(cargo_span, len(edited))
    if cargo_dna and cargo_span and cargo_span[1] > cargo_span[0]:
        anns.append(Annotation(cargo_span[0], cargo_span[1], "Inserted cargo / payload", "insert"))
        if edited[cargo_span[0] : cargo_span[1]] != cargo_dna:
            notes.append(
                "cargo block in the edited allele does not match the supplied cargo exactly."
            )

    arm_e = _clamp_span(arm_e, len(edited))
    orig_arm = genome[min(s, nick) : max(s, nick)]
    if arm_e and arm_e[1] > arm_e[0]:
        edited_arm = edited[arm_e[0] : arm_e[1]]
        anns.append(Annotation(arm_e[0], arm_e[1], "Recoded replacement arm", "recode"))
        if len(orig_arm) == len(edited_arm):
            for k, (a, b) in enumerate(zip(orig_arm, edited_arm)):
                if a != b:
                    pos = arm_e[0] + k
                    anns.append(Annotation(pos, pos + 1, f"recoded base {a}->{b}", "recode_edit"))
        else:
            notes.append(
                f"recoded-arm length ({len(edited_arm)}) != original arm length "
                f"({len(orig_arm)}); per-base change highlighting skipped."
            )
    return AnnotatedSequence("3. Fully-edited genomic DNA (predicted)", edited, anns, notes)


def _components(result, cargo: str) -> List[ComponentRow]:
    rows: List[ComponentRow] = [
        ComponentRow("spacer", result.spacer, "targets + nicks at S", "pegRNA"),
        ComponentRow("scaffold", STANDARD_SGRNA_SCAFFOLD, "sgRNA scaffold (binds nCas9)", "pegRNA"),
        ComponentRow(
            "PBS2", result.pbs2, "primer-binding site 2 (2nd-strand capture)", "pegRNA (tjPE jump)"
        ),
    ]
    if result.rtt_components:
        for label, seq in result.rtt_components:
            if label == "insert":
                rows.append(
                    ComponentRow("insert", seq, "insert / cargo payload (in RTT)", "pegRNA (RTT)")
                )
            else:
                rows.append(
                    ComponentRow(
                        "replacement_gDNA", seq, "recoded replacement arm (in RTT)", "pegRNA (RTT)"
                    )
                )
    else:
        rows.append(
            ComponentRow("RTT", result.rtt, "reverse-transcription template", "pegRNA (RTT)")
        )
    rows.append(
        ComponentRow(
            "PBS1", result.pbs1, "primer-binding site 1 (pegRNA nick)", "pegRNA (tjPE jump)"
        )
    )
    rows.append(ComponentRow("PolyU", result.polyu, "3' U tail", "pegRNA"))
    rows.append(
        ComponentRow(
            "full pegRNA", result.pegRNA_sequence, "assembled winning pegRNA oligo", "pegRNA"
        )
    )
    rows.append(
        ComponentRow(
            "ngRNA protospacer", result.ngrna_protospacer, "second-nick guide spacer", "ngRNA"
        )
    )
    if cargo:
        rows.append(
            ComponentRow("cargo (sense)", cargo, "cargo as supplied (sense orientation)", "input")
        )
    return rows


def build_verification_report(result, *, genome: str, cargo: str) -> VerificationReport:
    """Build the (pure, renderer-agnostic) verification model for a winning design."""
    metadata = {
        "site_s": result.site_s,
        "nick_position": result.nick_position,
        "deletion_size": result.deletion_size,
        "recode_edits": result.recode_edits,
        "pegRNA_length": len(result.pegRNA_sequence),
        "candidates_tried": result.candidates_tried,
    }
    return VerificationReport(
        title="Recoded TJ-PE pegRNA -- manual verification report",
        genomic_context=_genomic_context_section(result, genome),
        pegrna=_pegrna_section(result),
        edited_allele=_edited_allele_section(result, genome, cargo),
        components=_components(result, cargo),
        metadata=metadata,
    )


# --------------------------------------------------------------------- renderer


def _priority(cat: str) -> int:
    return CATEGORY_STYLE.get(cat, (cat, "GRAY_25", 10))[2]


def _base_category(section: AnnotatedSequence, pos: int) -> Optional[str]:
    """Highest-priority non-zero-width annotation covering base ``pos``."""
    best_cat, best_pri = None, -1
    for a in section.annotations:
        if a.end > a.start and a.start <= pos < a.end:
            pri = _priority(a.category)
            if pri > best_pri:
                best_cat, best_pri = a.category, pri
    return best_cat


def render_verification_docx(report: VerificationReport, out_path) -> str:
    """Render a :class:`VerificationReport` to a colour-coded DOCX (python-docx).

    Raises :class:`ImportError` if python-docx is not installed (it is CX3/report
    only; the model is testable without it). Returns the written path.
    """
    try:
        from docx import Document  # noqa: PLC0415
        from docx.enum.text import WD_COLOR_INDEX  # noqa: PLC0415
        from docx.shared import Pt, RGBColor  # noqa: PLC0415
    except Exception as exc:  # pragma: no cover - exercised only without python-docx
        raise ImportError(
            "render_verification_docx requires python-docx; install it in the report "
            f"environment. Original error: {exc}"
        ) from exc

    def hl(name: str):
        return getattr(WD_COLOR_INDEX, name, WD_COLOR_INDEX.GRAY_25)

    doc = Document()
    doc.add_heading(report.title, level=0)
    doc.add_paragraph(
        "Cross-check the pipeline output by eye: every region the pegRNA binds or "
        "changes is colour-coded below. Verify each component is present, correctly "
        "placed, and that only the intended bases change."
    )

    # Legend
    doc.add_heading("Legend", level=1)
    seen = set()
    for cat, (label, color, _pri) in sorted(CATEGORY_STYLE.items(), key=lambda kv: -kv[1][2]):
        if label in seen:
            continue
        seen.add(label)
        p = doc.add_paragraph()
        run = p.add_run("  SAMPLE  ")
        run.font.highlight_color = hl(color)
        run.font.name = "Consolas"
        p.add_run(f"  {label}")

    # Metadata
    doc.add_heading("Design summary", level=1)
    for k, v in report.metadata.items():
        doc.add_paragraph(f"{k}: {v}", style="List Bullet")

    def _emit_section(section: AnnotatedSequence):
        doc.add_heading(section.title, level=1)
        # annotation table
        real = [a for a in section.annotations if a.end > a.start]
        markers = [a for a in section.annotations if a.end == a.start]
        if real or markers:
            tbl = doc.add_table(rows=1, cols=4)
            tbl.style = "Light Grid Accent 1"
            hdr = tbl.rows[0].cells
            hdr[0].text, hdr[1].text, hdr[2].text, hdr[3].text = (
                "region",
                "len",
                "category",
                "label",
            )
            for a in sorted(section.annotations, key=lambda x: (x.start, x.end)):
                row = tbl.add_row().cells
                if a.end == a.start:
                    row[0].text = f"^{a.start}"
                    row[1].text = "0"
                else:
                    row[0].text = f"{a.start}-{a.end}"
                    row[1].text = str(a.end - a.start)
                row[2].text = a.category
                row[3].text = a.label
        for note in section.notes:
            doc.add_paragraph(f"Note: {note}", style="Intense Quote")
        # coloured sequence (monospace, wrapped), highlight = top-priority category
        seq = section.sequence
        marker_at = {a.start for a in markers}
        para = doc.add_paragraph()
        para.paragraph_format.space_after = Pt(6)
        width = 60
        for line_start in range(0, len(seq), width):
            chunk = seq[line_start : line_start + width]
            for off, base in enumerate(chunk):
                pos = line_start + off
                if pos in marker_at:
                    m = para.add_run("|")
                    m.font.name = "Consolas"
                    m.font.bold = True
                    m.font.color.rgb = RGBColor(0xC0, 0x00, 0x00)
                r = para.add_run(base)
                r.font.name = "Consolas"
                r.font.size = Pt(9)
                cat = _base_category(section, pos)
                if cat and CATEGORY_STYLE.get(cat, (None, None, 0))[2] >= 20:
                    r.font.highlight_color = hl(CATEGORY_STYLE[cat][1])
            para.add_run("\n")

    _emit_section(report.genomic_context)
    _emit_section(report.pegrna)
    _emit_section(report.edited_allele)

    # Section 4: component table
    doc.add_heading("4. Component sequences (insert + winning pegRNA assembly)", level=1)
    tbl = doc.add_table(rows=1, cols=5)
    tbl.style = "Light Grid Accent 1"
    hdr = tbl.rows[0].cells
    for c, t in zip(hdr, ("component", "len", "part of", "role", "sequence (5'->3')")):
        c.text = t
    for comp in report.components:
        row = tbl.add_row().cells
        row[0].text = comp.name
        row[1].text = str(len(comp.sequence))
        row[2].text = comp.part_of
        row[3].text = comp.role
        cell = row[4]
        cell.text = ""
        run = cell.paragraphs[0].add_run(comp.sequence)
        run.font.name = "Consolas"
        run.font.size = Pt(8)

    out = str(out_path)
    doc.save(out)
    return out
