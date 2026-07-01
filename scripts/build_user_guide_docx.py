#!/usr/bin/env python3
"""Generate docs/USER_GUIDE.docx - a non-technical, biologist-readable manual.

Run with a python that has python-docx installed:
    python scripts/build_user_guide_docx.py
Produces docs/USER_GUIDE.docx. Edit THIS file (not the .docx) and regenerate.
"""
from __future__ import annotations

from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt, RGBColor

OUT = Path(__file__).resolve().parents[1] / "docs" / "USER_GUIDE.docx"

doc = Document()

# ---- base styles ----
normal = doc.styles["Normal"]
normal.font.name = "Calibri"
normal.font.size = Pt(11)


def h1(t):
    doc.add_heading(t, level=1)


def h2(t):
    doc.add_heading(t, level=2)


def p(t="", **kw):
    para = doc.add_paragraph(t)
    if kw.get("italic"):
        para.runs[0].italic = True
    return para


def bullet(t):
    doc.add_paragraph(t, style="List Bullet")


def num(t):
    doc.add_paragraph(t, style="List Number")


def code(t):
    para = doc.add_paragraph()
    run = para.add_run(t)
    run.font.name = "Consolas"
    run.font.size = Pt(9.5)
    return para


# ---- title ----
title = doc.add_heading("pegRNA Design Pipeline", level=0)
sub = doc.add_paragraph("User Guide - designing template-jumping prime-editing pegRNAs")
sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
sub.runs[0].italic = True
sub.runs[0].font.size = Pt(12)
doc.add_paragraph()

# ---- 1. What this tool does ----
h1("1. What this tool does")
p(
    "This tool designs a prime-editing guide RNA (a pegRNA) that inserts a piece of "
    "DNA you choose (the cargo) at a precise position in a gene, in human iPSCs "
    "(the KOLF2.1J line). It uses template-jumping prime editing (tjPE), which lets "
    "you insert larger pieces than ordinary prime editing."
)
p(
    "Crucially, it does not just design something that will insert - it checks that "
    "the edit is safe. Before it gives you a design, it confirms the edit will not "
    "break the protein, will not disrupt splicing, will not collapse into a deletion, "
    "and will not misbehave elsewhere in the genome. If it cannot find a design that "
    "passes every check, it tells you why instead of handing you a bad one."
)
p(
    "You give it three things - a gene, a position, and a cargo - and it returns a "
    "ready-to-order pegRNA (plus its partner nicking guide) with a full breakdown of "
    "every part.",
)

# ---- 2. A few terms, in plain language ----
h1("2. A few terms, in plain language")
terms = [
    ("pegRNA", "the guide RNA that both targets the DNA and carries the sequence to be written in."),
    ("tjPE (template-jumping)", "a prime-editing variant that uses two nicks so it can write a larger insert."),
    ("cargo / insert", "the DNA you want to add - for example a tag, a selection marker, or a promoter."),
    ("replacement gDNA", "a short stretch next to the insert that the tool rewrites in a synonymous (protein-silent) way so the edit does not accidentally get deleted by the cell."),
    ("PBS1 and PBS2", "two short 'primer binding' regions the pegRNA uses to grip the DNA during the edit."),
    ("spacer and scaffold", "the standard parts every guide RNA has: the spacer sets the target, the scaffold binds the editing enzyme."),
    ("nick", "a controlled single-strand cut the editor makes; tjPE uses two of them."),
]
for name, desc in terms:
    para = doc.add_paragraph(style="List Bullet")
    r = para.add_run(name + " - ")
    r.bold = True
    para.add_run(desc)

# ---- 3. What each safety check protects ----
h1("3. What the safety checks protect against")
p("Each design is put through a series of checks (\"gates\"). In plain terms:")
gates = [
    ("Deletion protection", "stops the insert from being lost to the cell's repair machinery, by rewriting a short neighbouring region so it no longer matches."),
    ("Protein protection", "makes sure any rewriting is silent - the protein sequence is unchanged."),
    ("Splicing protection", "checks the edit does not create or destroy a splice site (using SpliceAI)."),
    ("Structure protection", "checks the edit does not land on a structurally critical part of the protein (using AlphaFold)."),
    ("Guide-folding check", "makes sure the pegRNA folds properly and does not tangle up its own working parts."),
    ("Enzyme-binding check", "makes sure the editing enzyme can still dock onto the guide."),
    ("Off-target check", "makes sure the nicking guide is specific and does not cut elsewhere in the genome (using Cas-OFFinder)."),
]
for name, desc in gates:
    para = doc.add_paragraph(style="List Bullet")
    r = para.add_run(name + " - ")
    r.bold = True
    para.add_run(desc)
p(
    "If a design fails a check, the tool moves on to the next candidate. If nothing "
    "passes, it stops and reports which check was the problem - it never returns an "
    "unsafe design.",
    italic=True,
)

# ---- 4. What you provide ----
h1("4. What you provide (the design request)")
p(
    "You describe your edit in a small text file (a \"config\"). A complete, ready-to-run "
    "example ships with the tool at examples/ctcf_cterm/config.json (a CTCF C-terminal "
    "tag knock-in). You can copy it and change the values. The main things you set:"
)
inputs = [
    ("Gene / transcript", "which gene and transcript you are editing (for the records and the protein-structure check)."),
    ("Target site", "the exact position where the cargo goes in."),
    ("Cargo", "the DNA to insert - e.g. a GS linker, a FLAG tag, a T7 element, a P2A peptide, and a Zeocin resistance marker."),
    ("Strand and side", "which strand you are targeting and which side of the insert the rewritten region sits on."),
    ("The local DNA sequence", "the stretch of genome around your target, and which bases are coding."),
]
for name, desc in inputs:
    para = doc.add_paragraph(style="List Bullet")
    r = para.add_run(name + " - ")
    r.bold = True
    para.add_run(desc)
p(
    "You do not need to design any of the pegRNA parts yourself - the tool works those "
    "out. You only describe the biological goal."
)

# ---- 5. How to run it ----
h1("5. How to run it")
h2("Option A - a quick dry run (any computer, no special hardware)")
p(
    "This runs the full design logic and produces a real report, using stand-in "
    "versions of the heavy checks. It is the fastest way to see a design and confirm "
    "your request is set up correctly:"
)
code("pegrna design-recoded --config examples/ctcf_cterm/config.json --mode fake --output-dir results/dry")
h2("Option B - a full run with all checks (on the HPC)")
p(
    "The full biological checks use large models that need the research computer "
    "(the HPC). Setup is one command (see the installation guide), after which you "
    "submit a run:"
)
code("bash scripts/hpc/install.sh")
code("qsub -v PEGRNA_DATA=$PEGRNA_DATA,CONFIG=examples/ctcf_cterm/config.json \\\n     scripts/hpc/run_recoded_pipeline.pbs")
p(
    "Full step-by-step install and run instructions are in scripts/hpc/README.md, and "
    "common problems and fixes are in docs/TROUBLESHOOTING.md.",
    italic=True,
)

# ---- 6. How to read the output ----
h1("6. How to read the output")
p("Each run writes two files to your output folder:")
bullet("recoded_tjpe_report.txt - a readable summary of the design.")
bullet(
    "recoded_tjpe_components.csv - the same design as a table (one row per part), "
    "easy to open in Excel."
)
p("The design is reported 5'->3' in this order:")
code("Spacer . Scaffold . PBS2 . RTT . PBS1 . PolyU")
p("The RTT (the written-in part) is broken down into two labelled pieces:")
bullet("insert - your cargo.")
bullet(
    "replacement gDNA - the short neighbouring region, rewritten silently to prevent "
    "deletion."
)
p(
    "The report also gives the target position, the nicking guide, how many silent "
    "changes were made, and any advisory warnings. The table's \"part_of\" column tells "
    "you which pieces are the template-jumping core (PBS2, RTT, PBS1) and which complete "
    "the full guide (spacer, scaffold)."
)

# ---- 7. Before you order DNA ----
h1("7. Before you order DNA")
p(
    "The tool's safety thresholds are set from the published literature but have not "
    "yet been tuned against wet-lab results, and the heavy checks should be confirmed "
    "on the HPC. Treat the output as a strong, checked starting point rather than a "
    "final validated design, and confirm with your own judgement before ordering.",
)
p("Two principles the tool follows that are worth knowing:")
bullet(
    "Any DNA it rewrites near your insert is changed silently (the protein never "
    "changes)."
)
bullet(
    "Validation comes first - a prediction is only trusted once it has passed the "
    "checks."
)

# ---- 8. Getting help ----
h1("8. Getting help")
bullet("Installation and running on the HPC: scripts/hpc/README.md")
bullet("Common problems and fixes: docs/TROUBLESHOOTING.md")
bullet("What the tool can do and how it works (more technical): docs/FEATURES.md")
bullet("Why it is built the way it is (design history): docs/DESIGN_HISTORY.md")

doc.save(str(OUT))
print("wrote", OUT, f"({OUT.stat().st_size} bytes)")
