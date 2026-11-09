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
sub = doc.add_paragraph("User Guide - designing prime-editing pegRNAs (tjPE and Prime Assembly)")
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
p(
    "There are two ways to carry the insert. A tag or a short insert goes inside the "
    "pegRNA (tjPE, which most of this guide describes). A long insert, such as a whole "
    "cassette or gene, would make the pegRNA too long to check, so the tool also offers "
    "Prime Assembly: two short pegRNAs plus a separate piece of donor DNA that carries "
    "the insert. Section 8 explains it."
)

# ---- 2. A few terms, in plain language ----
h1("2. A few terms, in plain language")
terms = [
    ("pegRNA", "the guide RNA that both targets the DNA and carries the sequence to be written in."),
    ("tjPE (template-jumping)", "a prime-editing variant that uses two nicks so it can write a larger insert."),
    ("cargo / insert", "the DNA you want to add - for example a tag, a selection marker, or a promoter."),
    ("replacement gDNA", "a short stretch next to the insert that the tool rewrites in a synonymous (protein-silent) way so the edit does not accidentally get deleted by the cell."),
    ("PBS1 and rcPBS2", "two short regions of the tjPE pegRNA. PBS1 grips the DNA at the pegRNA's nick. rcPBS2 is the reverse complement of PBS2: the editor copies it onto the new DNA, where it captures the strand opened by the second nick (the 'template jump')."),
    ("spacer and scaffold", "the standard parts every guide RNA has: the spacer sets the target, the scaffold binds the editing enzyme."),
    ("nick", "a controlled single-strand cut the editor makes; tjPE uses two of them."),
    ("Prime Assembly", "a way to insert a long piece of DNA: two short pegRNAs make two nicks, and a separate donor DNA carries the insert."),
    ("donor DNA", "the piece of DNA you order or make that carries the insert in Prime Assembly. Its two ends match the two flaps."),
    ("flap", "the short stretch (about 35 nucleotides) that each Prime Assembly pegRNA writes at its nick, and that one end of the donor sticks to."),
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
p(
    "You do not need an AI assistant or any special agent to run this tool - it is an "
    "ordinary command you type (or that a colleague can run for you). There are three "
    "ways to run it, from simplest to most automated."
)
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
h2("Option C - let the tool drive the HPC for you (optional)")
p(
    "If you would rather not log in to the research computer yourself, the tool can do "
    "it for you: it checks your request, submits the job, waits for it to finish, and "
    "downloads the results to your machine. You set this up once, then start runs with "
    "a single command:"
)
code("pegrna setup                      # once - saves how to reach the HPC")
code('pegrna run --gene CTCF --construct "<your construct>" --mode tjpe --output-dir results')
p("And to keep track of runs:")
code("pegrna status      # see your runs\npegrna tail <run-id>   # watch progress\npegrna cancel <run-id> # stop a run")
p(
    "This is just convenience automation (it logs in and submits for you) - it is not "
    "an AI making decisions. It currently drives the tjPE / PE7 workflows; the "
    "recoded design is run with Option A or B above, and Prime Assembly has its own "
    "command (section 8).",
    italic=True,
)

# ---- 6. How to read the output ----
h1("6. How to read the output")
p("Each tjPE run writes these files to your output folder:")
bullet("recoded_tjpe_report.txt - a readable summary of the design.")
bullet(
    "recoded_tjpe_components.csv - the same design as a table (one row per part), "
    "easy to open in Excel."
)
bullet(
    "recoded_tjpe_changes.csv - every single nucleotide that differs from what you gave "
    "the tool, with its part and codon, so you can check each change."
)
bullet("recoded_tjpe_gates.json - the result of every check, for scripts and records.")
bullet(
    "recoded_tjpe_verification.docx - a colour-annotated manual cross-check (written "
    "only when the optional python-docx package is installed)."
)
p("The design is reported 5'->3' in this order:")
code("Spacer . Scaffold . rcPBS2 . RTT . PBS1 . PolyU")
p("The RTT (the written-in part) is broken down into two labelled pieces:")
bullet("insert - your cargo.")
bullet(
    "replacement gDNA - the short neighbouring region, rewritten silently to prevent "
    "deletion."
)
p(
    "The report also gives the target position, the nicking guide, how many silent "
    "changes were made, and any advisory warnings. The table's \"part_of\" column tells "
    "you which pieces are the template-jumping core (rcPBS2, RTT, PBS1) and which complete "
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

# ---- 8. Prime Assembly ----
h1("8. Prime Assembly - for a cassette or a gene")
p(
    "A tag or a short insert is written by the pegRNA itself (tjPE, above). A long "
    "insert is different: the pegRNA would need several hundred extra nucleotides, "
    "longer than the enzyme-binding check can judge. Prime Assembly gets around this by "
    "moving the insert out of the pegRNA."
)
p(
    "Two short pegRNAs (about 147 nucleotides each, whatever the size of the insert) "
    "make two nicks on opposite strands, one on each side of your insertion point. Each "
    "writes a short flap. A separate piece of donor DNA carries your insert, and its two "
    "ends match the two flaps, so the cell stitches the donor in. In this tool the donor "
    "also puts back the DNA that lies between the two nicks (rewritten silently where it "
    "codes for protein), so nothing is deleted - the same promise as tjPE."
)
h2("When to use it")
bullet(
    "Your insert is at least two flaps long: 70 nucleotides with the default "
    "35-nucleotide flap. A shorter insert is refused, with a pointer back to tjPE."
)
bullet("For a tag or a short insert, tjPE is the better fit.")
h2("What it does for you")
p(
    "You do not choose the guides. The tool looks at every pair of guides that face each "
    "other across your insertion point, drops the ones that break the guide rules or "
    "match elsewhere in the genome, ranks the rest by how few known risks each carries, "
    "and tries them best first. Both pegRNAs of a pair go through the same kinds of "
    "checks as in section 3 (folding, enzyme binding, splicing), the donor is checked "
    "for silent recoding and for deletion, and the first pair that passes everything is "
    "returned. If none does, the tool lists, step by step, what stopped each pair."
)
h2("How to run it")
p("A dry run on any computer, to see a design and confirm your request is set up correctly:")
code("pegrna design-pa --config examples/ctcf_cterm/config_pa.json --mode fake --output-dir results/dry_pa")
p(
    "A full run with the real models (--mode api) needs a machine that has them "
    "installed. The web page runs on one, so it is the easiest way: in the form, step 3 "
    "asks how the insert should be written - choose \"As a separate donor DNA\"."
)
h2("What you get")
p("Each Prime Assembly run writes five files to your output folder:")
bullet(
    "recoded_pa_report.txt - a readable summary: both pegRNAs, the donor, where the two "
    "nicks are, and the result of every check."
)
bullet("recoded_pa_components.csv - the same design as a table, one row per part.")
bullet(
    "recoded_pa_changes.csv - every single nucleotide that differs from what you gave "
    "the tool. Only the donor's genomic arms can change; your insert is used exactly as "
    "given."
)
bullet(
    "recoded_pa_donor.fasta - the three sequences to order: pegRNA A, pegRNA B and the "
    "donor. A pegRNA is RNA, so read each T as U."
)
bullet("recoded_pa_gates.json - the result of every check, for scripts and records.")
h2("Things to know before you rely on it")
bullet(
    "It is untested at the bench. The published method replaced the DNA between the "
    "nicks with other DNA. A donor that puts the original DNA back, as here, is this "
    "project's own extension and nobody has tested it. Every result says so."
)
bullet(
    "A guide's recognition sequence (its PAM) can come back with the restored DNA, so the "
    "guide could cut the edited allele again. The tool breaks a PAM with a silent change "
    "where it can, and reports the rest: \"weakened\" means only a weak PAM is left (the "
    "enzyme still cuts it, at a fraction of its normal rate); \"regenerated\" means the "
    "PAM is back as it was."
)
bullet(
    "The strict off-target setting (no near-match allowed anywhere) can leave very few "
    "guide pairs. In the CTCF example only 1 of 24 survives. The \"Off-target hits "
    "allowed\" setting admits more, and every hit it allows is listed in the result."
)
bullet(
    "A full run takes a few minutes. The result page and the report also carry notes on "
    "making the donor, taken from the published method (Liu et al., Nature 2026)."
)

# ---- 9. Getting help ----
h1("9. Getting help")
bullet("Installation and running on the HPC: scripts/hpc/README.md")
bullet("Common problems and fixes: docs/TROUBLESHOOTING.md")
bullet("What the tool can do and how it works (more technical): docs/FEATURES.md")
bullet("Why it is built the way it is (design history): docs/DESIGN_HISTORY.md")
bullet(
    "Prime Assembly - how it works and what the real runs showed: "
    "docs/decisions/0030-prime-assembly-design-pa.md"
)

doc.save(str(OUT))
print("wrote", OUT, f"({OUT.stat().st_size} bytes)")
