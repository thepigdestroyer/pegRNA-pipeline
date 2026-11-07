"""Reports for a Prime Assembly design: the text report, the components CSV, the FASTA, the JSON the web UI reads.

A reader orders DNA from these, so three things are never left to chance: a stand-in gate is never passed off as
a verdict (the banner comes before the order sheet), the standing caveat that a donor restoring the span is this
project's extension of the method is near the top and not buried in the warnings, and every sequence to order is
printed in full.
"""

from __future__ import annotations

import csv
import io
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from pegrna_pipeline.hpc.recoded_pa_workflow import STANDING_CAVEAT, PAGuideDesign, RecodedPAResult
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD

__all__ = ["render_pa_report", "render_pa_csv", "render_pa_fasta", "pa_gates_payload"]

_STRAND_WORD = {"top": "top strand", "bottom": "bottom strand"}

_PAM_WORDS = {
    "broken": "broken in the donor by one synonymous change, so the guide cannot cut the product",
    "broken-by-recode": "broken by the recoding of the arm, so the guide cannot cut the product",
    "weakened": "WEAKENED: it reads {pam} in the product, which SpCas9 still cuts at a fraction of its NGG rate (about a fifth for NAG, less for NGA), so the guide may still nick the product",
    "intact-noncoding": "REGENERATED: it lies in non-coding DNA, which is kept exactly, so the guide can cut the product again",
    "intact-no-synonymous-change": "REGENERATED: it lies in coding sequence with no synonymous change that breaks it, so the guide can cut the product again",
}

_DONOR_PREP = (
    "The paper (Liu et al., Nature 2026) used a double-stranded donor made by PCR, whose two ends match the two "
    "flaps. A 0.8 kb donor gave about 50% insertion at AAVS1 in HEK293T cells (ddPCR); 2.2 kb and 4 kb donors also "
    "worked there, and a 0.8 kb donor worked at TRAC in Jurkat cells. The flaps can only anneal once the donor's "
    "5' ends are displaced or resected.",
    "A donor with 3' overhangs (PCR, then lambda exonuclease, with five internal phosphorothioate linkages next to "
    "the region that matches the flap so the exonuclease does not resect too far) gave 20-40% insertion: slightly "
    "less than the plain PCR product, and slightly more precise.",
    "A DNA-PK inhibitor, AZD-7648, about doubled insertion and lowered the indel rate. The editor was PE6c with two "
    "epegRNAs, each with a 35 nt flap.",
    "None of those experiments used a donor that restores the DNA between the nicks, as this design does, so what "
    "to expect at your locus has to be measured.",
)


def _pam_sentence(g: PAGuideDesign) -> str:
    """What happens to this guide's PAM in the product, in words (naming the PAM it really has)."""
    return _PAM_WORDS.get(g.pam_status, g.pam_status).replace("{pam}", g.product_pam)


def _fasta_wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width)) or ""


def _guide_block(side: str, g: PAGuideDesign) -> List[str]:
    where = "5'" if side == "A" else "3'"
    end = f"the donor's {where} end"
    lines = [
        f"  pegRNA {side}  ({_STRAND_WORD[g.strand]}, nick between window nt {g.cut} and {g.cut + 1}; "
        f"Spacer.Scaffold.RTT.PBS.PolyU, {len(g.pegRNA_sequence)} nt)",
        f"    spacer  : {g.spacer}",
        f"    RTT     : {g.rtt}",
        f"              (writes flap {side}, {g.flap_len} nt = {end})",
        f"    PBS     : {g.pbs}",
        f"              ({g.pbs_len} nt, anneals to the 3' end of the nicked {g.strand} strand)",
        f"    PolyU   : {g.tail}",
        f"    full    : {g.pegRNA_sequence}",
    ]
    if g.nick_pairing_bp:
        k = g.nick_pairing_bp
        lines.insert(
            6,
            f"              (the RTT's last {k} nt pair with the spacer's last {k} nt, because the donor copies the "
            "genome at this nick: expected, and not counted as a misfold)",
        )
    return lines


def render_pa_report(result: RecodedPAResult, manifest: Mapping[str, Any]) -> str:
    """The human-readable report. ``manifest`` carries provenance (and ``stand_in_gates`` / ``scorer_mode``)."""
    ga, gb = result.guide_a, result.guide_b
    n1, n2 = result.span
    out: List[str] = ["=== Prime Assembly design: two pegRNAs and one donor ===", ""]

    stand_in = str(manifest.get("stand_in_gates", "") or "")
    if stand_in:
        out += [
            "*** NOT A VALIDATED DESIGN - STAND-IN GATES IN USE ***",
            f"  The following gates did NOT run a real model: {stand_in}",
            "  Their verdicts below are placeholders, not evidence.",
        ]
        if str(manifest.get("scorer_mode", "")) == "fake":
            out.append("  scorer_mode=fake: this whole run is a dry run of the design logic.")
        else:
            out.append("  Re-run with the missing credentials/resources before ordering DNA.")
        out.append("")

    out += [f"NOTE: {STANDING_CAVEAT}", ""]

    out += ["WHAT TO ORDER (all sequences 5'->3'; a pegRNA is RNA, so read each T as U)"]
    out += _guide_block("A", ga)
    out.append("")
    out += _guide_block("B", gb)
    out.append("")
    out.append(
        f"  donor DNA  ({len(result.donor)} nt, double-stranded; the sense strand is shown; "
        "arm5 + cargo + arm3)"
    )
    labels = {"arm5": "arm5  ", "cargo": "cargo ", "arm3": "arm3  "}
    notes = {
        "arm5": f"{len(result.arm5)} nt: the restored copy of window nt {n1 + 1}-{result.site_s}",
        "cargo": f"{len(result.cargo)} nt: the insert, used exactly as you gave it",
        "arm3": f"{len(result.arm3)} nt: the restored copy of window nt {result.site_s + 1}-{n2}",
    }
    for name, seq in result.donor_parts:
        out.append(f"    {labels[name]} : {seq}")
        out.append(f"               ({notes[name]})")
    out.append(f"    full   : {result.donor}")
    out.append(
        f"    ends   : the 5' {ga.flap_len} nt match flap A and the 3' {gb.flap_len} nt match flap B"
    )
    out.append("")

    out.append("HOW TO PREPARE THE DONOR (from the published method, not yet from this design)")
    out += [f"  - {line}" for line in _DONOR_PREP]
    out.append("")

    out.append("Coordinates & invariants:")
    out.append(f"  insert site S      : {result.site_s}")
    out.append(
        f"  nicks              : A between window nt {n1} and {n1 + 1}; "
        f"B between window nt {n2} and {n2 + 1}"
    )
    out.append(
        f"  span replaced      : {n2 - n1} bp (window nt {n1 + 1}-{n2}), restored by the donor's arms "
        f"(5' arm {len(result.arm5)} bp, 3' arm {len(result.arm3)} bp)"
    )
    out.append(f"  deletion size      : {result.deletion_size} (INV-0 asserted zero)")
    out.append(f"  recode edits       : {result.recode_edits}")
    out.append(f"  PAM of guide A     : {_pam_sentence(ga)}")
    out.append(f"  PAM of guide B     : {_pam_sentence(gb)}")
    out.append(
        f"  flaps              : A {ga.flap_len} nt (longest run shared with the strand it displaces: "
        f"{ga.genomic_run}); B {gb.flap_len} nt ({gb.genomic_run})"
    )
    out.append(f"  pairs considered   : {result.pairs_considered}")
    out.append(f"  pairs walked       : {result.candidates_tried}")
    out.append(
        f"  predicted product  : {len(result.product)} nt (window {len(result.product) - len(result.cargo)} nt "
        f"+ cargo {len(result.cargo)} nt)"
    )
    out.append("")

    if result.changes is not None:
        out.append("Changes from your input:")
        out += ["  " + line for line in result.changes.lines(changes_file="recoded_pa_changes.csv")]
        out.append("")

    out.append("Advisory warnings:")
    if result.warnings:
        out += [f"  - {w}" for w in result.warnings]
    else:
        out.append("  (none)")
    out.append("")

    out.append("Provenance manifest:")
    for key in sorted(manifest):
        out.append(f"  {key}: {manifest[key]}")
    return "\n".join(out)


def _component_rows(result: RecodedPAResult) -> List[Tuple[str, str, str, str]]:
    rows: List[Tuple[str, str, str, str]] = []
    for side, g in (("A", result.guide_a), ("B", result.guide_b)):
        part = f"pegRNA {side}"
        end = "5'" if side == "A" else "3'"
        rows += [
            (
                f"spacer {side}",
                f"targets and nicks the {g.strand} strand at window nt {g.cut}",
                part,
                g.spacer,
            ),
            (f"scaffold {side}", "sgRNA scaffold (binds nCas9)", part, STANDARD_SGRNA_SCAFFOLD),
            (
                f"RTT {side} (flap {side})",
                f"reverse-transcription template: writes the donor's {end} {g.flap_len} nt",
                part,
                g.rtt,
            ),
            (f"PBS {side}", f"primer-binding site ({g.pbs_len} nt) at the nick", part, g.pbs),
            (f"PolyU {side}", "3' U tail (La/PE7)", part, g.tail),
        ]
    roles = {
        "arm5": f"restored copy of window nt {result.span[0] + 1}-{result.site_s} (recoded where coding)",
        "cargo": "insert / cargo payload, used exactly as given",
        "arm3": f"restored copy of window nt {result.site_s + 1}-{result.span[1]} (recoded where coding)",
    }
    for name, seq in result.donor_parts:
        rows.append((name, roles[name], "donor", seq))
    return rows


def render_pa_csv(result: RecodedPAResult, manifest: Optional[Mapping[str, Any]] = None) -> str:
    """One row per labelled component, in construct order (pegRNA A, pegRNA B, donor), with the same columns as
    the tj-PE CSV: ``index,component,role,part_of,sequence,length``. A trailing block of ``#`` rows carries the
    full sequences and the coordinates."""
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["index", "component", "role", "part_of", "sequence", "length"])
    for i, (component, role, part_of, seq) in enumerate(_component_rows(result), start=1):
        w.writerow([i, component, role, part_of, seq, len(seq)])
    meta = {
        "full_pegRNA_A": result.guide_a.pegRNA_sequence,
        "full_pegRNA_B": result.guide_b.pegRNA_sequence,
        "donor": result.donor,
        "donor_length": len(result.donor),
        "site_s": result.site_s,
        "nick_A": result.span[0],
        "nick_B": result.span[1],
        "deletion_size": result.deletion_size,
        "recode_edits": result.recode_edits,
    }
    for k, v in {**meta, **dict(manifest or {})}.items():
        w.writerow([f"# {k}", v])
    return buf.getvalue()


def render_pa_fasta(result: RecodedPAResult) -> str:
    """Three records: pegRNA A, pegRNA B, the donor. Sequences are wrapped at 80 and are exactly the ones in the
    report."""
    ga, gb = result.guide_a, result.guide_b
    n1, n2 = result.span
    records = [
        (
            f"pegRNA_A {_STRAND_WORD[ga.strand]}, nick between window nt {n1} and {n1 + 1}, "
            f"flap {ga.flap_len} nt, PBS {ga.pbs_len} nt, {len(ga.pegRNA_sequence)} nt",
            ga.pegRNA_sequence,
        ),
        (
            f"pegRNA_B {_STRAND_WORD[gb.strand]}, nick between window nt {n2} and {n2 + 1}, "
            f"flap {gb.flap_len} nt, PBS {gb.pbs_len} nt, {len(gb.pegRNA_sequence)} nt",
            gb.pegRNA_sequence,
        ),
        (
            f"donor arm5 + cargo + arm3, {len(result.donor)} nt, replaces window nt {n1 + 1}-{n2}",
            result.donor,
        ),
    ]
    return "\n".join(f">{head}\n{_fasta_wrap(seq)}" for head, seq in records) + "\n"


def _guide_payload(g: PAGuideDesign) -> Dict[str, Any]:
    return {
        "side": g.side,
        "strand": g.strand,
        "spacer": g.spacer,
        "pam": g.pam,
        "pam_status": g.pam_status,
        "product_pam": g.product_pam,
        "cut": g.cut,
        "flap_len": g.flap_len,
        "pbs_len": g.pbs_len,
        "pbs": g.pbs,
        "rtt": g.rtt,
        "tail": g.tail,
        "pegRNA": g.pegRNA_sequence,
        "pegRNA_length": len(g.pegRNA_sequence),
        "offtargets": g.offtargets,
        "on_target": g.on_target,
        "genomic_run": g.genomic_run,
        "nick_pairing_bp": g.nick_pairing_bp,
    }


def pa_gates_payload(
    result: RecodedPAResult, stand_in_gates: Sequence[str], mode: str
) -> Dict[str, Any]:
    """The JSON the web UI reads: the generic fields every design carries, plus a ``pa`` object that describes
    both guides, the donor and its parts, and the spans that were recoded or copied exactly."""
    parts: List[Dict[str, Any]] = []
    cursor = 0
    for name, seq in result.donor_parts:
        parts.append({"name": name, "start": cursor, "end": cursor + len(seq), "sequence": seq})
        cursor += len(seq)
    return {
        "strategy": "pa",
        "gate_results": [list(x) for x in result.gate_results],
        "recode_edits": result.recode_edits,
        "deletion_size": result.deletion_size,
        "candidates_tried": result.candidates_tried,
        "pairs_considered": result.pairs_considered,
        "warnings": list(result.warnings),
        "stand_in_gates": list(stand_in_gates),
        "scorer_mode": mode,
        "changes": result.changes.to_dict() if result.changes is not None else None,
        "pa": {
            "site_s": result.site_s,
            "span": list(result.span),
            "cargo_length": len(result.cargo),
            "donor_length": len(result.donor),
            "product_length": len(result.product),
            "guides": [_guide_payload(result.guide_a), _guide_payload(result.guide_b)],
            "donor": {
                "sequence": result.donor,
                "parts": parts,
                "flap_a": [0, result.guide_a.flap_len],
                "flap_b": [len(result.donor) - result.guide_b.flap_len, len(result.donor)],
            },
            "recode_spans": [list(s) for s in result.recode_spans],
            "locked_spans": [list(s) for s in result.locked_spans],
            "simulated": result.simulated,
        },
    }
