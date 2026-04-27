"""Per-component annotations for TJ-PE, PA, and PE7 results.

Takes Plan 2's TJPEResult / PAResult / PE7Result plus a GeneRecord and
produces an ordered list of ComponentAnnotation records that the report
writers render mode-aware.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import List

from pegrna_pipeline.hpc.data_types import GeneRecord
from pegrna_pipeline.hpc.pa_workflow import PAResult
from pegrna_pipeline.hpc.pe7_workflow import PE7Result
from pegrna_pipeline.hpc.tjpe_workflow import STANDARD_SGRNA_SCAFFOLD, TJPEResult


class ComponentSource(Enum):
    PRIDICT2 = "PRIDICT2"
    GENOMIC_CONTEXT = "genomic_context"
    CODON_ENUMERATION = "codon_enumeration"
    CODON_OPTIMISED_BULK = "codon_optimised_bulk"
    MOTIF_LIBRARY = "motif_library"
    INVARIANT = "invariant"
    CRISPR_SCORE = "crisprScore"


@dataclass(frozen=True)
class ComponentAnnotation:
    name: str
    sequence: str
    length_bp: int
    functional_role: str
    source: ComponentSource
    design_notes: str = ""


def annotate_tjpe_result(result: TJPEResult, gene: GeneRecord) -> List[ComponentAnnotation]:
    """Build per-component annotations for a TJ-PE design.

    Ordered: spacer, sgRNA_scaffold, insert, PBS1, RC-PBS2, ngRNA_spacer,
    ngRNA_scaffold.
    """
    annotations: List[ComponentAnnotation] = []
    annotations.append(
        ComponentAnnotation(
            name="spacer",
            sequence=result.spacer,
            length_bp=len(result.spacer),
            functional_role=(
                f"{len(result.spacer)}-nt guide that directs the prime editor "
                f"(Cas9-H840A nickase) to the {gene.gene_name} locus. Creates a "
                f"single-strand nick on the target strand."
            ),
            source=ComponentSource.PRIDICT2,
            design_notes=(
                f"PE2 {result.pegrna_pe2_cell_line.upper()} = " f"{result.pegrna_pe2_score:.2f}%"
            ),
        )
    )
    annotations.append(
        ComponentAnnotation(
            name="sgRNA_scaffold",
            sequence=STANDARD_SGRNA_SCAFFOLD,
            length_bp=len(STANDARD_SGRNA_SCAFFOLD),
            functional_role=(
                "Invariant hairpin recognised by Cas9; identical for every pegRNA "
                "and supplied by the cloning backbone."
            ),
            source=ComponentSource.INVARIANT,
        )
    )
    annotations.append(
        ComponentAnnotation(
            name="insert",
            sequence=result.full_insert_dna,
            length_bp=len(result.full_insert_dna),
            functional_role=(
                "Full DNA inserted by TJ-PE. The first 8 amino acids (24 bp) "
                "were chosen from codon variants ranked by PRIDICT 2.0 PE2 "
                "score; remainder is human-codon-optimised."
            ),
            source=ComponentSource.CODON_ENUMERATION,
        )
    )
    annotations.append(
        ComponentAnnotation(
            name="PBS1",
            sequence=result.pbs1,
            length_bp=len(result.pbs1),
            functional_role=(
                f"{len(result.pbs1)}-nt primer binding site that hybridises to "
                f"the free 3' end exposed by the pegRNA nick, priming the first "
                f"reverse-transcription extension."
            ),
            source=ComponentSource.PRIDICT2,
        )
    )
    annotations.append(
        ComponentAnnotation(
            name="RC-PBS2",
            sequence=result.rc_pbs2,
            length_bp=result.rc_pbs2_length,
            functional_role=(
                f"Reverse complement of PBS2 ({result.rc_pbs2_length} nt). Sits at "
                f"the 3' end of the TJ-pegRNA. After RT extension, the new top-strand "
                f"flap carries PBS2 at its 3' end, which hybridises to the bottom "
                f"strand 3' end exposed by the ngRNA to enable template jumping. "
                f"Derived from the mRNA-frame position "
                f"{result.ngrna_nick_in_mrna}."
            ),
            source=ComponentSource.GENOMIC_CONTEXT,
            design_notes=(
                "VERIFY against a published Zheng 2023 supplementary TJ-pegRNA "
                "before first synthesis (spec v1.1 checkpoint)."
            ),
        )
    )
    annotations.append(
        ComponentAnnotation(
            name="ngRNA_spacer",
            sequence=result.ngrna_spacer,
            length_bp=len(result.ngrna_spacer),
            functional_role=(
                f"20-nt nicking-sgRNA spacer. Directs Cas9-H840A to nick the "
                f"opposite strand {abs(result.ngrna_distance_to_edit)} bp from "
                f"the edit (signed {result.ngrna_distance_to_edit:+d})."
            ),
            source=ComponentSource.PRIDICT2,
            design_notes=f"DeepCas9 score = {result.ngrna_deepcas9_score:.2f}",
        )
    )
    annotations.append(
        ComponentAnnotation(
            name="ngRNA_scaffold",
            sequence=STANDARD_SGRNA_SCAFFOLD,
            length_bp=len(STANDARD_SGRNA_SCAFFOLD),
            functional_role="Invariant sgRNA scaffold (same as pegRNA scaffold).",
            source=ComponentSource.INVARIANT,
        )
    )
    return annotations


def annotate_pe7_result(result: PE7Result, gene: GeneRecord) -> List[ComponentAnnotation]:
    """Build per-component annotations for a PE7 (standard PE3-strategy) design.

    Ordered: spacer, sgRNA_scaffold, RT_template, PBS, ngRNA_spacer,
    ngRNA_scaffold. The insert is *encoded inside* the RT template — there is
    no separate ``insert`` component as in TJ-PE.
    """
    return [
        ComponentAnnotation(
            name="spacer",
            sequence=result.spacer,
            length_bp=len(result.spacer),
            functional_role=(
                f"{len(result.spacer)}-nt guide that directs the prime editor "
                f"(Cas9-H840A nickase) to the {gene.gene_name} locus on the "
                f"{result.pegrna_strand} strand. Creates the pegRNA-side nick."
            ),
            source=ComponentSource.PRIDICT2,
            design_notes=(
                f"PE2 {result.pegrna_pe2_cell_line.upper()} = " f"{result.pegrna_pe2_score:.2f}%"
            ),
        ),
        ComponentAnnotation(
            name="sgRNA_scaffold",
            sequence=STANDARD_SGRNA_SCAFFOLD,
            length_bp=len(STANDARD_SGRNA_SCAFFOLD),
            functional_role=(
                "Invariant hairpin recognised by Cas9; identical for every "
                "pegRNA and supplied by the cloning backbone."
            ),
            source=ComponentSource.INVARIANT,
        ),
        ComponentAnnotation(
            name="RT_template",
            sequence=result.rt,
            length_bp=len(result.rt),
            functional_role=(
                f"{len(result.rt)}-nt reverse-transcription template. Encodes "
                f"the {len(result.full_insert_dna)}-bp insert plus a few "
                f"flanking nucleotides for hybridisation."
            ),
            source=ComponentSource.PRIDICT2,
        ),
        ComponentAnnotation(
            name="PBS",
            sequence=result.pbs,
            length_bp=len(result.pbs),
            functional_role=(
                f"{len(result.pbs)}-nt primer binding site. Hybridises to the "
                f"3' end exposed by the pegRNA nick, priming RT extension."
            ),
            source=ComponentSource.PRIDICT2,
        ),
        ComponentAnnotation(
            name="ngRNA_spacer",
            sequence=result.ngrna_spacer,
            length_bp=len(result.ngrna_spacer),
            functional_role=(
                f"20-nt nicking-sgRNA spacer. Directs Cas9-H840A to nick the "
                f"opposite ({result.ngrna_strand}) strand "
                f"{abs(result.ngrna_distance_to_edit)} bp from the edit "
                f"(signed {result.ngrna_distance_to_edit:+d}). Drives PE3-style "
                f"repair bias."
            ),
            source=ComponentSource.PRIDICT2,
            design_notes=f"DeepCas9 score = {result.ngrna_deepcas9_score:.2f}",
        ),
        ComponentAnnotation(
            name="ngRNA_scaffold",
            sequence=STANDARD_SGRNA_SCAFFOLD,
            length_bp=len(STANDARD_SGRNA_SCAFFOLD),
            functional_role="Invariant sgRNA scaffold (same as pegRNA scaffold).",
            source=ComponentSource.INVARIANT,
        ),
    ]


def annotate_pa_result(
    result: PAResult,
    gene: GeneRecord,
    score_cell_line: str = "k562",
) -> List[ComponentAnnotation]:
    """Build per-component annotations for a Prime Assembly design.

    ``score_cell_line`` controls which PRIDICT PE2 column the per-pegRNA
    notes display ("k562" or "hek293t").
    """
    from pegrna_pipeline.hpc.data_types import _score_attr

    score_attr = _score_attr(score_cell_line)
    up = result.winning_pair.upstream
    dn = result.winning_pair.downstream
    cl_label = score_cell_line.upper()
    up_score = float(getattr(up, score_attr))
    dn_score = float(getattr(dn, score_attr))
    return [
        ComponentAnnotation(
            name="upstream_pegRNA_spacer",
            sequence=up.spacer,
            length_bp=len(up.spacer),
            functional_role=(
                f"20-nt spacer for the upstream pegRNA. Targets the "
                f"{up.target_strand} strand of {gene.gene_name} at editing "
                f"position {up.editing_position}. Installs the upstream flap."
            ),
            source=ComponentSource.PRIDICT2,
            design_notes=f"PE2 {cl_label} = {up_score:.2f}%",
        ),
        ComponentAnnotation(
            name="upstream_sgRNA_scaffold",
            sequence=STANDARD_SGRNA_SCAFFOLD,
            length_bp=len(STANDARD_SGRNA_SCAFFOLD),
            functional_role="Invariant scaffold for the upstream pegRNA.",
            source=ComponentSource.INVARIANT,
        ),
        ComponentAnnotation(
            name="upstream_PBS",
            sequence=up.pbs,
            length_bp=up.pbs_length,
            functional_role="PBS for the upstream pegRNA; primes RT on the upstream nick.",
            source=ComponentSource.PRIDICT2,
        ),
        ComponentAnnotation(
            name="upstream_RTT",
            sequence=up.rt,
            length_bp=up.rt_length,
            functional_role=(
                "RT template for the upstream pegRNA. Encodes the upstream flap "
                "(= 5' HA of the linear donor)."
            ),
            source=ComponentSource.PRIDICT2,
        ),
        ComponentAnnotation(
            name="upstream_flap_DNA",
            sequence=result.upstream_flap_dna,
            length_bp=len(result.upstream_flap_dna),
            functional_role=(
                "What the upstream pegRNA writes into the genome. Identical to "
                "the linear donor's 5' homology arm; they hybridise during "
                "in-cell assembly."
            ),
            source=ComponentSource.CODON_ENUMERATION,
        ),
        ComponentAnnotation(
            name="downstream_pegRNA_spacer",
            sequence=dn.spacer,
            length_bp=len(dn.spacer),
            functional_role=(
                f"20-nt spacer for the downstream pegRNA. Targets the "
                f"{dn.target_strand} strand of {gene.gene_name} at editing "
                f"position {dn.editing_position}. Installs the downstream flap."
            ),
            source=ComponentSource.PRIDICT2,
            design_notes=f"PE2 {cl_label} = {dn_score:.2f}%",
        ),
        ComponentAnnotation(
            name="downstream_sgRNA_scaffold",
            sequence=STANDARD_SGRNA_SCAFFOLD,
            length_bp=len(STANDARD_SGRNA_SCAFFOLD),
            functional_role="Invariant scaffold for the downstream pegRNA.",
            source=ComponentSource.INVARIANT,
        ),
        ComponentAnnotation(
            name="downstream_PBS",
            sequence=dn.pbs,
            length_bp=dn.pbs_length,
            functional_role="PBS for the downstream pegRNA.",
            source=ComponentSource.PRIDICT2,
        ),
        ComponentAnnotation(
            name="downstream_RTT",
            sequence=dn.rt,
            length_bp=dn.rt_length,
            functional_role=(
                "RT template for the downstream pegRNA. Encodes the downstream "
                "flap (= 3' HA of the donor)."
            ),
            source=ComponentSource.PRIDICT2,
        ),
        ComponentAnnotation(
            name="downstream_flap_DNA",
            sequence=result.downstream_flap_dna,
            length_bp=len(result.downstream_flap_dna),
            functional_role="What the downstream pegRNA writes into the genome.",
            source=ComponentSource.CODON_ENUMERATION,
        ),
        ComponentAnnotation(
            name="linear_donor_DNA",
            sequence=result.linear_donor_dna,
            length_bp=len(result.linear_donor_dna),
            functional_role=(
                "Linear dsDNA donor (gBlock). 5' matches upstream_flap_DNA "
                "verbatim; 3' matches downstream_flap_DNA verbatim. Middle is "
                "codon-optimised. Gibson-like in-cell assembly bridges the two "
                "installed flaps with this donor."
            ),
            source=ComponentSource.CODON_OPTIMISED_BULK,
            design_notes=(
                "Order as IDT gBlock or Twist Gene Fragment. Consider AZD-7648 "
                "(DNA-PK inhibitor) co-treatment for ~2x efficiency gain "
                "(Liu 2025)."
            ),
        ),
    ]
