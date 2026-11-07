"""Per-component annotations for TJ-PE and PE7 results.

Takes Plan 2's TJPEResult / PE7Result plus a GeneRecord and
produces an ordered list of ComponentAnnotation records that the report
writers render mode-aware.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import List

from pegrna_pipeline.hpc.data_types import GeneRecord
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
                f"{len(result.ngrna_spacer)}-nt nicking-sgRNA spacer. "
                f"Directs Cas9-H840A to nick the opposite strand "
                f"{abs(result.ngrna_distance_to_edit)} bp from the edit "
                f"(signed {result.ngrna_distance_to_edit:+d})."
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
                f"{len(result.ngrna_spacer)}-nt nicking-sgRNA spacer. "
                f"Directs Cas9-H840A to nick the opposite ({result.ngrna_strand}) "
                f"strand {abs(result.ngrna_distance_to_edit)} bp from the edit "
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
