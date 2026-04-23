"""Parse PRIDICT 2.0 CSV output files into typed dataclasses.

PRIDICT emits two CSVs per input sequence:

* ``<name>_pegRNA_Pridict_full.csv`` — all candidate pegRNAs with PE2 scores.
* ``<name>_nicking_guides.csv`` — candidate nicking guides for PE3/PE3b.

The column naming is *not* consistent between the two files: the pegRNA CSV
spells it ``Target-Strand`` (hyphen) whereas the ngRNA CSV uses
``Target_Strand`` (underscore). The ``_PEGRNA_REQUIRED`` and
``_NGRNA_REQUIRED`` lists document the exact columns this module needs —
missing any of them is a hard error so callers never silently process
partial data.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Union

import pandas as pd

from pegrna_pipeline.hpc.data_types import NgRNADesign, PegRNADesign

_PEGRNA_REQUIRED = [
    "sequence_name",
    "Spacer-Sequence",
    "PBSrevcomp",
    "RTrevcomp",
    "PBSlength",
    "RTlength",
    "Target-Strand",
    "Editing_Position",
    "PRIDICT2_0_editing_Score_deep_HEK",
]

_NGRNA_REQUIRED = [
    "Nicking-Protospacer",
    "Nicking-Position-to-edit",
    "Target_Strand",
    "DeepCas9score",
]


PathLike = Union[str, Path]


def _require_columns(df: pd.DataFrame, required: List[str], source: PathLike) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{source}: missing required column(s): {', '.join(missing)}")


def _coerce_bool(val: object) -> bool:
    if isinstance(val, bool):
        return val
    if isinstance(val, (int, float)) and not isinstance(val, bool):
        return bool(val)
    s = str(val).strip().lower()
    return s in {"true", "1", "yes", "t", "y"}


# ---------------------------------------------------------------- parsers


def parse_pegrna_csv(path: PathLike) -> List[PegRNADesign]:
    """Parse a ``*_pegRNA_Pridict_full.csv`` into a list of PegRNADesign.

    Raises
    ------
    ValueError
        If any required column is missing.
    """
    path = Path(path)
    df = pd.read_csv(path)
    _require_columns(df, _PEGRNA_REQUIRED, path)

    # Optional columns — fall back to sentinel values if absent.
    has_k562 = "PRIDICT2_0_editing_Score_deep_K562" in df.columns
    has_full = "pegRNA_full" in df.columns

    designs: List[PegRNADesign] = []
    for _, row in df.iterrows():
        designs.append(
            PegRNADesign(
                sequence_name=str(row["sequence_name"]),
                spacer=str(row["Spacer-Sequence"]),
                pbs=str(row["PBSrevcomp"]),
                rt=str(row["RTrevcomp"]),
                pbs_length=int(row["PBSlength"]),
                rt_length=int(row["RTlength"]),
                target_strand=str(row["Target-Strand"]),  # "Fw" or "Rv"
                editing_position=int(row["Editing_Position"]),
                pe2_hek293t=float(row["PRIDICT2_0_editing_Score_deep_HEK"]),
                pe2_k562=float(row["PRIDICT2_0_editing_Score_deep_K562"]) if has_k562 else 0.0,
                pegrna_full=str(row["pegRNA_full"]) if has_full else "",
            )
        )
    return designs


def parse_ngrna_csv(path: PathLike) -> List[NgRNADesign]:
    """Parse a ``*_nicking_guides.csv`` into a list of NgRNADesign.

    Raises
    ------
    ValueError
        If any required column is missing.
    """
    path = Path(path)
    df = pd.read_csv(path)
    _require_columns(df, _NGRNA_REQUIRED, path)

    has_pe3b = "PE3b" in df.columns
    has_pam_disrupt = "PAM_disrupted" in df.columns

    designs: List[NgRNADesign] = []
    for _, row in df.iterrows():
        designs.append(
            NgRNADesign(
                spacer=str(row["Nicking-Protospacer"]),
                distance_to_edit=int(row["Nicking-Position-to-edit"]),
                target_strand=str(row["Target_Strand"]),  # "Fw" or "Rv"
                deepcas9_score=float(row["DeepCas9score"]),
                pe3b=_coerce_bool(row["PE3b"]) if has_pe3b else False,
                pam_disrupt=_coerce_bool(row["PAM_disrupted"]) if has_pam_disrupt else False,
            )
        )
    return designs
