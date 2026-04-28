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

# PRIDICT 2.0's nicking-guides CSV currently spells the PAM-disrupt
# column ``Nicking-PAMdisrupt``; pre-2.0 builds used ``PAM_disrupted``.
# Look for either, in priority order.
_PAM_DISRUPT_COLS = ("Nicking-PAMdisrupt", "PAM_disrupted")


PathLike = Union[str, Path]


def _require_columns(df: pd.DataFrame, required: List[str], source: PathLike) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{source}: missing required column(s): {', '.join(missing)}")


def _coerce_bool(val: object) -> bool:
    """Map a PRIDICT 2.0 boolean-ish value to Python ``bool``.

    PRIDICT emits these flags as strings ("PE3b" / "No_PE3b" for PE3b,
    "Nicking_PAM_disrupt" / "No_nicking_PAM_disrupt" for the PAM-disrupt
    flag) rather than canonical booleans, so we map any value with a
    leading ``no_`` (case-insensitive) to ``False`` and any non-empty
    other string to ``True``.
    """
    if isinstance(val, bool):
        return val
    if isinstance(val, (int, float)) and not isinstance(val, bool):
        return bool(val)
    s = str(val).strip().lower()
    if not s or s in {"false", "0", "no", "f", "n", "na", "nan"}:
        return False
    if s.startswith("no_"):
        return False
    if s in {"true", "1", "yes", "t", "y"}:
        return True
    # Anything else non-empty (e.g. "PE3b", "Nicking_PAM_disrupt") → True.
    return True


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
    # PRIDICT 2.0 emits the assembled oligo as ``pegRNA`` (column 32 in
    # the v2.0 output schema). Older builds used ``Scaffold_Optimized``
    # or ``pegRNA_full`` — fall back through the chain.
    full_col = next(
        (c for c in ("pegRNA", "Scaffold_Optimized", "pegRNA_full") if c in df.columns),
        None,
    )

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
                pegrna_full=str(row[full_col]) if full_col else "",
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
    pam_col = next((c for c in _PAM_DISRUPT_COLS if c in df.columns), None)

    designs: List[NgRNADesign] = []
    for _, row in df.iterrows():
        designs.append(
            NgRNADesign(
                spacer=str(row["Nicking-Protospacer"]),
                distance_to_edit=int(row["Nicking-Position-to-edit"]),
                target_strand=str(row["Target_Strand"]),  # "Fw" or "Rv"
                deepcas9_score=float(row["DeepCas9score"]),
                pe3b=_coerce_bool(row["PE3b"]) if has_pe3b else False,
                pam_disrupt=_coerce_bool(row[pam_col]) if pam_col else False,
            )
        )
    return designs
