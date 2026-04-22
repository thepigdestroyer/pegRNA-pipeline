"""Document and guard the PRIDICT ``Editing_Position`` indexing convention.

PRIDICT 2.0 emits ``Editing_Position`` as a **0-indexed** offset from the
start of the input window. The fixtures we ship use literal ``50`` — that
same integer should appear verbatim on the parsed :class:`PegRNADesign`
objects. If PRIDICT ever switches to 1-indexed output, *this* test will
fail first and make the assumption obvious.
"""

from __future__ import annotations

from pathlib import Path

from pegrna_pipeline.hpc.pridict_parser import parse_pegrna_csv


FIXTURE = Path(__file__).parent / "fixtures" / "pridict" / "ctcf_tjpe_v1_pegRNA_Pridict_full.csv"


def test_editing_position_is_preserved_as_zero_indexed():
    designs = parse_pegrna_csv(FIXTURE)
    for d in designs:
        assert d.editing_position == 50, (
            "Fixture emits literal 50 for every row; Editing_Position must "
            "round-trip unchanged (0-indexed convention)."
        )


def test_editing_position_is_python_int():
    """A float here would silently break integer arithmetic downstream."""
    designs = parse_pegrna_csv(FIXTURE)
    for d in designs:
        assert isinstance(d.editing_position, int)
        assert not isinstance(d.editing_position, bool)
