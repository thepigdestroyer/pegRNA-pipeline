"""Human-frequency codon tiebreak in recode_arm (ADR 0009, within min-edits).

recode_arm breaks microhomology with the fewest synonymous edits; among the
CHANGED codon positions it now prefers the highest-frequency human synonymous
codon that still keeps the longest common run <= max_run. Unchanged positions
stay byte-exact (min-edits + no-op preserved).
"""

from Bio.Seq import Seq

from pegrna_pipeline.hpc.codon_expander import codon_frequency, codons_for_aa
from pegrna_pipeline.hpc.microhomology import longest_common_run
from pegrna_pipeline.hpc.recode_arm import _prefer_human_codons, recode_arm


def _codons(dna, frame=0):
    pad = (3 - frame) % 3
    body = dna[pad:]
    return [body[i : i + 3] for i in range(0, len(body) - len(body) % 3, 3)]


def _max_freq_codon(aa):
    return max(codons_for_aa(aa), key=lambda c: codon_frequency(aa, c))


# --- _prefer_human_codons unit ----------------------------------------------


def test_prefer_leaves_unchanged_positions_byte_exact():
    seq = "CTCGAA"
    assert _prefer_human_codons(seq, original=seq, genome="AAAAAAAAA", frame=0, max_run=3) == seq


def test_prefer_upgrades_changed_codon_to_max_freq_when_feasible():
    original = "CTC"  # Leu
    recoded = "TTA"  # Leu, low-frequency -> a changed position
    out = _prefer_human_codons(recoded, original=original, genome="AAAAAAAAA", frame=0, max_run=3)
    assert str(Seq(out).translate()) == "L"
    assert out == _max_freq_codon("L")


def test_prefer_preserves_protein_and_frame_prefix():
    # frame=1: the single pad base is outside any codon and must be untouched.
    original = "GCTCCTC"  # 'G' pad + CTC CTC after? build explicitly below
    # Simpler: pad base 'G', then two Leu codons; recoded changed both to low-freq.
    original = "G" + "CTC" + "CTC"
    recoded = "G" + "TTA" + "TTG"
    out = _prefer_human_codons(recoded, original=original, genome="AAAAAAAAAA", frame=1, max_run=3)
    assert out[0] == "G"  # pad base untouched
    assert str(Seq(out[1:]).translate()) == "LL"


def test_prefer_respects_run_constraint():
    original = "CTC"
    recoded = "TTA"
    top = _max_freq_codon("L")
    genome = top * 3  # a run of the top codon makes it infeasible at max_run=2
    out = _prefer_human_codons(recoded, original=original, genome=genome, frame=0, max_run=2)
    assert str(Seq(out).translate()) == "L"
    assert longest_common_run(out, genome) <= 2
    feasible = [c for c in codons_for_aa("L") if longest_common_run(c, genome) <= 2]
    best_feasible = max(feasible, key=lambda c: codon_frequency("L", c))
    assert codon_frequency("L", out) == codon_frequency("L", best_feasible)


# --- recode_arm end-to-end ---------------------------------------------------


def test_recode_arm_changed_codons_are_locally_freq_optimal():
    cds_span = "CTC" * 6  # Leu x6, self-repetitive -> forces a microhomology-breaking recode
    max_run = 4
    res = recode_arm(cds_span, genome=cds_span, frame=0, max_run=max_run)
    assert res is not None
    assert str(Seq(res.recoded).translate()) == str(Seq(cds_span).translate())
    assert longest_common_run(res.recoded, cds_span) <= max_run
    assert 0 < res.n_edits <= len(cds_span)

    recoded = res.recoded
    ocodons, rcodons = _codons(cds_span), _codons(recoded)
    for i, (oc, rc) in enumerate(zip(ocodons, rcodons)):
        if rc == oc:
            continue
        aa = str(Seq(rc).translate())
        for alt in codons_for_aa(aa):
            if codon_frequency(aa, alt) <= codon_frequency(aa, rc):
                continue
            s = 3 * i
            trial = recoded[:s] + alt + recoded[s + 3 :]
            assert (
                longest_common_run(trial, cds_span) > max_run
            ), f"changed codon {rc} at {i} should have been upgraded to higher-freq {alt}"


def test_recode_arm_no_change_when_below_run_stays_byte_exact():
    # Unrelated genome -> no microhomology to break -> no edits -> byte-exact,
    # so the human-codon pass must NOT introduce changes (min-edits preserved).
    cds_span = "CTCGAAAAG"
    res = recode_arm(cds_span, genome="TTTTTTTTTTTT", frame=0, max_run=6)
    assert res is not None
    assert res.n_edits == 0
    assert res.recoded == cds_span
