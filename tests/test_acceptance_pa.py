"""The acceptance check for a real-model Prime Assembly run must be able to FAIL.

``scripts/analysis/acceptance_pa.py`` reads the files a reader would order from and re-derives the design's
invariants from them. A check that cannot fail proves nothing, so this runs it on a fake-mode design (which is valid
except that every gate is a stand-in) and then corrupts one thing at a time: each corruption must turn exactly the
check that is about it red.
"""

import importlib.util
import json
import shutil
from pathlib import Path

import pytest
from pa_helpers import make_cargo, make_locus

from pegrna_pipeline.hpc.recoded_pa_driver import PA_FILES, design_pa

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "analysis" / "acceptance_pa.py"
S = 150
CARGO = make_cargo(120)
#: A fake-mode run is a correct design whose gates are all stand-ins, so exactly these checks cannot pass.
FAKE_ONLY = {
    "no gate was served by a stand-in",
    "real models were used",
    "the report carries no stand-in banner",
    "docking_a is a verdict, not an advisory",
    "docking_b is a verdict, not an advisory",
}


@pytest.fixture(scope="module")
def acceptance():
    spec = importlib.util.spec_from_file_location("acceptance_pa", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def config():
    loc = make_locus(
        [("top", 123, "AGG"), ("top", 147), ("bottom", 164), ("bottom", 190)],
        n=300,
        site_s=S,
        seed=7,
        coding_end=S,
    )
    return dict(cargo=CARGO, genome=loc.genome, site_s=S, cds_mask=list(loc.cds_mask))


@pytest.fixture(scope="module")
def good(tmp_path_factory, config):
    out = tmp_path_factory.mktemp("good")
    design_pa(**config, mode="fake", output_dir=out)
    return out


@pytest.fixture
def run(tmp_path, good, acceptance, config):
    """Copy the good output, let the test corrupt it, then return {check text: ok}."""
    work = tmp_path / "out"
    shutil.copytree(good, work)

    def go(cfg=None):
        return acceptance.run_checks(cfg or config, work)

    go.dir = work
    return go


def failed(results):
    """The full text of every check that failed."""
    return [text for ok, text in results if not ok]


def has(texts, prefix):
    """True if some failed check starts with ``prefix``."""
    return any(text.startswith(prefix) for text in texts)


def test_a_correct_design_passes_every_check_except_the_ones_only_real_models_can_pass(run):
    results = run()
    assert len(results) >= 24
    got = failed(results)
    assert all(any(text.startswith(prefix) for prefix in FAKE_ONLY) for text in got), got
    assert all(has(got, prefix) for prefix in FAKE_ONLY), got


def test_it_stops_early_and_says_why_when_a_file_is_missing(run):
    (run.dir / PA_FILES["fasta"]).unlink()
    results = run()
    assert len(results) == 1 and results[0][0] is False and "missing" in results[0][1]


def test_a_failed_run_is_not_an_acceptance(run):
    gates = json.loads((run.dir / PA_FILES["gates"]).read_text())
    gates.update(failed=True, error="NoValidPADesign: nothing survived")
    (run.dir / PA_FILES["gates"]).write_text(json.dumps(gates))
    results = run()
    assert len(results) == 2 and results[1][0] is False and "NoValidPADesign" in results[1][1]


def _edit_gates(run, **changes):
    path = run.dir / PA_FILES["gates"]
    gates = json.loads(path.read_text())
    gates.update(changes)
    path.write_text(json.dumps(gates))
    return gates


def test_a_stand_in_gate_is_caught(run):
    _edit_gates(run, stand_in_gates=["docking"])
    assert has(failed(run()), "no gate was served by a stand-in")


def test_a_fake_scorer_mode_is_caught(run):
    assert has(failed(run()), "real models were used")  # fake mode: the good fixture is fake


def test_a_missing_gate_is_caught(run):
    gates = json.loads((run.dir / PA_FILES["gates"]).read_text())
    _edit_gates(run, gate_results=[g for g in gates["gate_results"] if g[0] != "splice"])
    assert has(failed(run()), "every gate ran, in pipeline order")


def test_gates_out_of_order_are_caught(run):
    gates = json.loads((run.dir / PA_FILES["gates"]).read_text())["gate_results"]
    _edit_gates(run, gate_results=list(reversed(gates)))
    assert has(failed(run()), "every gate ran, in pipeline order")


def test_an_advisory_docking_is_not_a_verdict(run):
    gates = json.loads((run.dir / PA_FILES["gates"]).read_text())["gate_results"]
    advisory = [
        [k, "advisory only: interface ipTM 0.7 at 669 nt" if k == "docking_a" else v]
        for k, v in gates
    ]
    _edit_gates(run, gate_results=advisory)
    assert has(failed(run()), "docking_a is a verdict, not an advisory")


def _edit_fasta(run, edit):
    path = run.dir / PA_FILES["fasta"]
    records = [r for r in path.read_text().split(">") if r.strip()]
    heads = [r.splitlines()[0] for r in records]
    seqs = ["".join(r.splitlines()[1:]) for r in records]
    seqs = edit(seqs)
    path.write_text("".join(f">{h}\n{s}\n" for h, s in zip(heads, seqs)))


def _swap(i, old, new):
    def edit(seqs):
        assert old in seqs[i]
        seqs[i] = seqs[i].replace(old, new, 1)
        return seqs

    return edit


def test_a_wrong_base_in_the_donor_is_caught_twice_over(run):
    """The donor in the FASTA no longer equals the parts the JSON lists, and the ordered pegRNAs no longer
    reproduce it: two independent reasons, both reported."""
    base = {"A": "C", "C": "A", "G": "T", "T": "G"}

    def edit(seqs):
        seqs[2] = base[seqs[2][0]] + seqs[2][1:]
        return seqs

    _edit_fasta(run, edit)
    got = failed(run())
    assert has(got, "the donor is arm5 + cargo + arm3")
    assert has(got, "simulate_pa rejected the ordered pegRNAs and donor")


def test_pegrnas_given_in_the_wrong_order_are_caught_twice_over(run):
    _edit_fasta(run, lambda s: [s[1], s[0], s[2]])
    got = failed(run())
    assert has(got, "pegRNA A starts with its spacer") and has(
        got, "pegRNA B starts with its spacer"
    )
    assert has(got, "simulate_pa rejected the ordered pegRNAs and donor")


def test_a_flap_that_does_not_match_the_donor_is_caught(run):
    gates = json.loads((run.dir / PA_FILES["gates"]).read_text())
    rtt = gates["pa"]["guides"][0]["rtt"]

    def edit(seqs):
        assert rtt in seqs[0]
        seqs[0] = seqs[0].replace(rtt, rtt[:-1] + ("A" if rtt[-1] != "A" else "C"), 1)
        return seqs

    _edit_fasta(run, edit)
    got = failed(run())
    assert has(got, "simulate_pa rejected the ordered pegRNAs and donor")
    assert not has(
        got, "the donor is arm5 + cargo + arm3"
    ), "only the pegRNA is wrong, not the donor"


def test_an_overlong_pegrna_is_caught(run):
    _edit_fasta(run, lambda s: [s[0][:-3] + "A" * 40 + "UUU", s[1], s[2]])
    got = failed(run())
    assert has(got, "pegRNA A is") and "at most 166" in next(
        g for g in got if g.startswith("pegRNA A is")
    )
    assert not has(got, "pegRNA B is")


def test_a_pegrna_without_its_polyu_tail_is_caught(run):
    _edit_fasta(run, lambda s: [s[0][:-3] + "GGG", s[1], s[2]])
    got = failed(run())
    assert has(got, "pegRNA A starts with its spacer and ends in UUU")
    assert not has(got, "pegRNA B starts with its spacer")


def test_a_terminator_in_the_extension_is_caught_even_across_the_junction(run):
    path = run.dir / PA_FILES["gates"]
    gates = json.loads(path.read_text())
    g = gates["pa"]["guides"][0]
    g["rtt"], g["pbs"] = "GCGCGCGCTT", "TTGCGCGC"  # neither part holds TTTT; the transcript does
    _edit_gates(run, pa=gates["pa"])
    assert has(failed(run()), "pegRNA A: no TTTT")


def test_the_scaffolds_own_run_of_t_is_not_mistaken_for_a_terminator(run):
    assert not any("no TTTT" in g for g in failed(run())), "the standard scaffold holds GTTTTAG"


def test_a_terminator_in_the_spacer_is_caught(run):
    """The spacer holds TTTT; the FASTA and the JSON agree on it, so only the terminator check can object."""
    gates = json.loads((run.dir / PA_FILES["gates"]).read_text())
    spacer = gates["pa"]["guides"][0]["spacer"]
    gates["pa"]["guides"][0]["spacer"] = "TTTT" + spacer[4:]
    _edit_gates(run, pa=gates["pa"])
    _edit_fasta(run, lambda s: ["TTTT" + s[0][4:], s[1], s[2]])
    got = failed(run())
    assert has(got, "pegRNA A: no TTTT")
    assert not has(
        got, "pegRNA A starts with its spacer"
    ), "the FASTA and the JSON agree on the spacer"


def test_a_guide_that_is_not_where_it_says_it_is_fails_cleanly_and_not_with_a_crash(run):
    gates = json.loads((run.dir / PA_FILES["gates"]).read_text())
    gates["pa"]["guides"][1]["cut"] += 1
    _edit_gates(run, pa=gates["pa"])
    got = failed(run())
    assert has(got, "both guides are real NGG guides of the window, at the stated cuts")
    assert not has(got, "simulate_pa"), "the simulation is skipped when the guides cannot be found"


def test_a_cargo_that_is_not_at_the_insertion_site_is_caught(run, config):
    moved = dict(config, site_s=config["site_s"] + 1)
    assert has(failed(run(moved)), "the cargo sits at the insertion site S")


def test_a_cargo_that_is_not_the_configs_is_caught(run, config):
    other = dict(config, cargo=config["cargo"][:-1] + ("A" if config["cargo"][-1] != "A" else "C"))
    got = failed(run(other))
    assert has(got, "the cargo in the donor is exactly the config's cargo")


def test_a_change_missing_from_the_changes_file_is_caught(run):
    path = run.dir / PA_FILES["changes"]
    lines = path.read_text().splitlines()
    assert len(lines) > 2, "the design must have recoded something for this to mean anything"
    path.write_text("\n".join([lines[0]] + lines[2:]) + "\n")
    assert has(failed(run()), "every changed base is listed")


def test_a_listed_change_that_is_not_in_the_product_is_caught(run):
    path = run.dir / PA_FILES["changes"]
    lines = path.read_text().splitlines()
    ghost = lines[1].split(",")
    ghost[2] = str(int(ghost[2]) + 1)  # one position along: listed, but nothing changed there
    path.write_text("\n".join(lines + [",".join(ghost)]) + "\n")
    assert has(failed(run()), "every changed base is listed")


def test_a_change_that_alters_the_protein_is_caught(run, config):
    """Same files, but the config's reading frame is shifted by one: the recoded bases are no longer synonymous."""
    mask = list(config["cds_mask"])
    shifted = dict(config, cds_mask=[False] + mask[:-1])
    got = failed(run(shifted))
    assert has(got, "every changed base is protein-coding and synonymous")


def test_a_change_outside_the_coding_sequence_is_caught(run, config):
    shifted = dict(config, cds_mask=[False] * len(config["cds_mask"]))
    assert has(failed(run(shifted)), "every changed base is protein-coding and synonymous")


def test_a_deletion_is_caught(run):
    """The donor replaces five bases more than its arms restore: the product is five bases short of the window plus
    the cargo, and that is exactly what INV-0 forbids."""
    gates = json.loads((run.dir / PA_FILES["gates"]).read_text())
    n1, n2 = gates["pa"]["span"]
    gates["pa"]["span"] = [n1, n2 + 5]
    _edit_gates(run, pa=gates["pa"])
    got = failed(run())
    assert has(got, "the product is the window plus the cargo")
    assert has(got, "nothing is deleted")
    assert has(
        got, "simulate_pa rebuilds the product"
    ), "the nicks no longer bound the span the donor replaced"


def test_a_report_without_the_standing_caveat_is_caught(run):
    path = run.dir / PA_FILES["report"]
    path.write_text(path.read_text().replace("extension of Prime Assembly", "work on assembly"))
    assert has(failed(run()), "the report carries the standing caveat")


def test_a_report_that_carries_the_stand_in_banner_is_caught(run):
    path = run.dir / PA_FILES["report"]
    path.write_text("*** NOT A VALIDATED DESIGN - STAND-IN GATES IN USE ***\n" + path.read_text())
    assert has(failed(run()), "the report carries no stand-in banner")


def test_the_command_line_entry_prints_each_check_and_exits_nonzero_on_failure(
    run, acceptance, tmp_path, config, capsys
):
    cfg_path = tmp_path / "cfg.json"
    cfg_path.write_text(json.dumps(config))
    code = acceptance.main(["acceptance_pa.py", str(cfg_path), str(run.dir)])
    out = capsys.readouterr().out
    assert code == 1 and "PASS  all five output files are present" in out and "FAIL  " in out
    assert "checks passed" in out
    assert acceptance.main(["acceptance_pa.py"]) == 2
