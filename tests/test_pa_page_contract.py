"""The page and the pipeline agree on what a failed Prime Assembly run said.

The page reads the walk's diagnostic lines (``web/frontend/src/pipeline.ts``) to say which step rejected a design and
how often, and checks the form with limits written in ``web/frontend/src/pa.ts``. Nothing in the walk's own tests would
notice if a message were reworded or a limit changed: the page would mislabel a failure, or let through a request the
server then refuses. This closes that gap from the Python side. It reads the page's OWN source (its patterns and its
numbers), builds every kind of rejection with the REAL walk (stand-in scorers), and checks that

* every line the walk writes about a pair is claimed by exactly one step, and by the right one, so the failure view
  neither drops a reason nor blames the wrong step;
* the form's defaults and limits are the engine's, so what the page accepts the server accepts.

The browser-side half (what the page does with the lines) is ``scripts/qa-frontend-logic.mjs``.
"""

import dataclasses
import importlib
import re
import sys
from pathlib import Path

import pytest
from pa_helpers import make_cargo
from test_recoded_pa_workflow import (
    CARGO,
    PLANTS,
    CrossRegionFold,
    DonorsInTheCargo,
    FlagBlockedImportance,
    S,
    _corrupting,
    _dirty_for,
    _missense,
    _plan_args,
    run,
)

from pegrna_pipeline.hpc.complex_dock import DockFake
from pegrna_pipeline.hpc.pa_geometry import MIN_NICK_SEPARATION_BP, find_guides
from pegrna_pipeline.hpc.pa_simulator import PAMismatch
from pegrna_pipeline.hpc.recoded_pa_driver import design_pa
from pegrna_pipeline.hpc.recoded_pa_workflow import NoValidPADesign, plan_pairs, validate_pa_request
from pegrna_pipeline.hpc.splice_screen import SpliceFake

REPO = Path(__file__).resolve().parents[1]
FRONTEND = REPO / "web" / "frontend" / "src"


# ---- reading the page's source -------------------------------------------------------------------------------


def _source(name):
    return (FRONTEND / name).read_text()


def _regexes(text):
    """The regex literals in ``text``, compiled. JavaScript's ``/source/flags`` is Python's source for these."""
    return [
        re.compile(source, re.IGNORECASE if "i" in flags else 0)
        for source, flags in re.findall(r"/((?:\\.|[^/\\\n])+)/([a-z]*)", text)
    ]


def page_steps():
    """``{step id: [failure patterns]}`` for Prime Assembly's steps, in the page's order."""
    text = _source("pipeline.ts")
    body = text[text.index("export const PA_PIPELINE") : text.index("/** The cascade to draw")]
    steps = {}
    for block in body.split("\n  {\n")[1:]:
        step = re.search(r'id: "(\w+)"', block).group(1)
        failure = re.search(r"failure: \[(.*?)\],\n", block, re.S)
        steps[step] = _regexes(failure.group(1)) if failure else []
    return steps


def page_pair_line():
    """The page's own pattern for 'a line about one pair'."""
    match = re.search(r"const ABOUT_ONE_PAIR = (/.*/[a-z]*);", _source("pipeline.ts"))
    return _regexes(match.group(1))[0]


def _object(text, name):
    """``{key: number-or-string}`` of the first ``export const NAME = { ... }`` in the page's ``pa.ts``."""
    body = re.search(rf"export const {name} = \{{(.*?)\}}", text, re.S).group(1)
    return {
        k: (int(v) if v.isdigit() else v.strip('"'))
        for k, v in re.findall(r'(\w+): ("?\w+"?),', body)
    }


def claimed_by(line, steps):
    return [step for step, patterns in steps.items() if any(p.search(line) for p in patterns)]


INFORMATIONAL = re.compile(r"\d+ PAM-in pair\(s\) from \d+ guide\(s\); \d+ pair\(s\) remain")


# ---- building each rejection with the real walk -------------------------------------------------------------


@pytest.fixture(scope="module")
def locus():
    from pa_helpers import make_locus

    return make_locus(PLANTS, n=300, site_s=S, seed=7, coding_end=S)


class _Boom:
    """Raises from whichever scorer method the walk calls, so the 'scorer error' lines are written."""

    def __init__(self, method):
        setattr(self, method, self._raise)

    def _raise(self, *a, **k):
        raise RuntimeError("it exploded")


class _FoldsOnly:
    """Confidently pairs the RTT with the PBS for exactly these pegRNAs, and nothing for the rest."""

    def __init__(self, bad):
        self.bad = set(bad)

    def pair_probs(self, seq):
        return {(100, len(seq) - 6): 0.999} if seq in self.bad else {}


class _DocksOnly(DockFake):
    """Docks every pegRNA well except these."""

    def __init__(self, bad, error=False):
        super().__init__(apo_iptm=0.85, candidate_iptm=0.82)
        self.bad, self.error = set(bad), error

    def interface_iptm(self, pegrna, cas9):
        if pegrna in self.bad:
            if self.error:
                raise RuntimeError("it exploded")
            return 0.1
        return 0.82


def _pegrnas(locus, side):
    plans = plan_pairs(**_plan_args(locus)).plans
    return {
        c.pegrna
        for plan in plans
        for c in (plan.candidates_a if side == "A" else plan.candidates_b)
    }


def _walk(locus, **over):
    with pytest.raises(NoValidPADesign) as exc:
        run(locus, **over)
    return list(exc.value.diagnostics)


def _plan(locus, **over):
    return plan_pairs(**_plan_args(locus, **over)).diagnostics


def _wf():
    from pegrna_pipeline.hpc import recoded_pa_workflow

    return recoded_pa_workflow


def scenarios(locus, monkeypatch):
    """``[(name, expected steps, diagnostic lines)]``: one rejection of each kind the walk can write.

    A scorer that fails every pegRNA is judged on both sides of each pair, so its lines belong to both pegRNAs' steps;
    the scorers that fail only one side are there to show each step alone.
    """
    wf = _wf()
    every = {g.spacer for g in find_guides(locus.genome)}
    found = []

    def add(name, expected, lines):
        found.append((name, {expected} if isinstance(expected, str) else set(expected), lines))

    add(
        "the insertion site is rejected",
        "importance",
        _walk(locus, importance_scorer=FlagBlockedImportance()),
    )
    add(
        "the importance scorer fails",
        "importance",
        _walk(locus, importance_scorer=_Boom("metrics")),
    )
    add("no pair of guides fits", "pairs", _plan(locus, min_span=95, max_span=100))

    real = wf.assess_spacer
    monkeypatch.setattr(
        wf,
        "assess_spacer",
        lambda s: dataclasses.replace(real(s), failures=("GC content 25% is outside 30-75%",)),
    )
    add("every spacer breaks a sequence rule", "rules", _plan(locus))
    monkeypatch.setattr(wf, "assess_spacer", real)

    add(
        "every guide fails the off-target screen",
        "offtarget",
        _plan(locus, offtarget_scorer=_dirty_for(every)),
    )

    real_arm = wf.build_arm
    monkeypatch.setattr(
        wf, "build_arm", lambda g, lo, hi, m, r: None if lo == 123 else real_arm(g, lo, hi, m, r)
    )
    add("an arm cannot be recoded", "recode", _plan(locus))
    monkeypatch.setattr(wf, "build_arm", real_arm)

    _corrupting(monkeypatch, 123, _missense)
    add("an arm changes the protein", "inv0", _plan(locus))
    monkeypatch.setattr(wf, "build_arm", real_arm)
    _corrupting(monkeypatch, 123, lambda sense: sense + "A")
    add("an arm changes the length", "inv0", _plan(locus))
    monkeypatch.setattr(wf, "build_arm", real_arm)

    add(
        "a terminator in every flap",
        "flaps",
        _plan(locus, cargo="ACGTACGTAC" + "AAAA" + make_cargo(106)),
    )
    add(
        "the donor's two ends are alike",
        "flaps",
        _plan(locus, cargo=locus.genome[170:182] + make_cargo(108)),
    )

    real_side, real_pa = wf.simulate_side, wf.simulate_pa

    def side_broken(**kw):
        raise PAMismatch("flap does not match the donor")

    def pa_broken(**kw):
        raise PAMismatch("the product is not the intended one")

    monkeypatch.setattr(wf, "simulate_side", side_broken)
    add("a pegRNA fails the strand check while being built", "strand", _plan(locus))
    monkeypatch.setattr(wf, "simulate_side", real_side)
    monkeypatch.setattr(wf, "simulate_pa", pa_broken)
    add("the chosen pair fails the final strand check", "strand", _walk(locus))
    monkeypatch.setattr(wf, "simulate_pa", real_pa)

    only = lambda side: _pegrnas(locus, side)  # noqa: E731 - a short name for the tables below
    add(
        "only pegRNA A folds badly",
        "structure_a",
        _walk(locus, full_structure_scorer=_FoldsOnly(only("A"))),
    )
    add(
        "only pegRNA B folds badly",
        "structure_b",
        _walk(locus, full_structure_scorer=_FoldsOnly(only("B"))),
    )
    add(
        "every pegRNA folds badly",
        {"structure_a", "structure_b"},
        _walk(locus, full_structure_scorer=CrossRegionFold()),
    )
    add(
        "the folder fails",
        {"structure_a", "structure_b"},
        _walk(locus, full_structure_scorer=_Boom("pair_probs")),
    )
    add("a cryptic splice site", "cryptic", _walk(locus, maxent_scorer=DonorsInTheCargo(CARGO)))
    add("the cryptic-splice scorer fails", "cryptic", _walk(locus, maxent_scorer=_Boom("score5")))
    add("SpliceAI rejects", "splice", _walk(locus, splice_scorer=SpliceFake(0.5)))
    add("SpliceAI fails", "splice", _walk(locus, splice_scorer=_Boom("delta")))
    add(
        "every pegRNA binds Cas9 too weakly",
        "docking_a",
        _walk(locus, dock_scorer=DockFake(0.85, 0.10)),
    )
    add(
        "only pegRNA B binds Cas9 too weakly",
        "docking_b",
        _walk(locus, dock_scorer=_DocksOnly(only("B"))),
    )
    add("docking A fails", "docking_a", _walk(locus, dock_scorer=_DocksOnly(only("A"), error=True)))
    add("docking B fails", "docking_b", _walk(locus, dock_scorer=_DocksOnly(only("B"), error=True)))
    return found


@pytest.fixture(scope="module")
def built(locus):
    """Every scenario, built once: each is a full walk, and several tests read them."""
    with pytest.MonkeyPatch.context() as monkeypatch:
        return scenarios(locus, monkeypatch)


# ---- the failure view reads what the walk wrote --------------------------------------------------------------


def test_every_rejection_the_walk_writes_is_filed_under_the_right_step(built):
    steps = page_steps()
    problems = []
    for name, expected, lines in built:
        reasons = [ln for ln in lines if not INFORMATIONAL.fullmatch(ln)]
        if not reasons:
            problems.append(f"{name}: the walk wrote no reason at all")
        reached = set()
        for line in reasons:
            claimed = claimed_by(line, steps)
            reached.update(claimed)
            if len(claimed) != 1 or not set(claimed) <= expected:
                problems.append(
                    f"{name}: wanted one of {sorted(expected)}, the page files it under {claimed}: {line[:140]}"
                )
        if reasons and reached != expected:
            problems.append(f"{name}: the lines reach {sorted(reached)}, not {sorted(expected)}")
    assert not problems, "\n".join(problems)


def test_the_scenarios_cover_every_step_that_can_reject_a_pair_or_a_run(built):
    covered = set().union(*(expected for _, expected, _ in built))
    can_reject = {step for step, patterns in page_steps().items() if patterns}
    assert covered == can_reject, (can_reject - covered, covered - can_reject)


#: The steps whose rejections are about the whole run or a whole guide, not about one pair of guides.
RUN_WIDE_STEPS = {"importance", "pairs", "rules", "offtarget"}


def test_the_page_tells_a_line_about_one_pair_from_a_line_about_the_run(built):
    """The page counts an unclaimed line as 'other' only if it is about a pair. Which lines those are is decided by
    what the step is about, not by how the walk happens to label a pair, so a new label fails here.
    """
    pair_line = page_pair_line()
    for name, expected, lines in built:
        about_a_pair = not expected <= RUN_WIDE_STEPS
        for line in lines:
            if INFORMATIONAL.fullmatch(line):
                continue
            assert bool(pair_line.search(line)) == about_a_pair, (name, line[:120])


def test_the_summary_line_carries_the_numbers_the_page_counts_from(locus, built):
    """One line reports both a number of guides that broke the sequence rules and a number the scan rejected."""
    text = _source("pipeline.ts")
    wanted = {m for m in re.findall(r"countFrom: (/.*/),", text)}
    assert wanted == {r"/(\d+) failed the sequence rules/", r"/(\d+) failed the off-target screen/"}
    lines = [ln for name, _, lines in built if name.startswith("every") for ln in lines]
    counted = {
        step: [int(m.group(1)) for ln in lines for m in [re.search(pattern, ln)] if m]
        for step, pattern in (
            ("rules", r"(\d+) failed the sequence rules"),
            ("offtarget", r"(\d+) failed the off-target screen"),
        )
    }
    assert counted["rules"] and counted["offtarget"], counted
    guides = len({g.spacer for g in find_guides(locus.genome)})
    assert all(n <= guides for ns in counted.values() for n in ns)


def test_no_step_of_the_page_is_dead(built):
    """Every pattern the page has for a step is one the walk's real lines match: none is left over from a wording
    that no longer exists."""
    steps = page_steps()
    lines = [ln for _, _, found in built for ln in found]
    for step, patterns in steps.items():
        for pattern in patterns:
            assert any(
                pattern.search(ln) for ln in lines
            ), f"{step}: nothing the walk writes matches {pattern.pattern}"


# ---- the form's numbers are the engine's ----------------------------------------------------------------------


def _valid(locus, **over):
    args = dict(cargo=CARGO, genome=locus.genome, site_s=S, cds_mask=list(locus.cds_mask))
    args.update(over)
    return validate_pa_request(**args)


def _refused(locus, **over):
    with pytest.raises((ValueError, NoValidPADesign)):
        _valid(locus, **over)


def test_the_forms_defaults_are_the_engines(locus):
    import inspect

    defaults = _object(_source("pa.ts"), "PA_DEFAULTS")
    params = inspect.signature(design_pa).parameters
    assert set(defaults) == {"flap_length", "min_span", "max_span", "max_offtargets", "max_pairs"}
    for name, value in defaults.items():
        assert str(params[name].default) == value, (name, params[name].default, value)
    app = _source("App.tsx")
    assert f'max_run: "{params["max_run"].default}"' in app
    assert 'pbs_lengths: "' + ",".join(str(n) for n in params["pbs_lengths"].default) + '"' in app


def test_what_the_page_accepts_at_its_limits_the_engine_accepts_and_just_outside_it_refuses(locus):
    limits = _object(_source("pa.ts"), "PA_LIMITS")
    _valid(locus, flap_length=limits["flapMin"], flap_min=limits["flapMin"])
    _valid(locus, flap_length=limits["flapMax"], flap_max=limits["flapMax"])
    _refused(locus, flap_length=limits["flapMin"] - 1, flap_min=limits["flapMin"] - 1)
    _refused(locus, flap_length=limits["flapMax"] + 1, flap_max=limits["flapMax"] + 1)
    assert limits["spanFloor"] == MIN_NICK_SEPARATION_BP
    _valid(locus, min_span=limits["spanFloor"])
    _refused(locus, min_span=limits["spanFloor"] - 1)
    _valid(locus, pbs_lengths=[limits["pbsMin"], limits["pbsMax"]])
    _refused(locus, pbs_lengths=[limits["pbsMin"] - 1])
    _refused(locus, pbs_lengths=[limits["pbsMax"] + 1])
    # the page is stricter than the engine here, never looser: what it lets through is accepted
    _valid(locus, min_span=limits["spanFloor"], max_span=limits["spanCap"])
    _valid(locus, max_offtargets=limits["maxOfftargets"])


@pytest.fixture()
def app_module(monkeypatch, tmp_path):
    pytest.importorskip("fastapi")
    pytest.importorskip("httpx")
    monkeypatch.setenv("PEGRNA_REPO", str(REPO))
    monkeypatch.setenv("PEGRNA_RUNS_DIR", str(tmp_path / "runs"))
    monkeypatch.delenv("PEGRNA_PA_MAX_PAIRS", raising=False)
    if str(REPO / "web") not in sys.path:
        sys.path.insert(0, str(REPO / "web"))
    sys.modules.pop("app", None)
    module = importlib.import_module("app")
    return importlib.reload(module)


def test_the_pages_caps_on_pairs_and_pbs_lengths_are_the_servers(app_module):
    limits = _object(_source("pa.ts"), "PA_LIMITS")
    assert limits["maxPairs"] == app_module.PA_MAX_PAIRS
    assert limits["pbsCount"] == app_module.PA_MAX_PBS_LENGTHS


def test_the_strategy_names_and_files_the_page_uses_are_the_servers(app_module):
    from pegrna_pipeline.hpc.recoded_pa_driver import PA_FILES

    panel = _source("components/ResultPanel.tsx")
    both = panel + _source("components/PaResult.tsx")
    assert (
        f'? "{app_module.STRATEGIES["pa"]["prefix"]}" : "{app_module.STRATEGIES["tjpe"]["prefix"]}"'
        in panel
    )
    # the report and the components are named from the prefix, the other two outright
    assert PA_FILES["report"] == f'{app_module.STRATEGIES["pa"]["prefix"]}_report.txt'
    assert PA_FILES["components"] == f'{app_module.STRATEGIES["pa"]["prefix"]}_components.csv'
    assert "`${prefix}_report.txt`" in panel and "`${prefix}_components.csv`" in panel
    assert f'"{PA_FILES["fasta"]}"' in both and f'"{PA_FILES["changes"]}"' in both
    assert set(app_module.STRATEGIES) == {"tjpe", "pa"}
    assert 'export type Strategy = "tjpe" | "pa";' in _source("types.ts")
