/**
 * Logic checks for the page's pure TypeScript: no browser, no server, milliseconds.
 *
 *   node scripts/qa-frontend-logic.mjs
 *
 * It imports web/frontend/src/{pa,pipeline}.ts directly, so what is tested is the source the page is built from.
 * Node must be able to strip types: that is the default from Node 23.6 (22.18); on 22.6 to 22.17 run it with
 * `node --experimental-strip-types`. The files have no JSX and import only types, which is what keeps this possible
 * (api.ts reads `import.meta.env`, which only Vite provides, and is covered by scripts/qa-ui-quick.mjs instead).
 *
 * What is here, because a browser run cannot reach it:
 *  - the Prime Assembly form rules, at their boundaries, against the numbers the server enforces;
 *  - which settings each strategy sends, and that neither sends the other's;
 *  - which strategy is advised for an insert of a given length, and that advice never forces a switch;
 *  - how the failure view reads the pipeline's diagnostic lines: which step each line belongs to, how many
 *    rejections each step accounts for, and that a line no step claims is counted and never dropped;
 *  - the live state of the Prime Assembly steps while a run streams, which a one-shot mocked response cannot show;
 *  - reading the progress stream, and giving up on one that has gone silent, with a fake reader on a timer.
 */
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const SRC = join(dirname(fileURLToPath(import.meta.url)), "..", "web", "frontend", "src");
const pa = await import(join(SRC, "pa.ts"));
const pl = await import(join(SRC, "pipeline.ts"));
const st = await import(join(SRC, "stream.ts"));

let failures = 0;
const fail = (m) => {
  console.log(`  FAIL  ${m}`);
  failures += 1;
  process.exitCode = 1;
};
const ok = (m) => console.log(`  PASS  ${m}`);
const check = (cond, name, detail = "") => (cond ? ok(name) : fail(`${name}${detail ? `: ${detail}` : ""}`));
const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
const show = (x) => JSON.stringify(x);

// ------------------------------------------------------------------------------------------------------------
console.log("\nPRIME ASSEMBLY FORM RULES");
{
  const good = { ...pa.PA_DEFAULTS, pbs_lengths: "13,15,17", max_run: "8" };
  const errs = (over = {}, cargo = 499) => pa.validatePa({ ...good, ...over }, cargo);

  check(same(errs(), {}), "the defaults and a 499 nt insert are valid", show(errs()));
  check(
    same(pa.PA_DEFAULTS, { flap_length: "35", min_span: "30", max_span: "100", max_offtargets: "0", max_pairs: "30" }),
    "the defaults are the engine's: flap 35, span 30-100, no off-targets, 30 pairs",
  );

  for (const bad of ["29", "51", "", "abc", "35.5", "-35", "3e1"])
    check("flap_length" in errs({ flap_length: bad }), `a flap length of "${bad}" is refused`);
  for (const fine of ["30", "35", "50"]) check(!("flap_length" in errs({ flap_length: fine })), `a flap length of ${fine} is accepted`);

  // the cargo rule: at least two flaps, so the two flaps cannot overlap
  check("cargo" in errs({}, 69) && /70 nt/.test(errs({}, 69).cargo), "an insert of 69 nt is refused for flap 35, naming 70", show(errs({}, 69)));
  check(!("cargo" in errs({}, 70)), "an insert of exactly 70 nt is accepted for flap 35");
  check("cargo" in errs({ flap_length: "50" }, 99) && !("cargo" in errs({ flap_length: "50" }, 100)), "the minimum follows the flap length (100 at flap 50)");
  check(!("cargo" in errs({}, 0)), "an empty insert is not called too short (it is called missing elsewhere)");
  check(!("cargo" in errs({ flap_length: "20" }, 30)) && "flap_length" in errs({ flap_length: "20" }, 30), "a bad flap length is reported as that, not as a cargo problem");

  check("min_span" in errs({ min_span: "22" }) && !("min_span" in errs({ min_span: "23" })), "the span floor is 23 bp, and 22 is refused");
  check("max_span" in errs({ min_span: "60", max_span: "50" }), "a maximum span below the minimum is refused");
  check(!("max_span" in errs({ min_span: "50", max_span: "50" })) && !("min_span" in errs({ min_span: "50", max_span: "50" })), "a maximum span equal to the minimum is accepted (one fixed distance)");
  check("max_span" in errs({ max_span: "201" }) && !("max_span" in errs({ max_span: "200" })), "the maximum span is capped at 200");
  check("min_span" in errs({ min_span: "x" }) && "max_span" in errs({ max_span: "" }), "a span that is not a number is refused");

  check(!("max_offtargets" in errs({ max_offtargets: "0" })) && !("max_offtargets" in errs({ max_offtargets: "10" })), "off-target allowances of 0 and 10 are accepted");
  for (const bad of ["-1", "11", "1.5", ""]) check("max_offtargets" in errs({ max_offtargets: bad }), `an off-target allowance of "${bad}" is refused`);

  check(!("max_pairs" in errs({ max_pairs: "1" })) && !("max_pairs" in errs({ max_pairs: "60" })), "1 and 60 pairs to try are accepted");
  check("max_pairs" in errs({ max_pairs: "0" }) && "max_pairs" in errs({ max_pairs: "61" }), "0 and 61 pairs to try are refused (the server caps at 60)");

  for (const fine of ["13,15,17", "13, 15 ,17", "6", "30", "8,9,10,11,12"])
    check(!("pbs_lengths" in errs({ pbs_lengths: fine })), `PBS lengths "${fine}" are accepted`);
  for (const bad of ["", "13,x", "5", "31", "8,9,10,11,12,13", "13;15", "1.5"])
    check("pbs_lengths" in errs({ pbs_lengths: bad }), `PBS lengths "${bad}" are refused`);
  check(same(pa.parsePbs("13, 15 ,17"), [13, 15, 17]) && pa.parsePbs("13,x") === null && pa.parsePbs("") === null, "PBS lengths parse as whole numbers or not at all");
  check(pa.parsePbs("13,,15") !== null && same(pa.parsePbs("13,,15"), [13, 15]), "an empty entry between commas is skipped");

  check("max_run" in errs({ max_run: "0" }) && "max_run" in errs({ max_run: "x" }) && !("max_run" in errs({ max_run: "8" })), "the longest allowed repeat must be a whole number from 1");
  check(Object.keys(errs({ flap_length: "1", min_span: "1", max_span: "999", max_offtargets: "99", max_pairs: "0", pbs_lengths: "", max_run: "0" })).length === 7, "every bad field is reported at once, not just the first");
}

// ------------------------------------------------------------------------------------------------------------
console.log("\nWHAT EACH STRATEGY SENDS");
{
  const form = { ...pa.PA_DEFAULTS, pbs_lengths: "13, 15", max_run: "7" };
  const base = {
    gene: "CTCF", transcript: "ENST1", uniprot_acc: "P49711", cargo: "ACGT", genome: "GGCC", site_s: 2, cds_mask: [true, true, false, false],
    rtt_aa: "L", pegrna_strand: "Rv", arm_side: "5prime", ngrna_distance_min: 40, ngrna_distance_max: 90, ngrna_max_offtargets: 0,
    spacer: "A".repeat(20), het_snps: [{}], het_edit: true, protein_length: 727, structure_threshold: 0.9,
  };
  const frozen = JSON.stringify(base);
  const out = pa.toPaConfig(base, form);

  check(JSON.stringify(base) === frozen, "building a Prime Assembly request does not change the form's own config");
  // written out here, not read from the code: a check that iterated the code's own list would shrink with it
  const tjOnly = ["rtt_aa", "pegrna_strand", "arm_side", "ngrna_distance_min", "ngrna_distance_max", "ngrna_max_offtargets", "spacer", "het_snps", "het_edit"];
  check(tjOnly.every((k) => k in base) && tjOnly.every((k) => !(k in out)), "a Prime Assembly request carries none of tj-PE's settings", show(Object.keys(out)));
  check(same([...pa.TJPE_ONLY_KEYS].sort(), [...tjOnly].sort()), "and the list the page uses for that is exactly those nine");
  check(["gene", "transcript", "uniprot_acc", "cargo", "genome", "site_s", "cds_mask", "protein_length", "structure_threshold"].every((k) => same(out[k], base[k])), "and keeps the shared ones untouched (including the transcript the splicing check needs)");
  check(out.flap_length === 35 && out.min_span === 30 && out.max_span === 100 && out.max_offtargets === 0 && out.max_pairs === 30 && out.max_run === 7, "the settings are numbers, not the form's text", show(out));
  check(same(out.pbs_lengths, [13, 15]), "the PBS lengths are a list of numbers");
  check(["flap_length", "min_span", "max_span", "max_offtargets", "max_pairs"].every((k) => typeof out[k] === "number"), "every numeric setting really is a number");

  const tj = pa.toTjpeConfig({ ...base, flap_length: 35, flap_min: 30, flap_max: 50, min_span: 30, max_span: 100, max_offtargets: 0, max_pairs: 30, max_run: 8, pbs_lengths: [13] });
  check(pa.PA_ONLY_KEYS.every((k) => !(k in tj)), "a tj-PE request carries none of Prime Assembly's settings", show(Object.keys(tj)));
  check(tj.max_run === 8 && same(tj.pbs_lengths, [13]) && tj.pegrna_strand === "Rv" && tj.spacer === "A".repeat(20), "and keeps the shared and its own ones, max_run and pbs_lengths among them");
  check(same([...pa.PA_ONLY_KEYS].sort(), ["flap_length", "flap_max", "flap_min", "max_offtargets", "max_pairs", "max_span", "min_span"]), "the Prime Assembly-only keys are exactly the settings only that engine reads");
  check(!pa.PA_ONLY_KEYS.includes("max_run") && !pa.PA_ONLY_KEYS.includes("pbs_lengths"), "max_run and pbs_lengths are shared by both strategies");
}

// ------------------------------------------------------------------------------------------------------------
console.log("\nWHICH STRATEGY FITS");
{
  const adv = pa.strategyAdvice;
  check(adv("pa", 0, 35) === null && adv("tjpe", 0, 35) === null, "an empty insert gets no advice");
  check(adv("pa", 70, 35) === null && adv("pa", 499, 35) === null, "Prime Assembly with an insert of 70 nt or more gets no advice");
  const warn = adv("pa", 69, 35);
  check(warn?.tone === "warn" && warn.switchTo === "tjpe" && /70 nt/.test(warn.text), "Prime Assembly with 69 nt is warned, naming 70, and offered the other way", show(warn));
  check(adv("pa", 99, 50)?.tone === "warn" && adv("pa", 100, 50) === null, "the warning threshold follows the flap length");
  check(adv("tjpe", 24, 35) === null && adv("tjpe", 69, 35) === null, "a tag or a short insert written inside the pegRNA gets no advice");
  const either = adv("tjpe", 85, 35);
  check(either?.tone === "info" && either.switchTo === "pa" && /either works/.test(either.text), "between 70 and 99 nt both are said to work, and neither is pushed", show(either));
  const longer = adv("tjpe", 100, 35);
  check(longer?.tone === "info" && /very long pegRNA/.test(longer.text) && longer.switchTo === "pa", "from 100 nt, writing it inside the pegRNA is advised against, with the reason");
  check(/roughly 649 nt/.test(adv("tjpe", 499, 35).text), "the pegRNA length it warns about is the insert plus about 150 nt (649 for the 499 nt cassette)");
  check(adv("tjpe", 100, 50)?.text.includes("very long") && !adv("tjpe", 99, 50)?.text.includes("very long"), "at flap 50 the 100 nt line still marks the same boundary");
  check([adv("pa", 69, 35), adv("tjpe", 85, 35), adv("tjpe", 499, 35)].every((a) => a.action && a.switchTo), "every piece of advice offers a button, and the page never switches by itself");

  const f = (status, pam = "AAG") => pa.pamFate({ pam_status: status, product_pam: pam });
  check(f("broken").tone === "ok" && f("broken-by-recode").tone === "ok", "a broken PAM is good news: the guide cannot cut the product");
  check(f("weakened", "AGA").tone === "warn" && f("weakened", "AGA").text.includes("AGA"), "a weakened PAM warns, and names the PAM the product has");
  check(f("intact-noncoding").tone === "warn" && f("intact-no-synonymous-change").tone === "warn", "a regenerated PAM warns that the guide can cut the product again");
  check(/non-coding/.test(f("intact-noncoding").text) && /no synonymous change/.test(f("intact-no-synonymous-change").text), "and says why, differently for each cause");
  check(pa.isStandingCaveat("A donor that restores the DNA between the nicks is this project's extension") && !pa.isStandingCaveat("guide A's PAM can only be weakened"), "the standing caveat is told apart from a finding about the design");
}

// ------------------------------------------------------------------------------------------------------------
console.log("\nTHE STEPS");
{
  const ids = (steps) => steps.map((s) => s.id);
  check(same(ids(pl.PIPELINE), ["input", "offtarget", "enum", "recode", "inv0", "cryptic", "splice", "importance", "structure", "docking", "report"]), "tj-PE's steps are exactly what they were");
  check(same(ids(pl.PA_PIPELINE), ["importance", "pairs", "rules", "offtarget", "recode", "inv0", "flaps", "strand", "structure_a", "structure_b", "cryptic", "splice", "docking_a", "docking_b", "report"]), "Prime Assembly's steps are in the order the walk does them");
  check(pl.pipelineFor("tjpe") === pl.PIPELINE && pl.pipelineFor("pa") === pl.PA_PIPELINE, "each strategy draws its own steps");
  check(new Set(ids(pl.PA_PIPELINE)).size === pl.PA_PIPELINE.length, "no step id is used twice");
  check(pl.PA_PIPELINE.every((s) => s.why.length > 40), "every step says what it protects against");
  check(pl.PA_PIPELINE.filter((s) => s.failure).every((s) => s.unit), "every step that can reject says what it rejects");

  // tj-PE's failure location is the page's original, by keyword
  const idxOf = (steps, id) => steps.findIndex((s) => s.id === id);
  check(pl.findBlockingStep("rejected by the off-target gate", "") === idxOf(pl.PIPELINE, "offtarget"), "tj-PE: a failure that mentions off-target is located at the off-target step");
  check(pl.findBlockingStep("nothing recognisable", "") === -1, "tj-PE: a failure that names nothing is located nowhere");
  check(pl.findBlockingStep("x", "REPORT: cryptic splice site") === idxOf(pl.PIPELINE, "cryptic"), "tj-PE: the report text is searched too");
}

// ------------------------------------------------------------------------------------------------------------
console.log("\nREADING A FAILED PRIME ASSEMBLY RUN");
{
  const steps = pl.PA_PIPELINE;
  const stepsFor = (line) => steps.filter((s) => s.failure?.some((re) => re.test(line))).map((s) => s.id);
  const summary = (n, m, extra) => `${n} PAM-in pair(s) from ${m} guide(s)${extra}; 0 pair(s) remain`;

  // each kind of line the walk writes, with the step it belongs to. These strings are the engine's own; the Python
  // test tests/test_pa_page_contract.py builds the same kinds with the real walk and checks them against the same
  // patterns, so a rewording on either side fails a test instead of silently mislabelling a failure.
  const kinds = [
    ["importance", "insertion site S residue(s) [725, 726] structurally important (buried/packed/special) -- inserting a tag/linker here risks disrupting the protein (this does not depend on the guide pair, so no pair was tried)"],
    ["importance", "importance scorer error: AlphaFold-DB unreachable"],
    ["pairs", "0 PAM-in pairs: the window has 4 top-strand and 4 bottom-strand NGG guides, but none make a top-strand guide cutting at or 5' of S and a bottom-strand guide cutting at or 3' of S that are 30-100 bp apart with their PAMs clear of each other"],
    ["rules", summary(24, 16, "; 16 failed the sequence rules (GC content 25% is outside 30-75%; contains TTTT: RNA polymerase III terminates there, so a U6-expressed guide is truncated)")],
    ["offtarget", summary(24, 16, "; 16 failed the off-target screen (each needs a unique on-target and at most 0 off-target site(s) within the gate mismatch tier)")],
    ["recode", "pair 123/164: arm partition failed (no synonymous CDS solution or non-monotone arm)"],
    ["inv0", "pair 123/164: INV-0: donor arm 5' is not protein-silent"],
    ["inv0", "pair 123/164: the product is 421 nt, expected 420"],
    ["flaps", "pair 147/190: guide A (cut 147): no usable extension: TTTT (a Pol III terminator) in the RTT or PBS at every flap length 30-50"],
    ["flaps", "pair 147/190: guide B (cut 190): no usable extension: the PBS window leaves the genome window, or the flap is longer than the donor"],
    ["flaps", "pair 147/190: the donor's two ends share a run of 12 bp (more than 8): the two flaps could pair with each other instead of with the donor"],
    ["strand", "pair 70/149: pegRNA A (PBS 13): PA strand mismatch: flap B does not match the 3' end of the donor"],
    ["strand", "pair 70/149: PA strand mismatch: the predicted product differs from the intended one"],
    ["structure_a", "pair 101/164: pegRNA A (PBS 13, flap 35): cross-region pairing (ViennaRNA): spacer[17..17] C x RTT[34..34] G (1 bp, p >= 0.999)"],
    ["structure_a", "pair 101/164: pegRNA A (PBS 13, flap 35): structure scorer error: ViennaRNA missing"],
    ["structure_b", "pair 101/164: pegRNA B (PBS 15, flap 35): cross-region pairing (ViennaRNA): spacer[17..18] AC x RTT[33..34] GT (2 bp, p >= 0.996)"],
    ["structure_b", "pair 101/164: pegRNA B (PBS 15, flap 35): structure scorer error: ViennaRNA missing"],
    ["cryptic", "pair 70/149: cryptic splice: new donor site MaxEnt 11.2 >= 8"],
    ["cryptic", "pair 70/149: cryptic-splice scorer error: no model"],
    ["splice", "pair 70/149: splice (SpliceAI cross-check): splice delta 0.500 > 0.2"],
    ["splice", "pair 70/149: splice cross-check error: model failed"],
    ["docking_a", "pair 70/149: pegRNA A (PBS 13, flap 35): docking: interface ipTM 0.31 is below the neutral extension's 0.58"],
    ["docking_a", "pair 70/149: pegRNA A (PBS 13, flap 35): docking scorer error: service error"],
    ["docking_b", "pair 70/149: pegRNA B (PBS 17, flap 35): docking: interface ipTM 0.31 is below the neutral extension's 0.58"],
    ["docking_b", "pair 70/149: pegRNA B (PBS 17, flap 35): docking scorer error: service error"],
  ];
  for (const [want, line] of kinds)
    check(same(stepsFor(line), [want]), `a "${want}" line belongs to that step and no other`, `${line.slice(0, 70)}... -> ${show(stepsFor(line))}`);

  check(same(stepsFor(summary(5, 6, "")), []) && same(stepsFor("0 pair(s) remain"), []), "a summary line that reports no rejection belongs to no step");
  check(same(stepsFor("pair 1/2: pegRNA AB (PBS 13): docking: x"), []), "'pegRNA A' is not mistaken for another name that starts with it");

  // the real failed run on DSI (CTCF cassette, strict off-target setting, before the nick-continuity fix)
  const real = [
    "24 PAM-in pair(s) from 16 guide(s); 1 failed the sequence rules (GC content 25% is outside 30-75%), 11 failed the off-target screen (each needs a unique on-target and at most 0 off-target site(s) within the gate mismatch tier); 1 pair(s) remain",
    ...[13, 15, 17].map((p) => `pair 101/164: pegRNA A (PBS ${p}, flap 35): cross-region pairing (ViennaRNA): spacer[17..17] C x RTT[34..34] G (1 bp, p >= 0.999)`),
    ...[13, 15, 17].map((p) => `pair 101/164: pegRNA B (PBS ${p}, flap 35): cross-region pairing (ViennaRNA): spacer[17..18] AC x RTT[33..34] GT (2 bp, p >= 0.996)`),
  ];
  const counts = pl.failureCounts(steps, real);
  check(same(counts.map((c) => [c.id, c.count, c.unit]), [["rules", 1, "guide"], ["offtarget", 11, "guide"], ["structure_a", 3, "pegRNA"], ["structure_b", 3, "pegRNA"]]), "the real run is counted: 1 guide by the rules, 11 by the scan, 3 + 3 pegRNAs by folding", show(counts));
  check(steps[pl.findBlockingStep("", "", steps, real)].id === "structure_b", "and the furthest step any candidate reached is the second pegRNA's folding");
  check(pl.findBlockingStep("", "", steps, ["a line no step claims"]) === -1, "a record with no recognisable line blocks no step");
  check(pl.findBlockingStep("", "", steps, []) === -1 && same(pl.failureCounts(steps, []), []), "an empty record counts nothing");

  // a count comes from the line's own number, and a line that has none is one rejection
  const one = pl.failureCounts(steps, [summary(10, 8, "; 7 failed the off-target screen (x)")]);
  check(one.length === 1 && one[0].id === "offtarget" && one[0].count === 7, "the off-target line counts the 7 guides it reports, not 1", show(one));
  const both = pl.failureCounts(steps, [summary(10, 8, "; 2 failed the sequence rules (x), 3 failed the off-target screen (y)")]);
  check(same(both.map((c) => [c.id, c.count]), [["rules", 2], ["offtarget", 3]]), "one summary line can carry two reasons, each with its own number", show(both));
  check(same(pl.failureCounts(steps, [kinds[8][1], kinds[9][1], kinds[10][1]]).map((c) => [c.id, c.count]), [["flaps", 3]]), "template and donor-end rejections are counted together under one step");

  // nothing is dropped: a pair line that no step claims is counted as 'other'
  const unknown = pl.failureCounts(steps, ["pair 5/9: something the page has never heard of", "pair 5/10: and another", kinds[5][1]]);
  check(same(unknown.map((c) => [c.id, c.count]), [["recode", 1], ["other", 2]]), "pair lines no step claims are counted as other reasons and never left out", show(unknown));
  check(!pl.failureCounts(steps, [summary(5, 6, "")]).some((c) => c.id === "other"), "a plain summary line is not counted as an unknown reason");
  check(steps[pl.findBlockingStep("", "", steps, [kinds[13][1], kinds[0][1]])].id === "structure_a", "the furthest step wins over an earlier one, wherever the lines are in the record");
}

// ------------------------------------------------------------------------------------------------------------
console.log("\nTHE STEPS WHILE A RUN STREAMS");
{
  const ev = (cand, gate, passed = true) => ({ cand, gate, info: "", passed, metrics: {} });
  const live = (...events) => pl.paLiveStates(Object.fromEntries(events.map((e) => [e.gate, e])));
  const of = (states) => Object.entries(states).filter(([, s]) => s !== "idle").map(([id, s]) => `${id}:${s}`);

  check(same(of(live()), ["importance:running"]), "before anything reports, the insertion-site check is the running step");
  check(same(of(live(ev(1, "importance"))), ["importance:passed", "pairs:passed", "rules:passed", "offtarget:running"]), "once the site is judged, the guide scan is the running step, not the instant one after it", show(of(live(ev(1, "importance")))));
  const afterFold = live(ev(1, "importance"), ev(1, "structure_a"));
  check(same(of(afterFold), ["importance:passed", "pairs:passed", "rules:passed", "offtarget:passed", "recode:passed", "inv0:passed", "flaps:passed", "strand:passed", "structure_a:passed", "structure_b:running"]), "after the first fold the planning steps are passed and the second pegRNA's fold is running", show(of(afterFold)));
  check(afterFold.cryptic === "idle" && afterFold.docking_a === "idle", "per-pair steps that have not reported are never assumed passed");
  const skipped = live(ev(1, "importance"), ev(1, "structure_a"), ev(1, "structure_b"), ev(1, "docking_a"));
  check(skipped.cryptic === "idle" && skipped.splice === "idle" && skipped.docking_a === "passed" && skipped.docking_b === "running", "a per-pair step that never reported is not claimed as passed because later ones did", show(of(skipped)));
  const failedFold = live(ev(1, "importance"), ev(1, "structure_a", false));
  check(failedFold.structure_a === "blocked" && !Object.values(failedFold).includes("running"), "a step that rejects is blocked, and nothing is running until the next pair starts", show(of(failedFold)));
  check(failedFold.offtarget === "passed" && failedFold.structure_b === "idle", "and the planning behind it is passed, the steps after it idle");
  const dockA = live(ev(1, "importance"), ev(1, "structure_a"), ev(1, "structure_b"), ev(1, "cryptic"), ev(1, "splice"), ev(1, "docking_a"));
  check(dockA.docking_b === "running" && dockA.docking_a === "passed", "while the second pegRNA docks it is the running step");
  const done = live(ev(1, "importance"), ev(1, "structure_a"), ev(1, "structure_b"), ev(1, "cryptic"), ev(1, "splice"), ev(1, "docking_a"), ev(1, "docking_b"));
  check(done.report === "running" && done.docking_b === "passed", "when the last check passes, the report is what is being written");
  const second = pl.paLiveStates(pl.mergeProgress(pl.mergeProgress({ importance: ev(1, "importance"), structure_a: ev(1, "structure_a"), docking_a: ev(1, "docking_a", false) }, ev(2, "structure_a"), "pa"), ev(2, "structure_b"), "pa"));
  check(second.docking_a === "idle" && second.structure_a === "passed" && second.cryptic === "running", "when the second pair starts, the first pair's results are gone and the second's are shown", show(of(second)));
  check(pl.paLatestEvent({}) === null, "with no events there is no newest event");
  check(pl.paLatestEvent({ structure_a: ev(1, "structure_a"), importance: ev(1, "importance") }).gate === "structure_a", "within a pair the furthest step is the newest");
  check(pl.paLatestEvent({ importance: ev(1, "importance"), structure_a: ev(2, "structure_a") }).cand === 2, "a later pair is newer than an earlier one's later step");
  check(live({ ...ev(1, "never_heard_of_it") }).importance === "running", "an event from a step this page does not know does not break it");

  // merging progress
  const p1 = { importance: ev(1, "importance"), structure_a: ev(1, "structure_a"), docking_a: ev(1, "docking_a", false) };
  check(same(Object.keys(pl.mergeProgress(p1, ev(1, "structure_b"), "pa")).sort(), ["docking_a", "importance", "structure_a", "structure_b"]), "a new step of the same pair is added to that pair's results");
  check(same(Object.keys(pl.mergeProgress(p1, ev(2, "structure_a"), "pa")).sort(), ["importance", "structure_a"]) && pl.mergeProgress(p1, ev(2, "structure_a"), "pa").structure_a.cand === 2, "a new pair drops the old pair's results and keeps the insertion-site check");
  check(pl.mergeProgress({ importance: ev(1, "importance"), cryptic: ev(3, "cryptic") }, ev(2, "structure_a"), "pa").cryptic?.cand === 3, "an event that arrives late does not erase a newer pair's results");
  check(same(Object.keys(pl.mergeProgress(p1, ev(2, "structure_a"), "tjpe")).sort(), ["docking_a", "importance", "structure_a"]) && pl.mergeProgress(p1, ev(2, "structure_a"), "tjpe").structure_a.cand === 2, "tj-PE just keeps the newest event of each gate, as it always did");
  check(pl.mergeProgress(p1, ev(2, "structure_a"), "pa") !== p1 && p1.structure_a.cand === 1, "merging never changes the progress it was given");
}

// ------------------------------------------------------------------------------------------------------------
console.log("\nTHE PROGRESS STREAM");
{
  const enc = new TextEncoder();
  const line = (o) => JSON.stringify(o) + "\n";
  const progress = (gate, cand = 1) => ({ cand, gate, info: gate, passed: true, metrics: {} });
  const beat = { type: "heartbeat" };
  const result = { type: "result", ok: true, run_id: "r1" };

  /** A reader that plays [delay in ms, text] steps (text null ends the stream), then goes silent, or fails. */
  function reader(script, failWith) {
    let i = 0;
    const r = { cancelled: 0 };
    r.read = () => {
      if (i >= script.length) return failWith ? Promise.reject(failWith) : new Promise(() => {});
      const [delay, text] = script[i++];
      return new Promise((go) => setTimeout(() => go(text === null ? { done: true } : { done: false, value: enc.encode(text) }), delay));
    };
    r.cancel = () => {
      r.cancelled += 1;
      return Promise.resolve();
    };
    return r;
  }
  /** Run a read and settle it, or report that it hung, so a missing watchdog fails a check and not the whole script. */
  const play = (rd, stallMs) => {
    const seen = [];
    const run = st.readDesignStream(rd, (ev) => seen.push(ev), stallMs).then((value) => ({ value }), (error) => ({ error }));
    const hung = new Promise((go) => setTimeout(() => go({ hung: true }), 1500));
    return Promise.race([run, hung]).then((out) => ({ ...out, seen }));
  };

  const text = line(progress("importance")) + line(beat) + line(progress("structure_a")) + line(result);
  const whole = await play(reader([[0, text.slice(0, 40)], [0, text.slice(40, 150)], [0, text.slice(150)], [0, null]]), 500);
  check(
    whole.value?.run_id === "r1" && same(whole.seen.map((e) => e.gate), ["importance", "structure_a"]),
    "a stream cut in the middle of its lines is read whole: the progress in order, the keepalive not shown as progress, then the result",
    show(whole),
  );
  const mixed = await play(reader([[0, "not json at all\n" + line(progress("importance")) + "{broken\n" + line(result)], [0, null]]), 500);
  check(mixed.value?.run_id === "r1" && mixed.seen.length === 1, "a line that is not JSON is skipped and does not end the stream", show(mixed));

  const silent = reader([[0, line(progress("importance"))]]);
  const t0 = Date.now();
  const stalled = await play(silent, 80);
  const waited = Date.now() - t0;
  check(!stalled.hung && /Lost contact with the design server/.test(stalled.error?.message ?? ""), "a stream that goes silent ends with a message that says so", show({ hung: stalled.hung, error: stalled.error?.message }));
  check(silent.cancelled === 1, "and the connection is closed, once");
  check(waited >= 70 && waited < 1000, "and only after the silence has lasted as long as allowed", `${waited} ms`);
  check(stalled.seen.length === 1, "what arrived before the silence is kept");
  const lateGap = reader([[0, line(progress("importance"))], [300, line(result)]]);
  const gap = await play(lateGap, 80);
  check(/Lost contact/.test(gap.error?.message ?? "") && lateGap.cancelled === 1, "a single gap longer than the limit in the middle of a stream is enough");

  const slow = reader([[40, line(progress("importance"))], [40, line(beat)], [40, line(beat)], [40, line(progress("structure_a"))], [40, line(result)], [40, null]]);
  const alive = await play(slow, 100);
  check(alive.value?.run_id === "r1", "a stream that is slow but never silent for the whole limit is never given up on, however long it runs", show(alive));
  check(Boolean(alive.value) && !alive.hung, "even though the whole read took longer than the limit");

  const noResult = await play(reader([[0, line(progress("importance"))], [0, null]]), 500);
  check(/ended without returning a result/.test(noResult.error?.message ?? ""), "a stream that closes without a result says so", show(noResult.error?.message));
  const aborted = Object.assign(new Error("The user aborted a request."), { name: "AbortError" });
  const cancelled = await play(reader([[0, line(progress("importance"))]], aborted), 500);
  check(cancelled.error === aborted, "a read that fails (an abort, a dropped connection) fails with its own error, unchanged");

  const timers = () => process.getActiveResourcesInfo().filter((x) => x === "Timeout").length;
  const before = timers();
  await st.readDesignStream(reader([[0, line(progress("importance")) + line(result)], [0, null]]), () => undefined);
  check(timers() === before, "a finished stream leaves no timer behind (a 90 s deadline left running would keep a process alive)");

  check(/90 seconds/.test(st.stalledMessage(90_000)) && /Run it again/.test(st.stalledMessage(90_000)), "the message gives the time and says what to do");
  check(st.STREAM_STALL_MS >= 4 * 15_000 && st.STREAM_STALL_MS <= 180_000, "the limit is several times the server's 15 s keepalive, and not unbounded");
  const api = readFileSync(join(SRC, "api.ts"), "utf8");
  check(api.includes("readDesignStream(r.body.getReader(), onProgress)") && !api.includes("reader.read()"), "the page reads its streams with this code and no other");
}

console.log(`\n${failures === 0 ? "all checks passed" : `${failures} check(s) FAILED`}`);
