/**
 * End-to-end browser check of the live public site.
 *
 * Exercises the real thing a user does, against the real backend: load the site, confirm it reached the pipeline
 * cross-origin, run a dry design, then run a REAL api design with live models and assert no gate was a placeholder.
 * There are two real designs, one per way of writing the insert, and each takes minutes, which is the point: it is
 * also the test that the NDJSON stream survives a long silent docking call.
 *
 *   node scripts/qa-browser.mjs [all|tjpe|pa]
 *
 *   all   (default) the public checks, then both real designs one after the other
 *   tjpe  the public checks, a dry run, and the real design with the insert inside the pegRNA (about 15 minutes)
 *   pa    the public checks, and the real Prime Assembly design of the 499 nt CTCF cassette (about 3 minutes)
 *
 * The designer is open to everyone, so this needs no credentials. (The quick checks that need no real run are in
 * scripts/qa-ui-quick.mjs.) Needs Playwright: `playwright` is imported if installed, else `playwright-core`; or set
 * PLAYWRIGHT_PATH to its entry file. PW_CHANNEL=chrome uses the installed Chrome. SHOTS_DIR says where the
 * screenshots go (default /tmp/peg-e2e/shots).
 */
import { mkdirSync } from "node:fs";

async function loadPlaywright() {
  for (const name of [process.env.PLAYWRIGHT_PATH, "playwright", "playwright-core"].filter(Boolean)) {
    try {
      return await import(name);
    } catch {
      /* try the next one */
    }
  }
  throw new Error("Playwright not found: npm i playwright, or set PLAYWRIGHT_PATH to its entry file.");
}
const { chromium } = await loadPlaywright();

const SITE = (process.env.SITE ?? "https://pegrna-designer.vercel.app").replace(/\/$/, "");
const WHICH = (process.argv[2] ?? "all").toLowerCase();
if (!["all", "tjpe", "pa"].includes(WHICH)) {
  console.error("usage: node scripts/qa-browser.mjs [all|tjpe|pa]");
  process.exit(2);
}
const SHOTS = process.env.SHOTS_DIR ?? "/tmp/peg-e2e/shots";
mkdirSync(SHOTS, { recursive: true });

const browser = await chromium.launch(process.env.PW_CHANNEL ? { channel: process.env.PW_CHANNEL } : {});
const ctx = await browser.newContext({ viewport: { width: 1440, height: 1050 } });
const page = await ctx.newPage();
const errors = [];
page.on("pageerror", (e) => errors.push(`pageerror: ${e.message}`));
page.on("console", (m) => m.type() === "error" && errors.push(`console: ${m.text()}`));

const fail = (m) => {
  console.log(`  FAIL  ${m}`);
  process.exitCode = 1;
};
const ok = (m) => console.log(`  PASS  ${m}`);

console.log(`PUBLIC SITE  (${SITE}, ${WHICH})`);
await page.goto(SITE, { waitUntil: "networkidle", timeout: 60000 });
(await page.locator("h1").first().innerText()) === "pegRNA designer"
  ? ok("app renders")
  : fail("app did not render");

await page.waitForTimeout(2500);
const chips = (await page.locator("header [aria-label='Service status']").allInnerTexts()).join(" ");
chips.includes("docking live") ? ok("backend reachable, docking live") : fail(`chips: ${chips}`);

const banner = await page.locator("text=The design server is not responding").isVisible().catch(() => false);
banner ? fail("backend-down banner is showing") : ok("no backend-down banner");

if (WHICH !== "pa") {
  console.log("\nDRY RUN (open to everyone)");
  await page.getByRole("button", { name: "Load CTCF example" }).click();
  await page.waitForFunction(() => document.querySelector("#genome")?.value?.length > 100, { timeout: 30000 });
  await page.getByRole("button", { name: "Design pegRNA", exact: true }).click();
  await page.waitForSelector("text=/Design found|No design/", { timeout: 180000 });
  (await page.locator("text=Design found").isVisible()) ? ok("fake design completes") : fail("fake design failed");
  (await page.locator("text=Dry run — this is not a validated design").isVisible())
    ? ok("dry-run clearly disclosed")
    : fail("dry-run disclosure missing");
  const segs = await page.locator("[aria-label='pegRNA parts, to scale'] button").count();
  segs === 7 ? ok(`construct map renders (${segs} parts)`) : fail(`map had ${segs} parts`);

  console.log("\nREAL MODELS, INSERT INSIDE THE PEGRNA (open to everyone)");
  await page.selectOption("#mode", "api");
  await page.waitForTimeout(600);
  (await page.locator("#access-token").count()) === 0
    ? ok("no token asked for in api mode")
    : fail("a token field is shown in api mode");

  await page.getByRole("button", { name: "Load example for real models" }).click();
  await page.waitForTimeout(1500);
  await page.selectOption("#mode", "api");
  await page.waitForTimeout(400);

  // The second nick's off-target tolerance is a visible, validated setting and defaults to strict.
  await page.locator("summary", { hasText: "Fine-tuning" }).click();
  const tolerance = page.locator("#ngrna_max_offtargets");
  (await tolerance.inputValue()) === "0"
    ? ok("off-target tolerance defaults to 0 (strict)")
    : fail(`off-target tolerance defaulted to ${await tolerance.inputValue()}`);
  await tolerance.fill("11");
  await tolerance.blur();
  (await page.locator("text=A whole number from 0 to 10.").isVisible())
    ? ok("an out-of-range tolerance is refused")
    : fail("no validation message for a tolerance of 11");
  await tolerance.fill("0");

  console.log("  running a real design (PRIDICT2 + cas-offinder + SpliceAI + AlphaFold + Boltz-2)...");
  const t0 = Date.now();
  await page.getByRole("button", { name: "Design pegRNA", exact: true }).click();

  // A real run takes ~15 min; the heartbeat must hold the connection open.
  await page.waitForSelector("text=/Design found|No design/", { timeout: 2_100_000 });
  const mins = ((Date.now() - t0) / 60000).toFixed(1);
  ok(`real design completed in ${mins} min without the stream dropping`);

  const verdict = (await page.locator("text=/Design found|No design/").first().innerText()).trim();
  console.log(`  verdict: ${verdict}`);

  // The shipped real-models example is the demo. It must produce a design: it failed once with
  // "no candidate survived after 0 tried" because of a pool-construction bug, and a browser
  // check that accepts either verdict would have passed straight through that.
  verdict === "Design found"
    ? ok("the shipped demo produces a design")
    : fail(`the shipped demo returned "${verdict}"`);

  // Whatever the verdict, no gate may be a placeholder.
  const standIn = await page.locator("text=did not run a real model").count();
  standIn === 0 ? ok("no stand-in gates reported") : fail("a gate was a placeholder");

  // The docking step must report its length-matched comparison, not the retired 0.6 floor.
  const docking = await page.locator("li", { hasText: "Cas9 binding" }).first().innerText();
  /neutral/i.test(docking) && /non-binder/i.test(docking)
    ? ok("docking is judged against a length-matched reference")
    : fail(`docking row does not show the calibrated comparison: ${docking.slice(0, 160)}`);

  // The off-target row states the chosen guide's actual count, not a fixed "0".
  const offtarget = await page.locator("li", { hasText: "Off-target check" }).first().innerText().catch(() => "");
  /off-target site/i.test(offtarget)
    ? ok("the off-target row reports the chosen guide's count")
    : fail(`off-target row does not state a count: ${offtarget.slice(0, 160)}`);

  // Advisories must be readable without opening the collapsed report. The chosen guide of the
  // shipped demo has a CCCCC run: that is reported, not hidden and not a reason to reject.
  if ((await page.locator("text=Check these before you order").count()) > 0) {
    ok("advisories are shown on the result");
    const advice = await page.locator("ul li", { hasText: "advisory" }).allInnerTexts();
    console.log(`  ${advice.length} advisory line(s): ${advice.map((a) => a.slice(0, 70)).join(" | ")}`);
  } else {
    console.log("  (no advisories on this design)");
  }

  const gateRows = await page.locator("ol li").count();
  console.log(`  pipeline steps shown: ${gateRows}`);

  await page.screenshot({ path: `${SHOTS}/12-real-design.png`, fullPage: true });
  console.log("  saved 12-real-design.png");
}

if (WHICH !== "tjpe") {
  console.log("\nREAL MODELS, PRIME ASSEMBLY (the 499 nt CTCF cassette as a separate donor; open to everyone)");
  // A clean page, so nothing of an earlier design is on it.
  await page.goto(SITE, { waitUntil: "networkidle", timeout: 60000 });
  await page.getByRole("button", { name: /Load cassette example/ }).click();
  await page.waitForFunction(() => document.querySelector("#flap_length") !== null && document.querySelector("#genome")?.value?.length > 100, { timeout: 30000 });
  await page.selectOption("#mode", "api");
  await page.waitForTimeout(500);
  (await page.getByRole("radio", { name: /separate donor DNA/ }).getAttribute("aria-checked")) === "true"
    ? ok("the cassette example selects Prime Assembly")
    : fail("the cassette example did not select Prime Assembly");

  const steps = [
    "Protein tolerance", "Pair the guides", "Guide sequence rules", "Off-target check", "Build the donor arms", "No deletion (INV-0)",
    "Write the flaps", "Strand check", "pegRNA A folding", "pegRNA B folding", "Cryptic splice sites", "Splicing impact",
    "pegRNA A binds Cas9", "pegRNA B binds Cas9", "Report",
  ];
  const stateOf = async (title) => {
    const row = page.locator("ol > li", { hasText: title }).first();
    return (await row.count()) ? (await row.locator("span.ml-auto").first().innerText().catch(() => "")).trim() || "idle" : "missing";
  };

  console.log("  running a real design (cas-offinder + AlphaFold + ViennaRNA + SpliceAI + Boltz-2)...");
  const t1 = Date.now();
  await page.getByRole("button", { name: /Design pegRNAs and donor/ }).click();

  // While it runs, watch the panel as a person would: which steps are shown as running. The docking calls are slow
  // (tens of seconds each), so a step that is really running must be caught by looking every second or two.
  const running = new Set();
  const verdictShown = () => page.locator("text=/Design found|No design/").first().isVisible().catch(() => false);
  while (!(await verdictShown())) {
    if (Date.now() - t1 > 1_200_000) {
      fail("the Prime Assembly design did not finish in 20 min");
      break;
    }
    for (const title of steps) if ((await stateOf(title)) === "running") running.add(title);
    await page.waitForTimeout(1200);
  }
  const secs = ((Date.now() - t1) / 1000).toFixed(0);
  ok(`real design completed in ${secs} s without the stream dropping`);
  console.log(`  steps seen running: ${[...running].join(", ") || "none"}`);
  [...running].some((t) => /binds Cas9/.test(t))
    ? ok("the panel showed a Cas9-binding step as running while the models ran")
    : fail(`the panel never showed a docking step running (saw: ${[...running].join(", ") || "nothing"})`);

  const verdict = (await page.locator("text=/Design found|No design/").first().innerText()).trim();
  console.log(`  verdict: ${verdict}`);
  verdict === "Design found" ? ok("the cassette example produces a design") : fail(`the cassette example returned "${verdict}"`);

  // The same guarantees as any real design: nothing a placeholder, and the verdicts are real ones.
  (await page.locator("text=did not run a real model").count()) === 0 && (await page.locator("text=Dry run — this is not a validated design").count()) === 0
    ? ok("no stand-in gates and no dry-run banner")
    : fail("a gate was a placeholder, or the result is flagged as a dry run");
  const states = await Promise.all(steps.map(stateOf));
  states.every((s) => s === "passed") ? ok("every one of the fifteen steps is shown as passed") : fail(`step states: ${steps.map((t, i) => `${t}=${states[i]}`).join(", ")}`);
  for (const side of ["A", "B"]) {
    const row = await page.locator("ol > li", { hasText: `pegRNA ${side} binds Cas9` }).first().innerText();
    /interface ipTM/i.test(row) && /within cutoff/i.test(row) ? ok(`docking of pegRNA ${side} is a measured verdict against a cutoff`) : fail(`docking row ${side}: ${row.slice(0, 160)}`);
  }

  // The order sheet.
  const headings = (await page.locator("h3, h4").allInnerTexts()).map((t) => t.trim());
  ["What to order", "pegRNA A", "pegRNA B", "Donor DNA", "Where it goes"].every((h) => headings.includes(h))
    ? ok("the result is an order sheet: both pegRNAs, the donor, and where it goes")
    : fail(`the result headings are ${JSON.stringify(headings)}`);
  (await page.getByText("This combination is untested").count()) === 1 ? ok("the untested-combination caveat is shown") : fail("the caveat callout is missing");
  const lengths = await page.locator("h4", { hasText: /^pegRNA [AB]$/ }).evaluateAll((hs) => hs.map((h) => Number(/(\d+) nt/.exec(h.parentElement.innerText)?.[1])));
  lengths.length === 2 && lengths.every((n) => n > 0 && n <= 166) ? ok(`both pegRNAs are within 166 nt (${lengths.join(" and ")})`) : fail(`pegRNA lengths: ${JSON.stringify(lengths)}`);

  // What a person would download and order from.
  const link = page.locator("a[href*='recoded_pa_donor.fasta']").first();
  const href = await link.getAttribute("href").catch(() => null);
  if (href) {
    const fasta = await (await page.request.get(href)).text();
    const records = Object.fromEntries(fasta.split(">").slice(1).map((r) => { const [head, ...seq] = r.trim().split("\n"); return [head.split(/\s/)[0], seq.join("")]; }));
    const shownDonor = Number(/Donor DNA\s+(\d+) nt/.exec(await page.locator("h4", { hasText: "Donor DNA" }).locator("xpath=..").innerText())?.[1]);
    Object.keys(records).join(",") === "pegRNA_A,pegRNA_B,donor" && records.donor.length === shownDonor
      ? ok(`the FASTA has the three sequences and its donor is the ${shownDonor} nt one on the page`)
      : fail(`the FASTA holds ${Object.keys(records).join(",")}, donor ${records.donor?.length} nt, the page says ${shownDonor}`);
  } else {
    fail("there is no FASTA download on the result");
  }
  (await page.getByText("How to prepare the donor", { exact: true }).count()) === 1 ? ok("how to prepare the donor is offered") : fail("no notes on preparing the donor");

  await page.screenshot({ path: `${SHOTS}/13-real-prime-assembly.png`, fullPage: true });
  console.log("  saved 13-real-prime-assembly.png");
}

console.log("\nJS errors:", errors.length ? errors : "none");
if (errors.length) process.exitCode = 1;
await browser.close();
