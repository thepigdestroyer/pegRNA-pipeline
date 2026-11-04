/**
 * End-to-end browser check of the live public site.
 *
 * Exercises the real thing a user does, against the real backend: load the
 * site, confirm it reached the pipeline cross-origin, run a dry design, then
 * run a REAL api design with live models and assert no gate was a placeholder.
 *
 * The real design takes minutes, which is the point -- it is also the test that
 * the NDJSON stream survives a long silent docking call.
 *
 *   node scripts/qa-browser.mjs
 *
 * The designer is open to everyone, so this needs no credentials. (The quick checks that
 * need no real run are in scripts/qa-ui-quick.mjs.)
 */
import { chromium } from "playwright";

const SITE = "https://pegrna-designer.vercel.app";

const browser = await chromium.launch();
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

console.log("PUBLIC SITE");
await page.goto(SITE, { waitUntil: "networkidle", timeout: 60000 });
(await page.locator("h1").first().innerText()) === "pegRNA designer"
  ? ok("app renders")
  : fail("app did not render");

await page.waitForTimeout(2500);
const chips = (await page.locator("header [aria-label='Service status']").allInnerTexts()).join(" ");
chips.includes("docking live") ? ok("backend reachable, docking live") : fail(`chips: ${chips}`);

const banner = await page.locator("text=The design server is not responding").isVisible().catch(() => false);
banner ? fail("backend-down banner is showing") : ok("no backend-down banner");

console.log("\nDRY RUN (open to everyone)");
await page.getByRole("button", { name: "Load CTCF example" }).click();
await page.waitForFunction(() => document.querySelector("#genome")?.value?.length > 100, { timeout: 30000 });
await page.getByRole("button", { name: "Design pegRNA" }).click();
await page.waitForSelector("text=/Design found|No design/", { timeout: 180000 });
(await page.locator("text=Design found").isVisible()) ? ok("fake design completes") : fail("fake design failed");
(await page.locator("text=Dry run — this is not a validated design").isVisible())
  ? ok("dry-run clearly disclosed")
  : fail("dry-run disclosure missing");
const segs = await page.locator("[aria-label='pegRNA parts, to scale'] button").count();
segs === 7 ? ok(`construct map renders (${segs} parts)`) : fail(`map had ${segs} parts`);

console.log("\nREAL MODELS (open to everyone)");
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
await page.getByRole("button", { name: "Design pegRNA" }).click();

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

await page.screenshot({ path: "/tmp/peg-e2e/shots/12-real-design.png", fullPage: true });
console.log("  saved 12-real-design.png");

console.log("\nJS errors:", errors.length ? errors : "none");
if (errors.length) process.exitCode = 1;
await browser.close();
