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
 * Needs the api-mode token: macOS keychain item `pegrna-api-mode-token`, or the
 * PEGRNA_API_MODE_TOKEN environment variable.
 */
import { chromium } from "playwright";
import { execSync } from "child_process";

const SITE = "https://pegrna-designer.vercel.app";
// The api-mode token lives in the macOS keychain (service: pegrna-api-mode-token), never in
// a file in the repo. PEGRNA_API_MODE_TOKEN overrides it, for CI or non-macOS hosts.
function loadToken() {
  if (process.env.PEGRNA_API_MODE_TOKEN) return process.env.PEGRNA_API_MODE_TOKEN;
  try {
    return execSync("security find-generic-password -s pegrna-api-mode-token -w", {
      stdio: ["ignore", "pipe", "ignore"],
    })
      .toString()
      .trim();
  } catch {
    throw new Error(
      "No api-mode token: set PEGRNA_API_MODE_TOKEN, or store one with\n" +
        "  security add-generic-password -a $USER -s pegrna-api-mode-token -w <token>",
    );
  }
}
const TOKEN = loadToken();

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

console.log("\nREAL MODELS (token-gated)");
await page.selectOption("#mode", "api");
await page.waitForTimeout(600);
(await page.locator("#access-token").isVisible()) ? ok("token required in api mode") : fail("no token field");
await page.fill("#access-token", TOKEN);

await page.getByRole("button", { name: "Load example for real models" }).click();
await page.waitForTimeout(1500);
await page.selectOption("#mode", "api");
await page.waitForTimeout(400);

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

const gateRows = await page.locator("ol li").count();
console.log(`  pipeline steps shown: ${gateRows}`);

await page.screenshot({ path: "/tmp/peg-e2e/shots/12-real-design.png", fullPage: true });
console.log("  saved 12-real-design.png");

console.log("\nJS errors:", errors.length ? errors : "none");
if (errors.length) process.exitCode = 1;
await browser.close();
