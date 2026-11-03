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
 * Needs the api-mode token in /tmp/peg-e2e/api-token.txt (or edit TOKEN below).
 */
import { chromium } from "playwright";
import { readFileSync } from "fs";

const SITE = "https://pegrna-designer.vercel.app";
const TOKEN = readFileSync("/tmp/peg-e2e/api-token.txt", "utf8").trim();

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

// Whatever the verdict, no gate may be a placeholder.
const standIn = await page.locator("text=did not run a real model").count();
standIn === 0 ? ok("no stand-in gates reported") : fail("a gate was a placeholder");

const gateRows = await page.locator("ol li").count();
console.log(`  pipeline steps shown: ${gateRows}`);

await page.screenshot({ path: "/tmp/peg-e2e/shots/12-real-design.png", fullPage: true });
console.log("  saved 12-real-design.png");

console.log("\nJS errors:", errors.length ? errors : "none");
if (errors.length) process.exitCode = 1;
await browser.close();
