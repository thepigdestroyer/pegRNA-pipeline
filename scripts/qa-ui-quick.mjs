/**
 * Quick UI checks: no token, no real-model run, a few seconds (the dry-run designs take about a second each).
 *
 *   node scripts/qa-ui-quick.mjs [url]      (default: the live site)
 *
 * THEME. The toggle decides the theme, not the operating system, and every part of the page
 * follows it. The bug this guards against: Tailwind's `dark:` utilities followed the OS
 * setting while the body and cards followed the toggle, so whenever the two disagreed the
 * header, footer and action bar stayed in the wrong theme. Each combination of OS setting and
 * chosen theme is loaded and the colours that make up the page are compared.
 *
 * Needs Playwright. `playwright` is imported if installed, else `playwright-core`; or set
 * PLAYWRIGHT_PATH to its entry file. PW_CHANNEL=chrome uses the installed Chrome.
 */
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

const SITE = (process.argv[2] ?? "https://pegrna-designer.vercel.app").replace(/\/$/, "") + "/";
const { chromium } = await loadPlaywright();
const browser = await chromium.launch(process.env.PW_CHANNEL ? { channel: process.env.PW_CHANNEL } : {});

const errors = [];
const fail = (m) => {
  console.log(`  FAIL  ${m}`);
  process.exitCode = 1;
};
const ok = (m) => console.log(`  PASS  ${m}`);

/**
 * Open the app with an OS colour scheme, optionally a previously saved theme choice, and
 * optionally other localStorage entries a returning visitor's browser might hold.
 */
async function openApp({ os, saved, storage = {} }) {
  const ctx = await browser.newContext({ colorScheme: os, viewport: { width: 1280, height: 900 } });
  if (saved) await ctx.addInitScript((t) => localStorage.setItem("pegrna.theme", t), saved);
  await ctx.addInitScript((kv) => {
    for (const [k, v] of Object.entries(kv)) localStorage.setItem(k, v);
  }, storage);
  const page = await ctx.newPage();
  page.on("pageerror", (e) => errors.push(`pageerror: ${e.message}`));
  page.on("console", (m) => m.type() === "error" && errors.push(`console: ${m.text()}`));
  await page.goto(SITE, { waitUntil: "domcontentloaded", timeout: 60000 });
  await page.waitForSelector("header button[aria-label^='Switch to']", { timeout: 30000 });
  return { ctx, page };
}

/** The colours that make up the page, read as the browser actually resolved them. */
const readLook = (page) =>
  page.evaluate(() => {
    const css = (el, prop) => (el ? getComputedStyle(el).getPropertyValue(prop) : null);
    const header = document.querySelector("header");
    const metas = [...document.querySelectorAll('meta[name="theme-color"]')].filter(
      (m) => !m.media || matchMedia(m.media).matches,
    );
    return {
      "dark class": String(document.documentElement.classList.contains("dark")),
      "native controls (color-scheme)": css(document.documentElement, "color-scheme"),
      "page background": css(document.body, "background-color"),
      "header background": css(header, "background-color"),
      "header title": css(header?.querySelector("h1"), "color"),
      "header subtitle": css(header?.querySelector("p"), "color"),
      "action bar background": css(document.querySelector("div.sticky.bottom-0"), "background-color"),
      "footer text": css(document.querySelector("footer"), "color"),
      "card background": css(document.querySelector(".card"), "background-color"),
      "browser theme colour": metas[0]?.content ?? null,
    };
  });

const differing = (a, b) => Object.keys(a).filter((k) => a[k] !== b[k]);

console.log(`THEME  (${SITE})`);
const THEMES = ["light", "dark"];
const look = {};
for (const os of THEMES) {
  for (const theme of THEMES) {
    const { ctx, page } = await openApp({ os, saved: theme });
    look[`${os}/${theme}`] = await readLook(page);
    await ctx.close();
  }
}

for (const theme of THEMES) {
  const off = differing(look[`light/${theme}`], look[`dark/${theme}`]);
  off.length === 0
    ? ok(`saved "${theme}" looks the same whether the OS is light or dark`)
    : fail(`saved "${theme}" changes with the OS setting: ${off.join(", ")}`);
}

// The look in each theme must be different in every part, or that part is not switching.
const stuck = Object.keys(look["light/light"]).filter(
  (k) => look["light/light"][k] === look["light/dark"][k],
);
stuck.length === 0
  ? ok("every part of the page changes between light and dark")
  : fail(`these parts look identical in light and dark: ${stuck.join(", ")}`);

// The toggle itself, from a first visit (no saved theme), under each OS setting.
for (const os of THEMES) {
  const other = os === "light" ? "dark" : "light";
  const { ctx, page } = await openApp({ os });
  const toggle = page.locator("header button[aria-label^='Switch to']");
  const settle = (want) =>
    page.waitForFunction((w) => document.documentElement.classList.contains("dark") === (w === "dark"), want);

  const first = await readLook(page);
  const offStart = differing(first, look[`${os}/${os}`]);
  offStart.length === 0
    ? ok(`OS ${os}: a first visit starts in ${os}`)
    : fail(`OS ${os}: a first visit does not match saved "${os}": ${offStart.join(", ")}`);

  await toggle.click();
  await settle(other);
  await page.waitForTimeout(250);
  const flipped = await readLook(page);
  // Compared with the look of the target theme when the OS agrees with it, so a part that still
  // follows the OS setting cannot hide by matching a baseline taken under the same setting.
  const offFlip = differing(flipped, look[`${other}/${other}`]);
  offFlip.length === 0
    ? ok(`OS ${os}: the toggle switches the whole page to ${other}`)
    : fail(`OS ${os}: after the toggle these parts are not in ${other}: ${offFlip.join(", ")}`);

  await toggle.click();
  await settle(os);
  await page.waitForTimeout(250);
  const back = await readLook(page);
  const offBack = differing(back, first);
  offBack.length === 0
    ? ok(`OS ${os}: toggling back restores ${os}`)
    : fail(`OS ${os}: toggling back leaves these parts behind: ${offBack.join(", ")}`);
  await ctx.close();
}

// ------------------------------------------------------------------------------------------
// OPEN ACCESS. Nothing asks for a token. A returning visitor whose browser still holds the
// token the old UI stored must have it cleared, and it must never be sent.
console.log(`\nOPEN ACCESS  (${SITE})`);
{
  const STALE = "pegrna.accessToken";
  const { ctx, page } = await openApp({ os: "light", storage: { [STALE]: "stale-token-from-the-old-ui" } });

  // The modes list arrives from the server after the first paint, so wait for it.
  const offered = await page
    .waitForSelector('#mode option[value="api"]', { state: "attached", timeout: 15000 })
    .then(() => true, () => false);
  offered ? ok("real models are offered") : fail("real models are not offered");
  await page.selectOption("#mode", "api");
  await page.waitForTimeout(600);

  (await page.locator("#access-token, input[type=password]").count()) === 0
    ? ok("no token field in real-models mode")
    : fail("a token field is shown in real-models mode");
  const words = (await page.locator("body").innerText()).match(/access token|needs a token|required for real models/i);
  words ? fail(`the page still talks about a token: "${words[0]}"`) : ok("the page never mentions an access token");
  (await page.locator("text=Runs the real models").count()) > 0
    ? ok("real-models mode explains what it does, with no condition")
    : fail("real-models mode has no plain explanation");
  (await page.evaluate((k) => localStorage.getItem(k), STALE)) === null
    ? ok("a token left in the browser by the old UI is cleared")
    : fail("a stale token is still stored in the browser");

  // Run a design in real-models mode with the request answered by a stub, so nothing real runs.
  let sent = null;
  const answered = new Promise((resolve) => {
    page.route("**/api/design/stream", async (route) => {
      sent = route.request();
      resolve();
      await route.fulfill({
        status: 200,
        contentType: "application/x-ndjson",
        headers: { "access-control-allow-origin": "*", "access-control-allow-headers": "*" },
        body:
          JSON.stringify({ type: "result", ok: false, run_id: "ui-check", mode: "api", returncode: 1,
            error: "answered by the UI check", stderr: "", report: "", components: [], has_csv: false,
            docking_stub: false, gates: {}, het_edit: null }) + "\n",
      });
    });
  });
  await page.getByRole("button", { name: "Load example for real models" }).click();
  await page.waitForFunction(() => document.querySelector("#genome")?.value?.length > 100, { timeout: 30000 });
  await page.selectOption("#mode", "api");
  await page.getByRole("button", { name: "Design pegRNA" }).click();
  await Promise.race([answered, page.waitForTimeout(15000)]);
  if (!sent) {
    fail("the real-models design request was never sent");
  } else {
    sent.postDataJSON().mode === "api" ? ok("the request is a real-models design") : fail("the request was not api mode");
    Object.keys(sent.headers()).some((h) => h.toLowerCase() === "x-access-token")
      ? fail("the design request carries an access token")
      : ok("the design request carries no token");
  }
  await ctx.close();
}

// ------------------------------------------------------------------------------------------
// AUTHOR LINK. The footer the visitor sees credits the author with a real link to the
// portfolio. The same link used to exist only in the static fallback that React replaces on
// load, so the rendered page, which is what Google indexes, had no link back to the author.
console.log(`\nAUTHOR LINK  (${SITE})`);
{
  const { ctx, page } = await openApp({ os: "light" });
  const link = page.locator("footer a[rel~='author']");
  const n = await link.count();
  n === 1 ? ok("the footer has one author link") : fail(`the footer has ${n} author links`);
  if (n === 1) {
    const href = await link.getAttribute("href");
    href === "https://mrinal.is-a.dev/" ? ok("it goes to the portfolio") : fail(`it goes to ${href}`);
    const text = (await link.innerText()).trim();
    text === "Mrinal Joshi" ? ok("it is named Mrinal Joshi") : fail(`it is named "${text}"`);
    (await link.isVisible()) ? ok("it is visible") : fail("it is hidden");
    const underline = await link.evaluate((el) => getComputedStyle(el).textDecorationLine);
    underline.includes("underline") ? ok("it is underlined, so it reads as a link") : fail("it is not marked as a link");
  }
  await ctx.close();
}

// ------------------------------------------------------------------------------------------
// HEADER ON A PHONE. The status chips arrive after the backend answers. They used to share one
// row with the title, which crushed the subtitle into a column one word wide (a header 400 px
// tall) and, when the chips appeared, pushed the whole form down by about 300 px (a layout shift
// of 0.39, "poor" in Core Web Vitals). The backend is slowed down here so the chips arrive late,
// as they do for a real visitor on a mobile connection.
console.log(`\nHEADER ON A PHONE  (${SITE})`);
{
  const ctx = await browser.newContext({ viewport: { width: 412, height: 823 }, isMobile: true, hasTouch: true, deviceScaleFactor: 2 });
  await ctx.addInitScript(() => {
    window.__cls = 0;
    new PerformanceObserver((list) => {
      for (const e of list.getEntries()) if (!e.hadRecentInput) window.__cls += e.value;
    }).observe({ type: "layout-shift", buffered: true });
  });
  const page = await ctx.newPage();
  await page.route("**/api/health", async (route) => {
    await new Promise((r) => setTimeout(r, 1500));
    await route.continue();
  });
  await page.goto(SITE, { waitUntil: "domcontentloaded", timeout: 60000 });
  const banner = page.getByRole("banner").first();
  await banner.waitFor({ timeout: 30000 });
  await page.waitForTimeout(400); // the app is drawn, the chips have not arrived
  const before = (await banner.boundingBox()).height;
  const chip = page.getByText(/service ready|pipeline not found/).first();
  await chip.waitFor({ state: "visible", timeout: 30000 });
  await page.waitForTimeout(800);
  const after = (await banner.boundingBox()).height;
  const sub = (await banner.locator("p").first().boundingBox()).width;
  const cls = await page.evaluate(() => window.__cls);

  Math.abs(after - before) <= 2
    ? ok(`the header does not change height when the status chips arrive (${Math.round(before)} px)`)
    : fail(`the header grows from ${Math.round(before)} px to ${Math.round(after)} px when the chips arrive`);
  after <= 200 ? ok(`the header is compact on a phone (${Math.round(after)} px)`) : fail(`the header is ${Math.round(after)} px tall on a phone`);
  sub >= 250 ? ok(`the subtitle has room to read (${Math.round(sub)} px wide)`) : fail(`the subtitle is squeezed to ${Math.round(sub)} px wide`);
  cls < 0.02 ? ok(`no layout shift while the page loads (${cls.toFixed(3)})`) : fail(`layout shift ${cls.toFixed(3)} while the page loads`);

  // The chips and the theme button are still there, and still reachable, on a phone.
  (await chip.isVisible()) ? ok("the status chips are shown") : fail("the status chips are missing");
  (await banner.locator("button[aria-label^='Switch to']").isVisible()) ? ok("the theme button is shown") : fail("the theme button is missing");
  await ctx.close();
}


// ------------------------------------------------------------------------------------------
// PRIME ASSEMBLY. The insert can be written inside the pegRNA (tj-PE, the default) or supplied as a separate
// donor DNA (Prime Assembly). The page must offer both, show each its own settings and nothing of the other's,
// advise (never force) which fits the insert, refuse a request that cannot work before sending it, send the
// strategy and only that strategy's settings, show the order sheet for a Prime Assembly design, and leave
// tj-PE exactly as it was. Dry runs only: a real-model run is scripts/qa-browser.mjs.
console.log(`\nPRIME ASSEMBLY  (${SITE})`);
{
  const { ctx, page } = await openApp({ os: "light" });
  const sent = [];
  page.on("request", (r) => r.url().includes("/api/design/stream") && sent.push(r.postDataJSON()));
  const radio = (re) => page.getByRole("radio", { name: re });
  const has = async (sel) => (await page.locator(sel).count()) > 0;
  const designButton = (name) => page.getByRole("button", { name });
  await radio(/Inside the pegRNA/).waitFor({ timeout: 30000 });

  // The two ways are offered, and tj-PE is the default with its own settings.
  (await page.getByRole("radio").count()) === 2 ? ok("two ways to write the insert are offered") : fail("the chooser does not offer exactly two ways");
  (await radio(/Inside the pegRNA/).getAttribute("aria-checked")) === "true"
    ? ok("inside the pegRNA is the default")
    : fail("inside the pegRNA is not the default");
  (await has("#pegrna_strand")) && !(await has("#flap_length"))
    ? ok("the default shows tj-PE's settings and none of Prime Assembly's")
    : fail("the default shows the wrong settings");

  // The cassette example is a Prime Assembly config: it selects Prime Assembly and fills its settings.
  await page.getByRole("button", { name: /Load cassette example/ }).click();
  await page.waitForFunction(() => document.querySelector("#genome")?.value?.length > 100, { timeout: 30000 });
  (await radio(/separate donor DNA/).getAttribute("aria-checked")) === "true"
    ? ok("the cassette example selects Prime Assembly")
    : fail("the cassette example did not select Prime Assembly");
  const filled = await page.evaluate(() => ({
    flap: document.querySelector("#flap_length")?.value,
    off: document.querySelector("#max_offtargets")?.value,
    cargo: document.querySelector("#cargo")?.value.length,
  }));
  filled.flap === "35" && filled.off === "0" && filled.cargo === 499
    ? ok("it fills the flap length (35), the off-target allowance (0) and the 499 nt cassette")
    : fail(`the example filled ${JSON.stringify(filled)}`);
  // And the other way round: a tj-PE example takes the page back to tj-PE and to its settings.
  await page.getByRole("button", { name: "Load CTCF example" }).click();
  await page.waitForFunction(() => document.querySelector("#pegrna_strand") !== null, { timeout: 5000 }).catch(() => undefined);
  (await radio(/Inside the pegRNA/).getAttribute("aria-checked")) === "true" && !(await has("#flap_length"))
    ? ok("a tj-PE example takes the page from Prime Assembly back to tj-PE and its settings")
    : fail("a tj-PE example left the page on Prime Assembly");
  await page.getByRole("button", { name: /Load cassette example/ }).click();
  await page.waitForFunction(() => document.querySelector("#flap_length") !== null, { timeout: 5000 }).catch(() => undefined);
  const tjOnly = (await has("#pegrna_strand")) || (await has("#arm_side")) || (await has("#spacer")) || (await has("#het_snps")) || (await has("#rtt_aa"));
  tjOnly ? fail("Prime Assembly still shows tj-PE's settings") : ok("Prime Assembly shows none of tj-PE's settings (strand, arm side, own guide, allele-specific, anchor)");
  (await page.locator("text=Fine-tuning").count()) > 0 ? ok("the less common settings are under Fine-tuning") : fail("there is no Fine-tuning panel");

  // The pipeline panel describes the Prime Assembly checks, not tj-PE's.
  const panel = await page.locator("ol").first().innerText();
  const order = ["Protein tolerance", "Pair the guides", "Guide sequence rules", "Off-target check", "Build the donor arms", "No deletion (INV-0)", "Write the flaps", "Strand check", "pegRNA A folding", "pegRNA B folding", "Cryptic splice sites", "Splicing impact", "pegRNA A binds Cas9", "pegRNA B binds Cas9", "Report"];
  const at = order.map((t) => panel.indexOf(t));
  at.every((i) => i >= 0) && at.every((i, k) => k === 0 || i > at[k - 1])
    ? ok("the pipeline panel lists Prime Assembly's fifteen steps in the order the walk does them")
    : fail(`the pipeline panel's steps are missing or out of order: ${JSON.stringify(at)}`);

  // Advice: it recommends, with a button, and never switches by itself.
  await radio(/Inside the pegRNA/).click();
  const advice = await page.getByText(/very long pegRNA/).count();
  advice === 1 ? ok("a 499 nt insert written inside the pegRNA is advised against, with the reason") : fail(`no advice for a long insert inside the pegRNA (${advice})`);
  await page.getByRole("button", { name: "Use Prime Assembly" }).click();
  (await radio(/separate donor DNA/).getAttribute("aria-checked")) === "true"
    ? ok("the advice button switches to Prime Assembly")
    : fail("the advice button did nothing");

  // A short insert cannot be assembled: said beside the choice, said at the cargo, refused before any request.
  const cassette = await page.locator("#cargo").inputValue();
  await page.locator("#cargo").fill("GACTACAAGGACGACGATGACAAG");
  (await page.getByText(/needs at least 70 nt/).count()) > 0 ? ok("a 24 nt insert is advised against for Prime Assembly (70 nt minimum)") : fail("no advice for a too-short insert");
  await designButton(/Design pegRNAs and donor/).click();
  await page.waitForTimeout(400);
  const refusal = await page.locator("[role=alert]").allInnerTexts();
  refusal.some((t) => /Fix the highlighted fields/.test(t)) && refusal.some((t) => /at least 70 nt/.test(t))
    ? ok("the request is refused on the page, naming the minimum")
    : fail(`a too-short insert was not refused: ${JSON.stringify(refusal)}`);
  sent.length === 0 ? ok("and nothing was sent to the server") : fail("a request that cannot work was sent");
  await page.getByRole("button", { name: "Write it inside the pegRNA" }).click();
  (await radio(/Inside the pegRNA/).getAttribute("aria-checked")) === "true" ? ok("the other advice button switches back") : fail("the switch-back button did nothing");
  await radio(/separate donor DNA/).click();
  await page.locator("#cargo").fill(cassette);

  // Settings are checked as they are typed.
  await page.locator("#flap_length").fill("20");
  await page.locator("#flap_length").blur();
  (await page.getByText("A whole number from 30 to 50.").count()) === 1 ? ok("a flap length outside 30-50 is flagged") : fail("an invalid flap length is not flagged");
  await page.locator("#flap_length").fill("35");
  await page.locator("#max_offtargets").fill("11");
  await page.locator("#max_offtargets").blur();
  (await page.getByText("A whole number from 0 to 10.").count()) > 0 ? ok("an off-target allowance above 10 is flagged") : fail("an invalid off-target allowance is not flagged");
  await page.locator("#max_offtargets").fill("0");

  // A dry-run design, end to end.
  await designButton(/Design pegRNAs and donor/).click();
  await page.getByText("Design found").first().waitFor({ timeout: 120000 });
  sent.length === 1 ? ok("one request was sent") : fail(`${sent.length} requests were sent`);
  const body = sent[0] ?? {};
  body.strategy === "pa" && body.mode === "fake" ? ok("it names the strategy (pa) and the mode") : fail(`the request says ${body.strategy}/${body.mode}`);
  const cfg = body.config ?? {};
  cfg.flap_length === 35 && cfg.max_offtargets === 0 && cfg.max_pairs === 30 && cfg.min_span === 30 && cfg.max_span === 100 && Array.isArray(cfg.pbs_lengths)
    ? ok("it carries Prime Assembly's settings, named as the engine takes them")
    : fail(`the settings sent were ${JSON.stringify({ f: cfg.flap_length, o: cfg.max_offtargets, p: cfg.max_pairs })}`);
  ["rtt_aa", "pegrna_strand", "arm_side", "spacer", "ngrna_distance_min", "ngrna_max_offtargets", "het_snps"].some((k) => k in cfg)
    ? fail("the request carries tj-PE-only settings")
    : ok("it carries none of tj-PE's settings");
  cfg.cargo?.length === 499 && Array.isArray(cfg.cds_mask) && cfg.cds_mask.length === cfg.genome.length
    ? ok("the cargo and the coding mask are sent whole")
    : fail("the cargo or the mask is wrong");

  const headings = (await page.locator("h3, h4").allInnerTexts()).map((t) => t.trim());
  ["What to order", "pegRNA A", "pegRNA B", "Donor DNA", "Where it goes"].every((h) => headings.includes(h))
    ? ok("the result is an order sheet: both pegRNAs, the donor, and where it goes")
    : fail(`the result headings are ${JSON.stringify(headings)}`);
  (await page.getByText("This combination is untested").count()) === 1 ? ok("the untested-combination caveat has its own callout") : fail("the caveat callout is missing");
  // The notes come from the server, so the paper's wording lives in one place. The heading is matched exactly: the
  // plain-text report in the page holds the same words in capitals, and a substring match finds that as well.
  const prepTitle = page.getByText("How to prepare the donor", { exact: true });
  if ((await prepTitle.count()) === 1) {
    await prepTitle.click();
    const prep = page.locator("details", { has: prepTitle });
    const notes = await prep.locator("li").allInnerTexts();
    notes.length === 4 && /lambda exonuclease/.test(notes[1]) && /AZD-7648/.test(notes[2]) && /has to be measured/.test(notes[3])
      ? ok("how to prepare the donor lists the paper's three findings and says this design's variant is unmeasured")
      : fail(`the donor notes were ${JSON.stringify(notes.map((n) => n.slice(0, 40)))}`);
    (await prep.locator("a[href*='PMC12262537']").count()) === 1 ? ok("and they cite the paper") : fail("the donor notes do not cite the paper");
  } else {
    fail("no notes on preparing the donor");
  }
  const links = await page.locator("a[href*='/api/runs/']").evaluateAll((as) => as.map((a) => a.getAttribute("href")));
  const names = links.map((h) => h.split("/").pop());
  ["recoded_pa_report.txt", "recoded_pa_components.csv", "recoded_pa_donor.fasta", "recoded_pa_changes.csv"].every((n) => names.includes(n))
    ? ok("the downloads are the Prime Assembly files (report, components, FASTA, changes)")
    : fail(`the download links are ${JSON.stringify(names)}`);
  names.some((n) => n.startsWith("recoded_tjpe")) ? fail("a tj-PE download is offered for a Prime Assembly design") : ok("and no tj-PE download is offered");
  const cardText = await page.locator("main").innerText();
  /5.\s*35 nt match flap A/.test(cardText) && /3.\s*35 nt match flap B/.test(cardText)
    ? ok("the donor map marks the ends the two flaps match")
    : fail("the donor map does not mark the flap ends");
  const copies = await page.getByRole("button", { name: "Copy as RNA" }).count();
  copies === 2 ? ok("each pegRNA can be copied as RNA (T to U)") : fail(`${copies} copy-as-RNA buttons`);

  // Select a part of a pegRNA: its sequence is shown and copyable.
  await page.getByRole("button", { name: /^spacer A, 20 nucleotides/ }).click();
  (await page.getByRole("button", { name: "Copy sequence" }).count()) === 1 ? ok("selecting a part shows it with a copy button") : fail("selecting a part shows nothing");

  // A result belongs to the strategy that made it: switching shows tj-PE's pipeline idle and keeps the result.
  await radio(/Inside the pegRNA/).click();
  const idle = await page.locator("main").innerText();
  /Candidates run cheapest check first/.test(idle) && idle.includes("What to order")
    ? ok("switching strategy keeps the result and shows the other pipeline idle")
    : fail("switching strategy muddled the result and the pipeline");
  await ctx.close();
}

// tj-PE must be exactly what it was: its own example, its own request, its own result view.
console.log(`\nTJ-PE UNCHANGED  (${SITE})`);
{
  const { ctx, page } = await openApp({ os: "light" });
  const sent = [];
  page.on("request", (r) => r.url().includes("/api/design/stream") && sent.push(r.postDataJSON()));
  await page.getByRole("button", { name: "Load CTCF example" }).click();
  await page.waitForFunction(() => document.querySelector("#genome")?.value?.length > 100, { timeout: 30000 });
  (await page.getByRole("radio", { name: /Inside the pegRNA/ }).getAttribute("aria-checked")) === "true"
    ? ok("the tj-PE example leaves (or puts) the page on tj-PE")
    : fail("the tj-PE example did not select tj-PE");
  await page.getByRole("button", { name: "Design pegRNA", exact: true }).click();
  await page.getByText("Design found").first().waitFor({ timeout: 120000 });
  const b = sent[0] ?? {};
  b.strategy === "tjpe" && b.config?.pegrna_strand && !("flap_length" in b.config) && !("max_offtargets" in b.config)
    ? ok("the request is tj-PE's: its settings, none of Prime Assembly's")
    : fail(`the tj-PE request was ${JSON.stringify({ s: b.strategy, keys: Object.keys(b.config ?? {}) })}`);
  const hrefs = await page.locator("a[href*='/api/runs/']").evaluateAll((as) => as.map((a) => a.getAttribute("href").split("/").pop()));
  hrefs.includes("recoded_tjpe_report.txt") && hrefs.includes("recoded_tjpe_components.csv") && !hrefs.some((h) => h.startsWith("recoded_pa"))
    ? ok("the downloads are tj-PE's")
    : fail(`tj-PE downloads: ${JSON.stringify(hrefs)}`);
  (await page.getByText("pegRNA construct map").count()) === 1 && (await page.getByText("What to order").count()) === 0
    ? ok("the result is tj-PE's construct map, not the Prime Assembly order sheet")
    : fail("the tj-PE result view changed");
  await ctx.close();
}

// A server that predates Prime Assembly does not list it: the page says so and does not pretend.
console.log(`\nOLDER SERVER  (${SITE})`);
{
  const { ctx, page } = await openApp({ os: "light" });
  await page.route("**/api/health", (route) =>
    route.fulfill({
      status: 200,
      contentType: "application/json",
      headers: { "access-control-allow-origin": "*" },
      body: JSON.stringify({ status: "ok", version: "0.1.0", modes: ["fake", "api"], default_mode: "fake", repo_ok: true, api_python_ok: true, nvidia_key: true, concurrency: { limit: 2 } }),
    }),
  );
  await page.reload({ waitUntil: "domcontentloaded" });
  await page.getByRole("radio", { name: /separate donor DNA/ }).waitFor({ timeout: 30000 });
  await page.waitForTimeout(800);
  (await page.getByRole("radio", { name: /separate donor DNA/ }).isDisabled()) && (await page.getByText(/not been updated for Prime Assembly/).count()) > 0
    ? ok("a server without Prime Assembly gets a disabled option that says why")
    : fail("the page offers Prime Assembly on a server that does not have it");
  (await page.getByRole("radio", { name: /Inside the pegRNA/ }).isEnabled())
    ? ok("tj-PE is still available there")
    : fail("tj-PE is disabled on an older server");
  await ctx.close();
}

// A Prime Assembly run that finds no design explains why, from the pipeline's own record, and marks the step.
console.log(`\nPRIME ASSEMBLY, NO DESIGN  (${SITE})`);
{
  const { ctx, page } = await openApp({ os: "light" });
  // Strings the engine really writes (tests/test_pa_page_contract.py builds each kind with the real walk and checks
  // the page's patterns against them), plus one line nothing on the page recognises.
  const diagnostics = [
    "24 PAM-in pair(s) from 16 guide(s); 1 failed the sequence rules (GC content 25% is outside 30-75%), 11 failed the off-target screen (each needs a unique on-target and at most 0 off-target site(s) within the gate mismatch tier); 4 pair(s) remain",
    "pair 101/164: pegRNA A (PBS 13, flap 35): cross-region pairing (ViennaRNA): spacer[5..5] A x RTT[10..10] T (1 bp, p >= 0.999)",
    "pair 101/164: pegRNA A (PBS 15, flap 35): cross-region pairing (ViennaRNA): spacer[5..5] A x RTT[10..10] T (1 bp, p >= 0.999)",
    "pair 61/149: guide B (cut 149): no usable extension: TTTT (a Pol III terminator) in the RTT or PBS at every flap length 30-50",
    "pair 56/149: something this page has never heard of",
  ];
  await page.route("**/api/design/stream", (route) =>
    route.fulfill({
      status: 200,
      contentType: "application/x-ndjson",
      headers: { "access-control-allow-origin": "*", "access-control-allow-headers": "*" },
      body:
        JSON.stringify({ type: "result", ok: false, strategy: "pa", run_id: "ui-check", mode: "fake", returncode: 1, error: "",
          stderr: "design-pa failed: NoValidPADesign: no pair survived all gates after 1 tried", report: "", components: [], has_csv: false,
          has_fasta: false, docking_stub: false, het_edit: null,
          gates: { strategy: "pa", failed: true, error: "NoValidPADesign: no pair survived", diagnostics, stand_in_gates: [], scorer_mode: "fake" } }) + "\n",
    }),
  );
  await page.getByRole("button", { name: /Load cassette example/ }).click();
  await page.waitForFunction(() => document.querySelector("#genome")?.value?.length > 100, { timeout: 30000 });
  await page.getByRole("button", { name: /Design pegRNAs and donor/ }).click();
  await page.getByText("No design").first().waitFor({ timeout: 30000 });
  (await page.getByText("Why there is no design").count()) === 1 ? ok("a failed run says why") : fail("a failed run does not explain itself");
  const text = await page.locator("main").innerText();
  const counts = [
    [/Guide sequence rules:\s*1 guide rejected/, "1 guide by the sequence rules"],
    [/Off-target check:\s*11 guides rejected/, "11 guides by the off-target scan"],
    [/Write the flaps:\s*1 pair rejected/, "1 pair at the flaps"],
    [/pegRNA A folding:\s*2 pegRNAs rejected/, "2 pegRNAs by folding"],
    [/Other reasons:\s*1 pair rejected/, "1 pair for a reason the page does not know"],
  ];
  const missing = counts.filter(([re]) => !re.test(text)).map(([, what]) => what);
  missing.length === 0 ? ok("it counts every rejection by step, in the right units, and the unrecognised one as 'other'") : fail(`the rejection counts miss: ${missing.join("; ")}`);
  /Try more PBS lengths, a different flap length, or more pairs/.test(text) ? ok("it says what to try next, for the step that stopped the furthest candidate") : fail("no next step for the blocking step");
  const rowState = async (title) => (await page.locator("ol > li", { hasText: title }).first().locator("span.ml-auto").first().innerText()).trim();
  (await rowState("pegRNA A folding")) === "blocked" ? ok("the step that stopped the furthest candidate is marked blocked in the pipeline") : fail("the blocking step is not marked");
  const [before, after] = [await rowState("Write the flaps"), await rowState("pegRNA B folding")];
  before === "passed" && after === "not reached" ? ok("the steps before it are passed and the steps after it are not reached") : fail(`the steps around the blocked one read '${before}' and '${after}'`);
  (await page.getByText("Every rejection, in the order the pipeline recorded it").count()) === 1 ? ok("every recorded line is available") : fail("the full record is not available");
  await ctx.close();
}

// When the server refuses a request that passed the page's own checks (a limit the operator set, say), the page says
// the server's reason in its words and does not show a bare status code.
console.log(`\nSERVER REFUSAL  (${SITE})`);
{
  const { ctx, page } = await openApp({ os: "light" });
  const reason = "max_pairs may be at most 20 on this server (each pair can need docking)";
  await page.route("**/api/design/stream", (route) =>
    route.fulfill({
      status: 400,
      contentType: "application/json",
      headers: { "access-control-allow-origin": "*" },
      body: JSON.stringify({ detail: reason }),
    }),
  );
  await page.getByRole("button", { name: /Load cassette example/ }).click();
  await page.waitForFunction(() => document.querySelector("#genome")?.value?.length > 100, { timeout: 30000 });
  await page.getByRole("button", { name: /Design pegRNAs and donor/ }).click();
  await page.waitForTimeout(600);
  const alerts = (await page.locator("[role=alert]").allInnerTexts()).join(" | ");
  alerts.includes("The server rejected these inputs") && alerts.includes(reason)
    ? ok("a request the server refuses is answered with the server's own reason")
    : fail(`the refusal was shown as: ${alerts}`);
  // The generic advice for a rejected request is about the target window and site S, which is not what this refusal
  // is about, so it must not be tacked on after a specific reason.
  !alerts.includes("Check the target window and site S")
    ? ok("and it is not followed by generic advice about the target window")
    : fail(`a specific refusal was followed by generic advice: ${alerts}`);
  (await page.getByRole("button", { name: /Design pegRNAs and donor/ }).isEnabled())
    ? ok("and the page is ready to try again")
    : fail("the page is left stuck after a refusal");
  // The browser itself logs a console error for any 4xx response. This one was provoked on purpose, so it is not a
  // fault of the page; discount that one line and leave every other error to fail the run.
  const provoked = errors.findIndex((e) => /status of 400 \(Bad Request\)/.test(e));
  if (provoked >= 0) errors.splice(provoked, 1);
  await ctx.close();
}

// A Prime Assembly run in progress. The panel must show the step that is really running (not the instant one beside
// it), keep the insertion-site result while pairs are tried one after another, and drop a pair's results when the next
// pair starts. This needs a response that arrives one event at a time, so a small server here streams it by hand.
console.log(`\nPRIME ASSEMBLY, LIVE  (${SITE})`);
{
  const { createServer } = await import("node:http");
  let live = null;
  const headers = { "content-type": "application/x-ndjson", "access-control-allow-origin": "*", "access-control-allow-headers": "*" };
  const server = createServer((req, res) => {
    if (req.method === "OPTIONS") {
      res.writeHead(204, headers);
      res.end();
      return;
    }
    live = res;
    res.writeHead(200, headers);
    res.write("\n");
  });
  await new Promise((go) => server.listen(0, "127.0.0.1", go));
  const { ctx, page } = await openApp({ os: "light" });
  await page.route("**/api/design/stream", (route) => route.continue({ url: `http://127.0.0.1:${server.address().port}/stream` }));
  const send = (o) => new Promise((done) => live.write(JSON.stringify(o) + "\n", done));
  const ev = (cand, gate, passed = true) => ({ cand, gate, info: `${gate} ${passed ? "ok" : "rejected"}`, passed, metrics: {} });
  const state = async (title) => (await page.locator("ol > li", { hasText: title }).first().locator("span.ml-auto").first().innerText()).trim() || "idle";
  const states = async (...titles) => (await Promise.all(titles.map(state))).join(",");
  const header = async () => (await page.locator("h2", { hasText: "How your design is checked" }).locator("xpath=following-sibling::p").innerText()).trim();
  const settle = () => page.waitForTimeout(350);

  await page.getByRole("button", { name: /Load cassette example/ }).click();
  await page.waitForFunction(() => document.querySelector("#genome")?.value?.length > 100, { timeout: 30000 });
  await page.getByRole("button", { name: /Design pegRNAs and donor/ }).click();
  for (let i = 0; i < 100 && !live; i++) await page.waitForTimeout(100);
  if (!live) {
    fail("the page's design request never reached the streaming server");
  } else {
    await settle();
    (await states("Protein tolerance", "Off-target check")) === "running,idle" ? ok("before anything reports, the insertion-site check is the running step") : fail(`at the start: ${await states("Protein tolerance", "Off-target check")}`);

    await send(ev(1, "importance"));
    await settle();
    const scan = await states("Protein tolerance", "Pair the guides", "Guide sequence rules", "Off-target check", "Build the donor arms");
    scan === "passed,passed,passed,running,idle" ? ok("once the site is judged, the genome-wide guide scan is the running step, not the instant one beside it") : fail(`after the site check: ${scan}`);
    /screening the guides against the genome/.test(await header()) ? ok("and the header says what it is waiting for") : fail(`the header reads: ${await header()}`);

    await send(ev(1, "structure_a"));
    await send(ev(1, "structure_b"));
    await settle();
    const planned = await states("Off-target check", "Build the donor arms", "No deletion (INV-0)", "Write the flaps", "Strand check", "pegRNA A folding", "pegRNA B folding", "Cryptic splice sites");
    planned === "passed,passed,passed,passed,passed,passed,passed,running" ? ok("after the first folds, the planning steps are passed and the next check is running") : fail(`after the folds: ${planned}`);
    /trying guide pair 1/.test(await header()) ? ok("the header names the pair being tried") : fail(`the header reads: ${await header()}`);

    await send(ev(1, "cryptic"));
    await send(ev(1, "splice"));
    await send(ev(1, "docking_a", false));
    await settle();
    const rejected = await states("pegRNA A binds Cas9", "pegRNA B binds Cas9");
    rejected === "blocked,idle" ? ok("a pegRNA that fails to dock is blocked and the second is not started") : fail(`after the docking rejection: ${rejected}`);

    // Loading the other kind of example in the middle of a run changes the form's strategy. The panel that is showing
    // the run has to keep describing the run, not switch to the other strategy's steps under it.
    await page.getByRole("button", { name: "Load CTCF example" }).click();
    await settle();
    const stillPa = (await page.locator("ol > li", { hasText: "Guide sequence rules" }).count()) === 1;
    stillPa && (await state("pegRNA A binds Cas9")) === "blocked"
      ? ok("loading a tj-PE example mid-run leaves the panel on the run's own steps and states")
      : fail(`loading an example mid-run changed the panel (Prime Assembly steps present: ${stillPa})`);
    await page.getByRole("button", { name: /Load cassette example/ }).click();
    await settle();

    await send(ev(2, "structure_a"));
    await settle();
    const next = await states("Protein tolerance", "pegRNA A folding", "pegRNA B folding", "pegRNA A binds Cas9", "Cryptic splice sites");
    next === "passed,passed,running,idle,idle" ? ok("when the next pair starts, the first pair's results are gone and the insertion-site result stays") : fail(`after the second pair starts: ${next}`);
    /trying guide pair 2/.test(await header()) ? ok("and the header moves to pair 2") : fail(`the header reads: ${await header()}`);

    await send({ type: "result", ok: false, strategy: "pa", run_id: "ui-live", mode: "fake", returncode: 1, error: "", stderr: "design-pa failed: NoValidPADesign", report: "", components: [], has_csv: false, has_fasta: false, docking_stub: false, het_edit: null,
      gates: { strategy: "pa", failed: true, error: "NoValidPADesign: no pair survived", diagnostics: ["pair 101/164: pegRNA A (PBS 13, flap 35): docking: interface ipTM 0.31 is below the neutral extension's 0.58"], stand_in_gates: [], scorer_mode: "fake" } });
    live.end();
    await page.getByText("No design").first().waitFor({ timeout: 30000 });
    (await state("pegRNA A binds Cas9")) === "blocked" ? ok("when the run ends in failure, the pipeline marks the step that stopped it") : fail(`at the end: ${await state("pegRNA A binds Cas9")}`);
  }
  await ctx.close();
  server.close();
}

console.log("\nJS errors:", errors.length ? errors : "none");
if (errors.length) process.exitCode = 1;
await browser.close();
