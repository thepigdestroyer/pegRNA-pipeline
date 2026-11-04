/**
 * Quick UI checks: no token, no design run, a few seconds.
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

/** Open the app with an OS colour scheme, and optionally a previously saved theme choice. */
async function openApp({ os, saved }) {
  const ctx = await browser.newContext({ colorScheme: os, viewport: { width: 1280, height: 900 } });
  if (saved) await ctx.addInitScript((t) => localStorage.setItem("pegrna.theme", t), saved);
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

console.log("\nJS errors:", errors.length ? errors : "none");
if (errors.length) process.exitCode = 1;
await browser.close();
