// Offline render check for the collapsible atlas.
//   cd tools/project_map && npm install && node verify_render.mjs "$PWD/project_map.html"
// Verifies: no network, no console errors, sparse collapsed default, expansion
// grows the node count, focus dims non-neighbours, search reveals a match.
//
// Playwright resolution (portable — no machine-specific path; Issue #712):
//   1. PLAYWRIGHT_PATH env var, if set — an existing playwright(-core) package
//      *directory* (e.g. a global or npx-cache install). Resolved via the
//      package's own exports map, so the inner entry filename never matters.
//   2. otherwise a bare `playwright-core` import — resolved from this dir's
//      node_modules after `npm install` (package.json pins it).
// playwright-core (no bundled browsers) suffices because we launch *system*
// Chrome via channel:'chrome' below — no browser download needed.
import { pathToFileURL } from 'node:url';
import { createRequire } from 'node:module';

async function loadChromium() {
  const override = process.env.PLAYWRIGHT_PATH;
  if (override) {
    // createRequire resolves the target package via its exports/main. It is a
    // CommonJS require, so it loads the package's CJS entry (playwright-core
    // ships one); it cannot load an ESM-only package.
    try {
      const require = createRequire(pathToFileURL(override + '/').href);
      const { chromium } = require(override);
      if (!chromium) throw new Error('package has no "chromium" export');
      return chromium;
    } catch (e) {
      console.error(
        `PLAYWRIGHT_PATH="${override}" failed to resolve: ${e.message}\n` +
        'Point it at a playwright(-core) package directory.');
      process.exit(3);
    }
  }
  try {
    return (await import('playwright-core')).chromium;
  } catch {
    console.error(
      'Could not resolve Playwright. From a fresh clone:\n' +
      '  cd tools/project_map && npm install && node verify_render.mjs "$PWD/project_map.html"\n' +
      'or set PLAYWRIGHT_PATH to an existing playwright(-core) install directory.');
    process.exit(3);
  }
}
const chromium = await loadChromium();

const url = pathToFileURL(process.argv[2]).href;
const browser = await chromium.launch({ channel: 'chrome', headless: true });
const ctx = await browser.newContext();
const blocked = [];
await ctx.route('**/*', r => r.request().url().startsWith('file:')
  ? r.continue() : (blocked.push(r.request().url()), r.abort()));
await ctx.setOffline(true);
const page = await ctx.newPage();
const errors = [];
page.on('console', m => { if (m.type() === 'error') errors.push(m.text()); });
page.on('pageerror', e => errors.push(e.message));

await page.goto(url, { waitUntil: 'load', timeout: 30000 });
await page.waitForTimeout(2000);   // let the force sim settle

const circ = () => page.$$eval('#svg-container svg circle', e => e.length).catch(() => 0);
// Dispatch a click on the SVG node matching `sel` (bypasses Playwright's
// actionability checks — the force sim keeps nodes moving + the svg intercepts).
const clickNode = sel => page.evaluate(s => {
  const el = document.querySelector(s);
  if (el) el.dispatchEvent(new MouseEvent('click', { bubbles: true, cancelable: true }));
  return !!el;
}, sel);

const collapsed = await circ();

// Expand: click the first collapsed-container node group.
await clickNode('.node-group[data-collapsed="1"]');
await page.waitForTimeout(900);
const expandedCount = await circ();

// Focus: click a true leaf (not a container, not a root) and confirm dimming.
await clickNode('.node-group:not([data-container]):not([data-node-type="root"])');
await page.waitForTimeout(500);
const dimmed = await page.$$eval('#svg-container svg .node-group circle',
  els => els.filter(e => parseFloat(e.getAttribute('opacity') || '1') < 0.5).length).catch(() => 0);

// Search reveals a match (expanding ancestors).
await page.fill('#search-input', 'mhcflurry');
await page.waitForTimeout(900);
const afterSearch = await circ();

const failed = await page.$eval('#svg-container', el => el.innerText.includes('Failed to load')).catch(() => false);
const d3type = await page.evaluate(() => typeof window.d3);
const issueOpts = await page.$$eval('#issue-filter option', o => o.length).catch(() => 0);
await browser.close();

const pass = blocked.length === 0 && errors.length === 0 && !failed && d3type === 'object'
  && collapsed > 0 && collapsed <= 20 && expandedCount > collapsed && dimmed > 0 && afterSearch > 0 && issueOpts > 1;
console.log(JSON.stringify({ blocked: blocked.length, errors, d3: d3type, collapsed,
  expandedCount, dimmedOnFocus: dimmed, afterSearch, issueOptions: issueOpts, failed, pass }, null, 2));
process.exit(pass ? 0 : 2);
