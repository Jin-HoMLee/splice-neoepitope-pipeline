# Project Map Atlas Redesign — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the project map from a force-directed hairball into a legible, collapsible whole-project atlas — fix the data model (classification, real cross-links, issue tags) and stop drawing folder hierarchy as edges.

**Architecture:** Two units. (1) `tools/project_map/extract_graph.py` produces correct, richly-connected `graph.json`. (2) `tools/project_map/project_map.html` renders a collapsible hierarchy — containment expressed by expand/collapse + position, only cross-cutting edges drawn, with focus/search/filter. `build_html.py` (inlines data + vendored D3) is unchanged.

**Tech Stack:** Python 3 (stdlib only) for extraction; D3 v7 (vendored, inlined) for the viz; Playwright + system Chrome for offline render verification.

**Spec:** [`docs/superpowers/specs/2026-06-06-project-map-atlas-redesign-design.md`](../specs/2026-06-06-project-map-atlas-redesign-design.md)

**Test runner:** `tools/project_map/test_extract_graph.py`, run with the project's pytest venv:
`workflow/tests/.venv/bin/python -m pytest tools/project_map/test_extract_graph.py -v`
(falls back to any `python -m pytest`; tests import only stdlib + `extract_graph`).

---

## File Structure

| File | Responsibility | Action |
|---|---|---|
| `tools/project_map/extract_graph.py` | repo → graph data (classification, cross-edges, issue tags, guards) | Modify |
| `tools/project_map/test_extract_graph.py` | unit + invariant tests for the extractor | Create |
| `tools/project_map/project_map.html` | collapsible-hierarchy render + interaction | Modify (template region only) |
| `tools/project_map/graph.json` | generated data | Regenerate (don't hand-edit) |
| `tools/project_map/verify_render.mjs` | offline Playwright render check | Create |

---

## Phase 1 — Data model (`extract_graph.py`), TDD

### Task 1: Fix `classify_path` (the leading-slash bug)

**Files:**
- Modify: `tools/project_map/extract_graph.py` (the `classify_path` function)
- Test: `tools/project_map/test_extract_graph.py`

- [ ] **Step 1: Write the failing test**

```python
# tools/project_map/test_extract_graph.py
import importlib.util, pathlib
_spec = importlib.util.spec_from_file_location(
    "extract_graph", pathlib.Path(__file__).parent / "extract_graph.py")
eg = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(eg)


def test_classify_path_relative_paths():
    # These are the relative paths the os.walk actually passes (no leading slash).
    assert eg.classify_path("research/manuscript/METHODS.md") == "manuscript"
    assert eg.classify_path("docs/google_cloud_guide.md") == "doc"
    assert eg.classify_path("research/evals/issue_188_boltz2/slides.qmd") == "slide"
    assert eg.classify_path("research/lab_notebook/pm.md") == "labnotebook"
    assert eg.classify_path("workflow/scripts/run_mhcflurry.py") == "script"
    assert eg.classify_path("workflow/envs/python.yaml") == "env"
    assert eg.classify_path("config/config.yaml") == "config"
    assert eg.classify_path("tools/ci/closure_audit.py") == "ci"
    assert eg.classify_path("research/notebooks/patient_001_results.ipynb") == "notebook"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest tools/project_map/test_extract_graph.py::test_classify_path_relative_paths -v`
Expected: FAIL — several asserts return `'resource'` (e.g. `METHODS.md`).

- [ ] **Step 3: Implement the fix**

In `classify_path`, the `is_dir=True` branch (exact `==` / `startswith`) is correct for relative paths and stays. Replace the **file branch** so its directory-substring checks operate on a slash-prefixed copy, and basename checks use the basename. Find the line after the `if is_dir:` block ends (the first `if p.endswith(".smk")`) and replace the whole file branch with:

```python
    # ── file classification ──
    # os.walk passes relative paths ("research/x.md"); the directory tests below
    # are written for absolute-style ("/research/"), so normalise with a leading
    # slash. basename checks use `base`.
    q = "/" + p
    base = p.rsplit("/", 1)[-1]

    if q.endswith(".smk"): return "rule"
    if q.endswith((".py", ".sh")) and "/workflow/scripts/" in q: return "script"
    if q.endswith(".py") and "/scripts/" in q and "/pm/" not in q: return "script"
    if q.endswith(".yaml") and "/workflow/envs/" in q: return "env"
    if q.endswith(".yaml") and "/config/" in q: return "config"
    if q.endswith(".tsv") and "/config/samples/" in q: return "sample"
    if q.endswith(".tsv"): return "data"
    if q.endswith(".md") and "/docs/" in q: return "doc"
    if q.endswith(".md") and "/research/manuscript/" in q: return "manuscript"
    if q.endswith(".md") and "/research/lab_notebook/" in q: return "labnotebook"
    if q.endswith("/lab_notebook.md"): return "labnotebook"
    if q.endswith(".ipynb") and "/research/notebooks/" in q: return "notebook"
    if q.endswith(".qmd") and ("/research/evals/" in q or "/research/experiments/" in q
                               or "/research/decisions/" in q or "/research/slides/" in q
                               or "/docs/features/" in q): return "slide"
    if q.endswith(".sh") and "/scripts/" in q: return "script"
    if q.endswith(".py") and "/tools/ci/" in q: return "ci"
    if q.endswith(".py") and "/tools/news/" in q: return "news"
    if q.endswith(".py") and "/workflow/tests/" in q: return "test"
    if q.endswith(".py") and "/scripts/tests/" in q: return "test"
    if q.endswith(".ini"): return "config"
    if base == "snakefile": return "rule"
    if base == "dockerfile.pipeline": return "docker"
    if q.endswith(".json") and "/scripts/" in q: return "config"
    if q.endswith((".csl", ".scss")): return "resource"
    if q.endswith(".bib"): return "resource"
    if base == "readme.md": return "doc"
    if q.endswith(".yaml") and "/tools/news/" in q: return "config"
    return "resource"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest tools/project_map/test_extract_graph.py::test_classify_path_relative_paths -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add tools/project_map/extract_graph.py tools/project_map/test_extract_graph.py
git commit -m "fix(project-map): classify_path matches relative walk paths (kills resource blob)"
```

---

### Task 2: Parse markdown→markdown links into `references` edges

**Files:**
- Modify: `tools/project_map/extract_graph.py` (add `resolve_link`, `parse_markdown_links`; call in `build_graph`)
- Test: `tools/project_map/test_extract_graph.py`

- [ ] **Step 1: Write the failing test (link resolution)**

```python
def test_resolve_link():
    assert eg.resolve_link("docs/a.md", "b.md") == "docs/b.md"
    assert eg.resolve_link("docs/a.md", "sub/c.md#section") == "docs/sub/c.md"
    assert eg.resolve_link("research/evals/x/README.md",
                           "../../slides/nature.csl") == "research/slides/nature.csl"
    assert eg.resolve_link("docs/a.md", "https://x.com/y.md") is None  # external
    assert eg.resolve_link("docs/a.md", "#anchor-only") is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest tools/project_map/test_extract_graph.py::test_resolve_link -v`
Expected: FAIL — `AttributeError: module 'extract_graph' has no attribute 'resolve_link'`.

- [ ] **Step 3: Implement `resolve_link` + `parse_markdown_links`**

Add near the other helpers (after `get_pipeline_env_deps`):

```python
MD_LINK_RE = re.compile(r'\]\(\s*([^)\s]+\.md)(?:#[^)]*)?\s*\)')


def resolve_link(src_rel, link):
    """Resolve a markdown link target to a repo-relative path, or None if external/anchor."""
    link = link.strip()
    if not link or link.startswith(("http://", "https://", "#", "mailto:")):
        return None
    link = link.split("#", 1)[0]
    if not link:
        return None
    joined = os.path.normpath(os.path.join(os.path.dirname(src_rel), link))
    return joined.replace(os.sep, "/")


def parse_markdown_links(node_ids):
    """Scan docs/ + research/ markdown for links to other repo .md files → references edges."""
    edges = []
    for base in ("docs", "research"):
        base_dir = PROJECT_ROOT / base
        if not base_dir.exists():
            continue
        for md in base_dir.rglob("*.md"):
            src_rel = str(md.relative_to(PROJECT_ROOT))
            src_id = slug(src_rel)
            if src_id not in node_ids:
                continue
            try:
                text = md.read_text()
            except Exception:
                continue
            for m in MD_LINK_RE.finditer(text):
                tgt_rel = resolve_link(src_rel, m.group(1))
                if not tgt_rel:
                    continue
                tgt_id = slug(tgt_rel)
                if tgt_id in node_ids and tgt_id != src_id:
                    edges.append({"source": src_id, "target": tgt_id,
                                  "type": "references", "label": "links to"})
    return edges
```

In `build_graph`, after the test-coverage edges block and **before** the degrees
section, add:

```python
    # ── Research/docs cross-links (markdown → markdown) ──
    edges.extend(parse_markdown_links({n["id"] for n in nodes}))
```

(The defense-in-depth dangling-edge filter that follows will drop any link whose
target isn't a node, so this is safe even if a link points outside the graph.)

- [ ] **Step 4: Run resolution test (passes) + add an integration assertion**

Run: `workflow/tests/.venv/bin/python -m pytest tools/project_map/test_extract_graph.py::test_resolve_link -v`
Expected: PASS.

Add this integration test (it walks the real repo, ~1s):

```python
def test_research_half_is_connected():
    g = eg.build_graph()
    ids = {n["id"] for n in g["nodes"]}
    research_ids = {n["id"] for n in g["nodes"] if n["group"] == "research"}
    real = [e for e in g["edges"]
            if e["type"] != "contained_in"
            and (e["source"] in research_ids or e["target"] in research_ids)]
    assert len(real) > 10, f"research still nearly disconnected: {len(real)} real edges"
    # every edge endpoint is a real node (D3 forceLink safety)
    assert all(e["source"] in ids and e["target"] in ids for e in g["edges"])
```

Run: `workflow/tests/.venv/bin/python -m pytest tools/project_map/test_extract_graph.py::test_research_half_is_connected -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add tools/project_map/extract_graph.py tools/project_map/test_extract_graph.py
git commit -m "feat(project-map): parse markdown cross-links into references edges"
```

---

### Task 3: Issue tagging from `issue_NNN/` folders

**Files:**
- Modify: `tools/project_map/extract_graph.py`
- Test: `tools/project_map/test_extract_graph.py`

- [ ] **Step 1: Write the failing test**

```python
def test_issue_tagging():
    g = eg.build_graph()
    by_path = {n.get("path"): n for n in g["nodes"] if n.get("path")}
    hit = [n for p, n in by_path.items() if "issue_218" in p and n.get("issue") == 218]
    assert hit, "expected nodes under issue_218_* to carry issue==218"
    # nodes outside any issue folder must not be tagged
    assert by_path.get("workflow/scripts/run_mhcflurry.py", {}).get("issue") is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `workflow/tests/.venv/bin/python -m pytest tools/project_map/test_extract_graph.py::test_issue_tagging -v`
Expected: FAIL — no `issue` key.

- [ ] **Step 3: Implement**

Add a helper near `slug`:

```python
ISSUE_RE = re.compile(r'(?:^|/)issue_(\d+)')


def issue_of(rel_path):
    m = ISSUE_RE.search(rel_path)
    return int(m.group(1)) if m else None
```

In `build_graph`, the directory-walk file loop builds each file node dict. Find
where it appends the walked file node (the `nodes.append({... "path": rel_path ...})`
in the walk) and set the issue tag right after computing `rel_path`:

```python
            _issue = issue_of(rel_path)
```

then add `**({"issue": _issue} if _issue else {})` into that node dict, e.g.:

```python
            nodes.append({
                "id": file_id, "label": f, "type": ftype,
                "group": group, "path": rel_path,
                "description": describe_file(rel_path, ftype), "size": sz,
                **({"issue": _issue} if _issue else {}),
            })
```

Do the same for the walked **directory** node dict (use `issue_of(rel_root)`).

- [ ] **Step 4: Run test to verify it passes**

Run: `workflow/tests/.venv/bin/python -m pytest tools/project_map/test_extract_graph.py::test_issue_tagging -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add tools/project_map/extract_graph.py tools/project_map/test_extract_graph.py
git commit -m "feat(project-map): tag nodes under issue_NNN folders with issue id"
```

---

### Task 4: Classification-health invariant + regenerate

**Files:**
- Modify: `tools/project_map/extract_graph.py` (`main` stats)
- Test: `tools/project_map/test_extract_graph.py`

- [ ] **Step 1: Write the failing test**

```python
def test_resource_blob_is_gone():
    g = eg.build_graph()
    from collections import Counter
    by_group = {}
    for n in g["nodes"]:
        by_group.setdefault(n["group"], []).append(n["type"])
    for grp, types in by_group.items():
        frac = types.count("resource") / len(types)
        assert frac < 0.40, f"group '{grp}' is {frac:.0%} unclassified 'resource'"
```

- [ ] **Step 2: Run test to verify it fails (before Task 1) / passes (after)**

Run: `workflow/tests/.venv/bin/python -m pytest tools/project_map/test_extract_graph.py::test_resource_blob_is_gone -v`
Expected after Tasks 1–3: PASS. (Documents the regression guard.)

- [ ] **Step 3: Add a one-line classification summary to `main()`**

After the existing stats prints in `main()`, add:

```python
    from collections import Counter
    res = sum(1 for n in graph["nodes"] if n["type"] == "resource")
    print(f"  Unclassified 'resource' nodes: {res}/{len(graph['nodes'])} "
          f"({100*res//len(graph['nodes'])}%)")
```

- [ ] **Step 4: Regenerate + run full test file**

Run:
```bash
conda activate snakemake
python tools/project_map/extract_graph.py
workflow/tests/.venv/bin/python -m pytest tools/project_map/test_extract_graph.py -v
```
Expected: all tests PASS; `resource` count drops sharply from 236.

- [ ] **Step 5: Commit**

```bash
git add tools/project_map/extract_graph.py tools/project_map/test_extract_graph.py tools/project_map/graph.json
git commit -m "test(project-map): guard against resource-blob regression; regenerate graph"
```

---

## Phase 2 — Visualization (`project_map.html`), collapsible hierarchy

> Verification oracle for Phase 2 is the **offline Playwright render** (Task 9),
> not unit tests — this is interactive D3. Build incrementally and re-run the
> render check after each task. Edit only the `<script>` region of
> `project_map.html`; rebuild with `python tools/project_map/build_html.py` after
> each change (it re-inlines data + D3).

### Task 5: Build a hierarchy model; stop drawing `contained_in`

**Files:** Modify `tools/project_map/project_map.html` (main `<script>`)

- [ ] **Step 1: Add hierarchy construction** after `graphData` is parsed in `loadGraph`/`init`. Containment edges are `{source: child, target: parent}` for file/dir nodes and `{source: parent, target: child}` for the 5 root containments — normalise both into a `parentOf` map keyed by child id, treating the 5 group roots as top-level (parent = null).

```javascript
// ── Hierarchy from contained_in edges (group roots are top-level) ──
const GROUP_ROOTS = new Set(graphData.nodes.filter(n => n.type === "root").map(n => n.id));
const childrenOf = new Map();   // id -> [childId]
const parentOf = new Map();     // id -> parentId
function buildHierarchy() {
  childrenOf.clear(); parentOf.clear();
  for (const e of graphData.edges) {
    if (e.type !== "contained_in") continue;
    // Orient: the endpoint that is NOT a group root (or is deeper) is the child.
    let child = e.source, parent = e.target;
    if (GROUP_ROOTS.has(e.source)) { child = e.target; parent = e.source; }
    parentOf.set(child, parent);
    if (!childrenOf.has(parent)) childrenOf.set(parent, []);
    childrenOf.get(parent).push(child);
  }
}
function isContainer(id) { return childrenOf.has(id); }
function ancestors(id) { const out = []; let p = parentOf.get(id);
  while (p != null) { out.push(p); p = parentOf.get(p); } return out; }
```

- [ ] **Step 2: Drop `contained_in` from drawn edges.** In `renderGraph`/`filterByView`, the drawn edge set must exclude `contained_in` (hierarchy is shown by expansion, not lines). Confirm via render that no faint grey containment lines remain.

- [ ] **Step 3: Rebuild + eyeball**

```bash
python tools/project_map/build_html.py
```
Open the file; the blob should already be far sparser (only ~84+ cross edges drawn).

- [ ] **Step 4: Commit**

```bash
git add tools/project_map/project_map.html
git commit -m "feat(project-map): model folder hierarchy; stop drawing containment edges"
```

### Task 6: Collapsible state — start at group level, expand on click

**Files:** Modify `tools/project_map/project_map.html`

- [ ] **Step 1:** Add an `expanded` Set (initially the group roots are *visible* but collapsed → `expanded` empty). Add `visibleNodes()`: a node is visible iff every ancestor is in `expanded`. A container that is visible but not expanded renders as a single aggregate node (with a `+` affordance + child count).

```javascript
const expanded = new Set();   // ids of containers whose children are shown
function visibleNodes() {
  return graphData.nodes.filter(n => {
    if (GROUP_ROOTS.has(n.id)) return true;            // groups always on screen
    return ancestors(n.id).every(a => expanded.has(a)); // else all ancestors expanded
  });
}
function toggle(id) {
  if (!isContainer(id)) return;
  expanded.has(id) ? collapse(id) : expanded.add(id);
  renderGraph();
}
function collapse(id) {                 // collapse id and all descendants
  expanded.delete(id);
  for (const c of (childrenOf.get(id) || [])) collapse(c);
}
```

- [ ] **Step 2:** Wire `toggle(d.id)` into the node click handler for containers (leaves keep the existing detail-panel behaviour). Re-run the force sim over `visibleNodes()` only.

- [ ] **Step 3:** Add a breadcrumb + "Collapse all" button (`expanded.clear(); renderGraph()`), and replace the old Full/Pipeline/Research tabs' behaviour to expand the corresponding group(s) instead of filtering a flat list.

- [ ] **Step 4:** Rebuild, then verify default shows ≤ ~15 nodes and clicking a group reveals children.

```bash
python tools/project_map/build_html.py
```

- [ ] **Step 5: Commit**

```bash
git add tools/project_map/project_map.html
git commit -m "feat(project-map): collapsible groups; force sim over visible frontier"
```

### Task 7: Aggregate cross-edges to the nearest visible ancestor

**Files:** Modify `tools/project_map/project_map.html`

- [ ] **Step 1:** Add `visibleEdges()`: for each non-`contained_in` edge, "lift" each endpoint to its nearest visible ancestor (itself if visible), drop self-loops, and aggregate duplicates into a single link carrying a `count`.

```javascript
function liftToVisible(id, visibleSet) {
  if (visibleSet.has(id)) return id;
  for (const a of ancestors(id)) if (visibleSet.has(a)) return a;
  return null;
}
function visibleEdges() {
  const vis = new Set(visibleNodes().map(n => n.id));
  const agg = new Map();   // "s|t|type" -> {source,target,type,count}
  for (const e of graphData.edges) {
    if (e.type === "contained_in") continue;
    const s = liftToVisible(e.source, vis), t = liftToVisible(e.target, vis);
    if (s == null || t == null || s === t) continue;
    const k = s + "|" + t + "|" + e.type;
    if (!agg.has(k)) agg.set(k, {source: s, target: t, type: e.type, count: 0});
    agg.get(k).count++;
  }
  return [...agg.values()];
}
```

- [ ] **Step 2:** Use `visibleEdges()` in `renderGraph`; set link stroke-width ∝ `count`, and show `count` (when >1) as an edge label or thickness. Collapsed groups now show "Research →N→ Pipeline" style aggregate links.

- [ ] **Step 3:** Rebuild + verify aggregated links appear collapsed and explode on expand.

- [ ] **Step 4: Commit**

```bash
git add tools/project_map/project_map.html
git commit -m "feat(project-map): aggregate cross-edges to nearest visible ancestor"
```

### Task 8: Focus mode, search-to-expand, filters, legend

**Files:** Modify `tools/project_map/project_map.html`

- [ ] **Step 1: Focus mode** — clicking a leaf highlights its real-edge neighbours (over the full edge set, lifting to visible) and dims the rest; click background to clear.
- [ ] **Step 2: Search-to-expand** — on a search hit, `expanded`-add all of the hit's `ancestors(id)` then `renderGraph()` and pan to it (reuse existing search box).
- [ ] **Step 3: Filters** — group + type + a new **issue** filter (populate from distinct `node.issue` values; selecting one dims/hides non-matching). Colour-by-issue toggle optional.
- [ ] **Step 4: Legend** — regenerate legend entries from the node `type`s actually present (now meaningful post-Task 1) with their colours; show the group palette.
- [ ] **Step 5:** Rebuild + verify each interaction.
- [ ] **Step 6: Commit**

```bash
git add tools/project_map/project_map.html
git commit -m "feat(project-map): focus mode, search-to-expand, issue filter, real legend"
```

### Task 9: Offline render verification

**Files:** Create `tools/project_map/verify_render.mjs`

- [ ] **Step 1: Write the render check** (system Chrome via Playwright, network fully blocked):

```javascript
// tools/project_map/verify_render.mjs  — run: node verify_render.mjs <abs path to project_map.html>
import pw from '/Users/jin-holee/.npm/_npx/e41f203b7505f1fb/node_modules/playwright/index.js';
import { pathToFileURL } from 'node:url';
const { chromium } = pw;
const url = pathToFileURL(process.argv[2]).href;
const browser = await chromium.launch({ channel: 'chrome', headless: true });
const ctx = await browser.newContext();
const blocked = [];
await ctx.route('**/*', r => r.request().url().startsWith('file:') ? r.continue()
                                                                   : (blocked.push(r.request().url()), r.abort()));
await ctx.setOffline(true);
const page = await ctx.newPage();
const errors = [];
page.on('console', m => { if (m.type()==='error') errors.push(m.text()); });
page.on('pageerror', e => errors.push(e.message));
await page.goto(url, { waitUntil: 'load' });
await page.waitForTimeout(1200);
const collapsed = await page.$$eval('#svg-container svg circle', e => e.length).catch(()=>0);
// expand the first group and re-count
const firstGroup = await page.$('[data-node-type="root"], #svg-container svg circle');
if (firstGroup) { await firstGroup.click(); await page.waitForTimeout(800); }
const expanded = await page.$$eval('#svg-container svg circle', e => e.length).catch(()=>0);
const failed = await page.$eval('#svg-container', el => el.innerText.includes('Failed to load')).catch(()=>false);
await browser.close();
const pass = blocked.length===0 && errors.length===0 && !failed && collapsed>0 && collapsed<=20 && expanded>collapsed;
console.log(JSON.stringify({ blocked: blocked.length, errors, collapsed, expanded, failed, pass }, null, 2));
process.exit(pass ? 0 : 2);
```

- [ ] **Step 2: Run it**

```bash
python tools/project_map/extract_graph.py && python tools/project_map/build_html.py
node tools/project_map/verify_render.mjs "$(pwd)/tools/project_map/project_map.html"
```
Expected: `pass: true` — 0 blocked network requests, 0 console errors, collapsed default ≤ 20 nodes, expanding a group increases the node count, no "Failed to load".

- [ ] **Step 3: Commit**

```bash
git add tools/project_map/verify_render.mjs tools/project_map/project_map.html tools/project_map/graph.json
git commit -m "test(project-map): offline render check for collapsible atlas"
```

---

## Self-Review

**Spec coverage:** Layout/collapsible (Tasks 5–7) ✓; navigation focus/search/filter/legend (Task 8) ✓; classify_path fix (Task 1) ✓; research cross-links (Task 2) ✓; issue tags (Task 3) ✓; data invariants (Tasks 2,4 + existing dangling/dup guards) ✓; offline render testing (Task 9) ✓. Out-of-scope items (import call-graph, issue-tracker) correctly absent.

**Type consistency:** `childrenOf`/`parentOf`/`expanded`/`visibleNodes()`/`visibleEdges()`/`liftToVisible()`/`ancestors()`/`isContainer()`/`toggle()`/`collapse()` are used consistently across Tasks 5–8. Python: `resolve_link`, `parse_markdown_links`, `issue_of`, `classify_path` consistent across tasks/tests.

**Placeholders:** none — every code step shows real code.

**Note for executor:** the Playwright path in Task 9 (`/Users/jin-holee/.npm/_npx/...`) is machine-specific; if absent, locate playwright via `find ~/.npm/_npx -path '*node_modules/playwright' -type d` and substitute, or `npx playwright`.
