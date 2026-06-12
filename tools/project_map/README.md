# Project Map

An interactive **2026 SOTA project map** — a collapsible D3 atlas of the whole
repository: pipeline rules/scripts/envs, the research program, infrastructure,
docs, and config, with the real dependency/dataflow/reference edges between them.

Open [`project_map.html`](project_map.html) directly in any browser. It is
**fully self-contained** — the graph data and the D3 library are both inlined, so
it needs no web server and no network connection (share it, email it, open it off
a USB stick — it just works).

## How to read it

Folder hierarchy is **structure, not edges** — you start zoomed-out at the 5 group
roots (Pipeline, Research, Infrastructure, Docs, Config) and their top-level
folders, and the only lines drawn are *cross-cutting relationships* (calls,
produces, depends-on, references, tests, …). So the graph shows how things relate,
not 300 "lives-in-this-folder" lines.

- **Click a folder** (`+`/`−`) to expand/collapse it; collapsed folders aggregate
  their crossing edges (e.g. "Research →3 references→ Pipeline").
- **Click a file** to focus its real relationships (callers, deps, tests) and dim
  the rest; the side panel lists them and lets you jump.
- **Search** expands the tree to reveal matches; the **issue filter** dims to one
  `issue_NNN`; the **Pipeline/Research** tabs jump to a pre-expanded slice; **⊖**
  collapses back to the overview.

## Files

| File | Role |
|------|------|
| `extract_graph.py` | Introspects the repo tree → `graph.json` (nodes + edges + issue tags). |
| `build_html.py` | Inlines `graph.json` + vendored D3 into `project_map.html`. |
| `graph.json` | Generated graph data (regenerate; don't hand-edit). |
| `project_map.html` | The artifact you open — **also the template** (hand-authored D3/CSS/interaction code lives here). |
| `vendor/d3.v7.min.js` | Pinned D3 v7.9.0, inlined at build time. |
| `test_extract_graph.py` | Pytest unit + invariant tests for the extractor. |
| `verify_render.mjs` | Offline Playwright render check (collapsible behaviour, no network). |

## Test

```bash
# data model (classification, cross-links, issue tags, invariants)
workflow/tests/.venv/bin/python -m pytest tools/project_map/test_extract_graph.py -v
# rendering (offline, system Chrome via Playwright)
node tools/project_map/verify_render.mjs "$(pwd)/tools/project_map/project_map.html"
```

The data-model suite is collected in CI by the `ci-tools-pytest` job
(`.github/workflows/tests.yml`, Issue #713) and runs on every PR. `build_graph()`
walks the committed tree without consulting `.gitignore`, so a clean checkout is
the canonical run environment: a local clone with populated gitignored artifacts
(`references/`, `logs/`, `data/`, `results/`) can inflate the `project` group and
trip `test_resource_blob_is_gone`. CI (a fresh checkout) is authoritative; locally,
re-run against `git worktree add --detach <tmp> HEAD` if that one test fails alone.
(The `verify_render.mjs` render check is not yet CI-wired — Issue #712.)

## Regenerate

Run from anywhere (the scripts resolve paths relative to themselves):

```bash
python tools/project_map/extract_graph.py   # repo tree → graph.json
python tools/project_map/build_html.py       # graph.json + D3 → project_map.html
```

Both scripts are **idempotent** — safe to re-run. `build_html.py` strips any prior
inlined data/D3 blocks before re-inserting, so the output never accumulates
duplicates.

## Notes

- `project_map.html` is the canonical source *and* the output: `build_html.py`
  edits it in place (swapping the inlined `<script id="graph-data">` and
  `<script id="d3-inlined">` blocks). Edit the D3/CSS/interaction code directly in
  this file; only the two inlined blocks are machine-managed.
- `graph.json`, `project_map.html`, `vendor/`, and Quarto's rendered slide output
  (`*_files/`) are excluded from the map itself — build artifacts and vendored
  deps aren't "project structure".
- `extract_graph.py` also: parses markdown→markdown links into `references` edges
  (so the research half is connected), tags nodes under `issue_NNN/` folders with
  their issue id, and builds a complete single-parent containment tree (group root
  → top folder → subfolder → file).
- `extract_graph.py` drops any edge whose endpoint has no node and de-duplicates
  node ids before writing, so a malformed edge can never blank the graph (D3's
  `forceLink` throws `missing: <id>` on dangling edges).
- To bump D3: replace `vendor/d3.v7.min.js` and re-run `build_html.py`.
