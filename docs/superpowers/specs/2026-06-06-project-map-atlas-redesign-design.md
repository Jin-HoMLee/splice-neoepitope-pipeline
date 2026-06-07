# Project Map — "Legible Whole-Project Atlas" Redesign

**Date:** 2026-06-06
**Status:** Approved design → implementation plan next
**Component:** [`tools/project_map/`](../../../tools/project_map/) (`extract_graph.py`, `build_html.py`, `project_map.html`)

## Goal

Turn the 2026 SOTA Project Map from a force-directed hairball into a **legible
whole-project atlas** — something where you can see *what's in the project* and
*how the parts relate* at a glance, and navigate 393 nodes without getting lost.

## Problem (from the 2026-06-06 audit)

| Finding | Evidence |
|---|---|
| It's a folder tree drawn as a force blob | **78% of nodes (309/393)** have no edge except `contained_in`; 268 are single-edge leaves; **309 of 393 edges are `contained_in`** |
| `classify_path` is broken (leading-slash bug) | **60% of nodes (236)** fall into the catch-all `resource` type — `METHODS.md → 'resource'`, but `/METHODS.md → 'manuscript'` |
| The research half is unmodeled | 170 nodes (43% of the map), only **10 real edges** |
| Real structure sits unused | **157** markdown→markdown links, **13** `issue_NNN/` folders + ~3,700 issue refs, **80** script imports — all modeled as **0** edges |

Two distinct problems: **(1)** a data model that is deep for the pipeline and a
flat file-dump everywhere else, and **(2)** a visual form (force graph of a
hierarchy) that fights comprehension.

## Core principle

**Containment is structure, not an edge.** Express folder hierarchy through
collapsible nesting + position; **draw only cross-cutting relationships**
(`produces`, `calls`, `depends_on`, `references`, `test_of`, `runs_in`,
`includes`, `defined_in`). This removes ~309 noise edges and surfaces the ~84
real ones (more, once research links are parsed).

## Design

### A. Layout — collapsible hierarchy (start zoomed-out)

- **Default state:** the 5 group containers (Pipeline, Research, Infrastructure,
  Docs, Config), sized by node count. Cross-group relationships shown as
  **aggregated edges** (e.g. "Research →3 references→ Pipeline").
- **Expand on click:** group → its folders/subgroups → files. Force simulation
  runs only over the **expanded frontier** (visible nodes), not all 393.
- **Containment is NOT drawn as edges** — it's expressed by nesting/expansion and
  group color. Only cross-cutting edges are drawn; collapsed groups aggregate
  their crossing edges to the group node.
- **Breadcrumb + "collapse all"** for orientation.
- Chosen over (A) group-anchored force clustering (still blobby) and (C) zoomable
  treemap/pack (great for hierarchy, weak for cross-links).

### B. Navigation primitives

- **Focus mode:** click a node → highlight its real relationships (callers, deps,
  tests, references), dim everything else; keep the existing detail panel.
- **Search:** expands/jumps to matching nodes (auto-expanding their ancestors).
- **Filters:** by group, by type, by issue.
- **Legend** that matches the (now-correct) node types.

### C. Data-model fixes (`extract_graph.py`)

1. **Fix `classify_path`** — normalize the path so its substring rules fire
   (the rules were written for leading-slash absolute paths but receive relative
   paths). Manuscripts, docs, slide decks (`.qmd`), notebooks, and lab-notebooks
   get correct types → meaningful colors + legend. Eliminates the 236-node
   `resource` blob. Add a regression assertion (no group should be >40% `resource`).
2. **Connect the research half** — parse markdown→markdown links across `docs/`
   and `research/` into `references` edges (target resolved to the linked file's
   node id; drop links whose target isn't a node, per the existing
   dangling-edge guard). Turns ~157 links into real edges.
3. **Light issue tagging** — stamp each node under an `issue_NNN/` folder with an
   `issue` attribute (e.g. `"issue": 218`) for color-by-issue / filter. An
   attribute, **not** a structural spine (purpose is atlas, not issue-driven).

### D. Out of scope (YAGNI)

- The 80-edge script-import call graph (pipeline is already well-connected).
- Any issue-tracker behaviour (that was the rejected "issue-driven" option).
- Server/runtime deps — stays one self-contained, offline HTML.

## Components & boundaries

| Unit | Responsibility | Output contract |
|---|---|---|
| `extract_graph.py` | Repo → graph data. Correct classification, real cross-edges, issue tags, dangling-edge/dup guards. | `graph.json`: `{nodes:[{id,label,type,group,issue?,…}], edges:[{source,target,type}]}` — every edge endpoint is a node id. |
| `project_map.html` | Render + interaction: collapsible hierarchy, aggregated edges, focus/search/filter, legend. | Reads `<script id="graph-data">`; no network. |
| `build_html.py` | Inline `graph.json` + vendored D3 → self-contained HTML. Idempotent. | unchanged (already done). |

The two render concerns — **hierarchy/expansion state** and **cross-edge
aggregation** — are the new logic in the HTML; keep them as separate functions so
each can be reasoned about alone.

## Testing / verification

- **Data invariants** (extend the existing check): 0 dangling edges, 0 duplicate
  ids, `resource` < 40% of any group, `contained_in` no longer the dominant *drawn*
  edge class, research group has > 10 real edges, issue tags present on
  `issue_NNN/` nodes.
- **Headless render (offline, Playwright + system Chrome):** default collapsed
  state renders ≤ ~15 top-level nodes; expanding a group reveals children;
  focus-on-click dims non-neighbors; search expands to a match; **no console
  errors; no network requests**; all verified with network blocked.
- Idempotent rebuild (counts stay 1/1/1) — already covered.

## Success criteria

A newcomer (or future-you) opens `project_map.html` and within ~30s can: see the
5 groups, expand Pipeline to read the dataflow, click a script to see exactly what
calls/tests/feeds it, and find any file by search — without ever facing a 393-node
hairball.
