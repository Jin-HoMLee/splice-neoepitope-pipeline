# Milestone Closure Report — Design Spec

- **Issue:** [Issue #752](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/752)
- **Status:** Design — approved 2026-06-16; pending implementation plan
- **Author:** PM
- **Arc:** `arc:board-governance`

## 1. Purpose

A per-milestone, self-contained **HTML closure report** produced *before* a milestone closes. It serves three audiences in layered sections — **PM retrospective**, **lab seminar**, **portfolio showcase** — from one artifact.

It formalizes our existing milestone-closure-routing convention (PM triggers a collaborative routing decision; Sci/Dev contribute domain expertise) into a concrete, reviewable artifact, authored via the **author-editor-critic** triad rather than written solely by PM.

### Why this shape (grounding)

- **Multi-agent best practice — author-editor-critic / generation-reflection.** The dominant collaborative-document pattern separates an *author* (proposes), a *critic* (flags flaws), and an *editor* (adjudicates, finalizes); author ≠ editor, and a distinct critic raises quality. Maps onto roles we already have.
- **Agile best practice — review ≠ retrospective.** The *sprint review* covers product/deliverables (stakeholder-facing); the *retrospective* covers process/health (what to improve). The whole team contributes; the facilitator is not the sole author. The report carries **both** layers.

| Triad role | Our mapping |
|---|---|
| Author (deliverables narrative) | the milestone's **lead role** (Scientist / Developer) |
| Editor (adjudicate, finalize) | **PM** |
| Critic (flag flaws) | **`@claude review`** on the PR (existing infra) |

## 2. Authoring flow (Report-as-PR, aggregate-assisted)

1. **Scaffold.** PM (at the milestone-health beat, when a milestone is ready to close) runs `scripts/pm/milestone_report.py <milestone>`. It emits the HTML with **data sections filled** and a **narrative first-draft** auto-seeded from the lead role's lab-notebook entries + issue closing comments for that milestone, into an author-owned markdown sidecar.
2. **Author.** The milestone's **lead role** writes/refines the Deliverables narrative on the branch.
3. **Edit.** **PM** edits, writes the Carried-forward & routing section + the Retrospective layer, and requests **`@claude review`** (critic).
4. **Finalize.** Address review → merge PR → close the milestone (the close is still the PM-coordinated routing act; the report is its evidence).

The flow is a **closure-ritual convention, not a hard gate** (see §7).

## 3. Architecture

One script plus a template and a per-milestone markdown sidecar. Follows the existing `workflow/scripts/generate_report.py` (Jinja → HTML) and `tools/project_map/build_html.py` (self-contained, vendored assets) precedents.

```
scripts/pm/milestone_report.py        # orchestrator + the four layers below
scripts/pm/templates/milestone_report.html.j2   # Jinja2 template (inline CSS)
docs/pm/milestone_reports/<slug>.html            # generated artifact (committed)
docs/pm/milestone_reports/<slug>.narrative.md    # author-owned sidecar (committed)
```

### Layers (within `milestone_report.py`)

1. **Data layer** — reuse the `scripts/board_open_items.py` paginated GraphQL helpers. Pull every issue in the milestone (open **and** closed): number, title, status, `role:*`, size, priority, `arc:*`, linked PRs (`closingIssuesReferences`), `createdAt`, `closedAt`, `subIssuesSummary`. Source of truth = board #9 + `gh issue list --milestone` + per-issue `gh issue view --json`.
2. **Metrics layer** — pure, unit-testable functions computing the headline numbers (§5). No I/O; takes the data-layer structs, returns numbers.
3. **Aggregation layer** — read `research/lab_notebook/<role>.md` entries dated within the milestone window + issue closing comments; assemble a first-draft narrative markdown. Best-effort: if nothing is found, emit the stub with a `<!-- no auto-seed found; author from scratch -->` marker.
4. **Render layer** — Jinja2 renders the data sections + injects the sidecar markdown (converted to HTML) into one self-contained file (inline CSS; vendor any JS, no CDN).

### Sidecar ownership

The `<slug>.narrative.md` sidecar holds **only** author-owned prose (Deliverables, plus PM's Carried-forward/Retrospective). The script **never overwrites** an existing sidecar — on re-run it regenerates the HTML from the (possibly hand-edited) sidecar + fresh board data. First run creates the sidecar from the auto-seed; subsequent runs preserve author edits. This keeps generated data and human narrative cleanly separated.

## 4. Report sections (content + ownership)

| # | Section | Content | Owner |
|---|---|---|---|
| 1 | **Header** | milestone name, stage S#, arc(s), opened→closed dates, duration | script |
| 2 | **At-a-glance** | N closed / M carried-forward, throughput, avg cycle time, per-role count bar | script |
| 3 | **Deliverables** (Review layer) | what shipped, grouped by deliverable/issue, with PR + lab-notebook + slide links | **lead role** (script seeds from lab-notebook + closing comments) |
| 4 | **Carried-forward & routing** | issues that didn't close + where they went (carve / arc) + the closure-routing decision (a/b/c/d) | **PM** (script seeds the carried-forward list) |
| 5 | **Retrospective** (process/health) | was it healthy, what to improve, WIP/aging observations | **PM** |
| 6 | **Inventory appendix** | full table of all milestone issues (status, role, size, PR, dates) | script |

Sections 1/2/6 are fully script-generated; 3/4/5 live in the sidecar.

## 5. Metrics (headline only)

Computed as pure functions over the data layer:

- **Duration** — earliest issue `createdAt` (or milestone creation) → milestone `closedAt`, in days.
- **Throughput** — closed issues per week over the duration.
- **Avg cycle time** — mean of (`closedAt − createdAt`) across closed issues, in days; also report median (robust to outliers).
- **Per-role counts** — closed-issue counts grouped by `role:*`.

Explicitly **excluded** (YAGNI): cumulative-flow / burn-down charts, cycle-time distributions, status-transition timing (board has no native transition timestamps; approximating is high-effort, low-value here).

## 6. Data sources

- Board #9 via `gh api graphql` (paginated; reuse `board_open_items.py` cursor loop + field helpers).
- `gh issue list --milestone "<full name>" --state all --json …`.
- Per-issue `gh issue view N --json comments,closingIssuesReferences,createdAt,closedAt,labels,…`.
- `research/lab_notebook/<role>.md` for the narrative auto-seed.

## 7. Process integration

- **Convention, not a hard gate.** Documented as a step in the morning-routine milestone-health beat (1c) — when a milestone is proposed for close, generate the report as the close-evidence step — plus a `feedback_*` memory rule. PM-coordinated; no machine block on milestone close (closing milestones is infrequent; a gate would be over-machinery for the frequency).
- **Linkage.** The merged report is linked from the milestone-close comment and the PM lab-notebook entry.

## 8. Where it lives / publishing

- Committed to the repo at `docs/pm/milestone_reports/<slug>.{html,narrative.md}`.
- `<slug>` = a sanitized milestone name (e.g. `pm-i6-pm-tooling-memory-methodology-ii`).
- GitHub Pages for the showcase audience is an **optional future** enhancement, not in this scope.

## 9. Testing

- **Unit tests** for the metrics pure-functions (duration, throughput, avg/median cycle time, per-role counts) with small synthetic fixtures — collected by the existing tools-pytest CI job.
- **Pilot dry-run** against an already-closed milestone — **pm-i6** (closed 2026-06-16, 21 closed issues) is a ready real-world fixture. Verify the HTML renders, metrics are sane, and the auto-seed pulls real lab-notebook/closing-comment material.

## 10. Scope boundaries

**In scope:** the script + template, the four layers, the six sections, headline metrics, the sidecar convention, the closure-ritual memory rule, unit tests, the pm-i6 pilot.

**Out of scope (YAGNI):** flow charts, a hard close-gate, a separate seminar deck (the experiment/eval/decision deck conventions already cover talks), a live-hosted dashboard, GitHub Pages publishing, status-transition timing metrics.

## 11. Open questions for the implementation plan

- Exact slug-sanitization rule for milestone names → report filename.
- Whether the narrative auto-seed converts markdown→HTML with an existing dependency already in the env (e.g. `markdown`/`jinja2`) or a minimal inline converter, to avoid adding a new dep to the `snakemake` env.
- Markdown→HTML for the sidecar: which env runs the script (likely `snakemake` conda or a small `scripts/pm` venv) — pin in the plan.
