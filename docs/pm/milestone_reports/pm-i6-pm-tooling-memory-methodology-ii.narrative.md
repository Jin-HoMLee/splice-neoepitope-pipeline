<!-- Author-owned narrative for pm-i6 - PM Tooling, Memory & Methodology II. Sections 3/4/5 only.
     The script regenerates the HTML from this file + fresh board data;
     it never overwrites this sidecar once it exists. -->

## Deliverables (Review layer)

<!-- Lead role: what shipped, grouped by deliverable. Auto-seed below. -->

- #696 Project map atlas — interactive whole-project D3 visualization (tools/project_map/)
- #642 feat(pm): thread issue timestamps into board_open_items.py (+ dormancy --stale-days sweep)
- #633 research(pm): board governance review — parent-as-milestone-anchor, Milestones vs Iteration field, deferral tracking, WIP limits
- #623 feat(roles): Sub 3+4 of #527 — personas-repo MM workspace (CLAUDE.md + .claude/ + memory_manager/)
- #618 chore(pm): dock — proves-out review for recheck_milestone hook (PM-local standing + capacity right-sizing)
- #587 Phase 1: late-commitment conventions rewrite (5 files)
- #584 docs: CLAUDE.md board-governance section — late-commitment model
- #583 chore: retime residual milestone-at-triage refs (best_next_issue / dependency_tracking / shared MEMORY)
- #582 feat: late-commitment morning-routine split + Ready-queue replenishment nudge
- #581 design: sub-issue milestone inheritance under late commitment
- #580 design: board status governance — migrate the left side to late-commitment Kanban (phased)
- #577 fix(ci): TestLiveIntegrationSmoke docstring says "skipped by default" but live tests run in ci-tools-pytest
- #569 migrate team coordination off team_standup.md → GitHub board + Discussions (retire the file)
- #567 design: Sub 5-7 of #527 — personas-repo governance (edit boundaries + MM-owned commit lifecycle)
- #561 fix(hooks): post_gh_pr_create misses `gh pr create` after a VAR=value prefix
- #550 infra: PostToolUse hook on gh pr create — automate project-add + Ready-for-review Status flip
- #549 infra: PreToolUse hook on gh issue develop — refuse parent-Issue targets
- #548 pm: session arc tracker 2026-05-28 — memory-slim epic carve lab notebook PR
- #530 feat(roles): Sub 2 of #527 — write + commit MM rollout implementation plan
- #528 design: Sub 1 of #527 — write + commit Memory Manager design doc
- #454 chore(pm): promote target-date sync hook from settings.local.json -> settings.json

### Auto-seed (lab-notebook + closing comments)

- **2026-06-15** (developer): — closure-audit reason-awareness + two Issue-lifecycle guardrails ([PR #744](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/744) closes [Issue #743](https://github.…
- **2026-06-12** (developer): — STAR sensitivity-flag benchmark sweep shipped ([PR #720](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/720) closes [Issue #411](https://github.com/Jin-HoMLee/spl…
- **2026-06-11** (developer): — Flatten `predictions/` wrapper dir merged: tool-named sibling folders + doc-sync ([PR #650](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/650) closes [Issue #435…
- **2026-06-11** (developer): — STAR annotated-flag cross-check merged: conflict resolution + merge-seam test ([PR #649](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/649) closes [Issue #375](h…
- **2026-06-11** (developer): — Cloud re-run robustness merged: pre-run snapshot + stale-log clear ([PR #666](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/666) closes [Issue #658](https://gith…
- **2026-06-11** (developer): — GTEx pan-tissue filter merged: sole-filter unit test + sign-off ([PR #653](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/653) closes [Issue #212](https://github.…
- **2026-06-10** (developer): — DataCite fallback for arXiv DOIs in zotero_add.py ([PR #648](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/648) closes [Issue #641](https://github.com/Jin-HoMLee…
- **2026-06-04** (developer): — Cloud re-run robustness + CI conda env-solve guard ([PR #666](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/666) closes [Issue #658](https://github.com/Jin-HoMLe…
- **2026-06-04** (developer): — STAR env re-broke overnight from upstream bioconda churn; pinned DOWN to 2.7.10b ([PR #661](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/661) closes [Issue #629…
- **2026-06-03** (developer): — Integrate the GTEx pan-tissue filter into filter_junctions ([PR #653](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/653) closes [Issue #212](https://github.com/J…
- **2026-06-03** (developer): — GTEx pan-tissue novel-junction blacklist (Snaptron gtexv2) ([PR #598](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/598) closes [Issue #211](https://github.com/J…
- **2026-06-03** (developer): — star.yaml conda prerelease-ordering trap ([PR #645](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/645) closes [Issue #629](https://github.com/Jin-HoMLee/splice-n…
- **2026-06-03** (developer): — dependency release-feed poller + pure-delta version watermark ([PR #640](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/640) closes [Issue #639](https://github.co…
- **2026-06-01** (developer): — REPO passthrough for the pre-merge lab-notebook gate ([PR #615](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/615) closes [Issue #607](https://github.com/Jin-HoM…
- **2026-05-31** (developer): — bot-review-offer pre-merge gate ([PR #600](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/600) closes [Issue #443](https://github.com/Jin-HoMLee/splice-neoepitope…
- **2026-05-31** (developer): — TestLiveIntegrationSmoke skip contract ([PR #589](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/589) closes [Issue #577](https://github.com/Jin-HoMLee/splice-neo…
- **2026-05-31** (developer): — surface pre-cap presenter total for the strong-presenter truncation notice ([PR #579](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/579) closes [Issue #226](http…
- **2026-05-30** (developer): — post_gh_pr_create matches `gh pr create` after a `VAR=value` prefix ([PR #575](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/575) closes [Issue #561](https://git…
- **2026-05-29** (developer): 21:30 UTC — Editor: Developer
- **2026-05-28** (developer): 20:38 UTC — Editor: Developer
- **2026-06-15** (pm): 17:53 UTC — Editor: PM
- **2026-06-13** (pm): 21:11 UTC — Editor: PM
- **2026-06-12** (pm): 18:35 UTC — Editor: PM
- **2026-06-11** (pm): 20:40 UTC — Editor: PM
- **2026-06-10** (pm): 18:21 UTC — Editor: PM
- **2026-06-05** (pm): 15:14 UTC — Editor: PM
- **2026-06-04** (pm): 20:31 UTC — Editor: PM
- **2026-06-03** (pm): 13:42 UTC — Editor: PM
- **2026-06-01** (pm): 15:36 UTC — Editor: PM
- **2026-05-31** (pm): 23:45 UTC — Editor: PM
- **2026-05-30** (pm): 15:20 UTC — Editor: PM
- **2026-05-29** (pm): 22:16 UTC — Editor: PM
- **2026-05-28** (pm): 19:38 UTC — Editor: PM

## Carried-forward & routing

<!-- PM: issues that didn't close + where they went (carve / arc) + the
     closure-routing decision (a/b/c/d). -->

- _(none carried forward)_

## Retrospective (process/health)

<!-- PM: was it healthy? what to improve? WIP/aging observations. -->

