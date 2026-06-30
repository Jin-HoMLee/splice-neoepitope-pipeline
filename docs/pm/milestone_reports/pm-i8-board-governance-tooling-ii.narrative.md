<!-- Author-owned narrative for pm-i8 - Board-Governance Tooling II. Sections 3/4/5 only.
     The script regenerates the HTML from this file + fresh board data;
     it never overwrites this sidecar once it exists. -->

## Deliverables (Review layer)

> **Delivery tally: 9 delivered (all COMPLETED), 0 descoped.** Unlike pm-i7 (which folded 2 `NOT_PLANNED` closes into its raw headline), every pm-i8 issue shipped real work — and the fix that makes this split machine-honest ([Issue #851](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/851)) shipped *inside* pm-i8, so this is the first closure report whose "delivered / descoped" line is computed, not hand-audited.

pm-i8 was the direct successor to pm-i7 and kept its character: a **board-governance tooling** pass that converted governance *rules* into enforced *mechanisms* and settled the remaining left-side model questions. Highlights by theme (full inventory in the appendix):

- **The parent/epic Status model — settled and mechanized.** [Issue #776](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/776) decided **Pattern A2** (park parents in a dedicated `Epic` Status, read progress off GitHub's native sub-issue bar) — eliminating the parent-Status drift *class* rather than policing it. [Issue #498](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/498) (Developer) shipped the GitHub Action that mirrors parent project Status from sub-issue lifecycle events, automating the leaf→parent rollup.

- **Commitment & accounting mechanics.** [Issue #754](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/754) replaced the count-only replenishment gate with the **per-role floor-5 / cap-18** gate (`check_ready_queue.sh`); [Issue #837](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/837) fixed the per-role counter to handle multi-role committed Issues correctly (fix the *counter*, not the label — a Scientist-caught premise correction). Together they made the Ready-queue gate role-aware.

- **Closure-evidence integrity.** [Issue #851](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/851) closed the masking gap the pm-i7 report itself surfaced — `milestone_report.py` now splits closed issues by reason, so a `NOT_PLANNED` close can no longer inflate the delivered count.

- **Branch / arc / cleanup tooling.** [Issue #578](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/578) shipped `scripts/new_branch.sh`, the canonical branch-name helper (this report's own branch was cut with it); [Issue #759](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/759) surfaced the active-arc slate as a pinned board card (the buried `arc_taxonomy.tsv` made glanceable); [Issue #790](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/790) wired the Friday branch-cleanup beat to auto-fire in the Scientist/Developer routines (per-clone self-clean).

- **Search eval.** [Issue #294](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/294) evaluated the GA'd GitHub Issues natural-language search and recorded the query patterns worth adopting.

## Carried-forward & routing

> **Note on the "Carried forward: 0" card.** The at-a-glance card counts only issues still *open and attached* to pm-i8 (`n_carried_forward = total − closed`); the three leaves below were *detached* (un-milestoned to Backlog on 2026-06-29), so the machine count is `0`. "Carried forward" in this section means *commitment* carried to a future iteration — the human sense the card doesn't capture. pm-i8 is the first report where a deliberate carve-forward makes the two senses diverge.

- **Three leaves decommitted to Backlog (2026-06-29), not failed delivery** — pm-i8 was emptied by a deliberate carve-forward when its committable work drained and no active WIP remained:
  - [Issue #864](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/864) (`awaiting-bot-review` auto-poll skill) — build was user-deferred; re-commits at a future `Backlog → Ready`.
  - [Issue #769](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/769) (extract board-governance conventions out of CLAUDE.md) and [Issue #745](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/745) (commitment-lag sweep) — both `arc:board-governance` at `arc-phase:later`; they re-commit when board-governance bumps back to `active`.
- **Routing decision: (c) extend the workstream — but the successor *milestone home* is itself under review.** Board-governance is a live, never-closing arc, so the work continues. However, **no `pm-i9` successor was opened**: the role-meta milestone-home question (does committing role-meta PM work even require a milestone?) is now an open governance item — [Issue #902](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/902) (filed 2026-06-30), cross-linked to [Issue #693](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/693). So the three carried-forward leaves sit **un-milestoned in Backlog** pending that decision, rather than being force-fitted into a new bucket.

## Retrospective (process/health)

- **Clean delivery, and the integrity fix is self-referential.** 9/9 COMPLETED, 0 descoped — and the one issue that makes that statement *trustworthy* ([Issue #851](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/851)) shipped in this same milestone. pm-i7 surfaced the masking bug; pm-i8 fixed it; this report is the first to benefit. The mechanism-over-memory loop closed on itself.
- **Character: mechanism-over-memory, again.** The majority of the 9 promoted a governance *rule* into an enforced *mechanism* — a replenishment gate (#754), a counter fix (#837), a branch helper (#578), a parent-Status Action (#498), an arc-slate card (#759), and the close-reason split (#851). The rung-3 escalation path discharged in a concentrated batch, continuing pm-i7's theme.
- **A mechanism shipped here already needs calibration — honestly.** The floor-5 gate (#754) drove **two over-commit incidents within days** of shipping (2026-06-29 i6-S3 to 8.5d; 2026-06-30 a false all-roles `[REPLENISH]` on a quiet board), spawning the revisit [Issue #902](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/902). This isn't a failed deliverable — it's evidence that **intuited parameters (floor = 5) need flow-data calibration**, the lesson web-confirmed against lean-queue best practice. A milestone that *builds* a gate closing the same week the gate's behavior gets a revisit Issue is the system learning at speed, not thrashing.
- **Cross-role review caught a premise error.** #837's framing ("multi-role is a violation") was corrected to "multi-role is legitimate — fix the counter, not the label" by the Scientist before it cemented — the establishing case for the `feedback_verify_premise_before_mechanizing` memory. Cross-role review working as the design intends.
- **Cycle time: avg 14.3d / median 7.2d** — median materially healthier than pm-i7's ~11.6d, dominated by Backlog dwell before commitment (late-commitment working as intended, not slow execution).
- **Cross-role split (delivered): 8 PM / 1 Developer (#498) / 1 Memory Manager (#759, dual-role).** PM specced the governance; the Developer carried the sub-issue→parent Status Action.
