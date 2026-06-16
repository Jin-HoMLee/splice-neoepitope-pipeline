# Prose-Dependency Reconciler — Design (Issue #722)

**Status:** design approved (brainstorm), pending spec review → writing-plans
**Issue:** [#722 — populate native blockedBy on the dependency graph](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/722)
**Arc:** `arc:board-governance` (active)
**Author:** PM
**Date:** 2026-06-16

## 1. Problem

Cross-Issue dependencies on board 9 live in **prose** — Issue bodies say "depends on #M" / "blocked on #M" — but the corresponding **native GitHub `blockedBy` edge** is usually never wired. So `is:blocked` audits return ~0 despite real dependency chains (e.g. #707→#708→#709 calibration; #680→PR #714), and "what's blocked?" can't be answered mechanically.

The cost is concrete and recent: on 2026-06-16, [Issue #745](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/745) was committed to Ready by hand while its body said "depends on #722" (open) — a prose-only dependency with zero native edge, so `is:blocked` passed it clean. The same blind spot caused [Issue #594](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/594) to sit un-pulled for 4 days after its blockers cleared.

The convention and the mechanism already exist (`shared/feedback_dependency_tracking.md`: "set the native edge at creation/discovery", via the `addBlockedBy` GraphQL mutation or the simpler REST `dependencies/blocked_by` POST). What's missing is (a) a **backfill** of the existing prose deps into native edges, and (b) a way to **sustain** graph accuracy so the convention stays honest.

## 2. Scope (decision: B — backfill + reusable reconciler)

Three options were weighed:
- **A. Backfill only** (AC-as-written): one-time scan + wire + verify; documented convention sustains it. Rejected — leaves the DoR/best-next blocked-check re-grepping by hand every commit.
- **B. Backfill + a reusable prose↔native drift-checker** *(chosen)*. The backfill *requires* writing the scan anyway; making it reusable is cheap and gives the convention a drift-detector — the same pattern every other board convention uses to stay honest.
- **C. B + a create/edit hook** that flags a prose-dep-without-edge. Rejected as premature — the mechanism-over-memory ladder reserves hooks for rules that have slipped ≥2× on the same shape; this is slip #1.

**Out of scope:** cross-repo dependencies (single-repo project, untested); auto-removal of closed-blocker edges (the native graph auto-reflects closure — no manual cleanup); a CI gate (the `--check` mode is built, but wiring it into a workflow is a separate follow-up if wanted).

## 3. Deliverable

`scripts/pm/scan_prose_deps.py` — a **prose↔native dependency reconciler**. Plain `python3` + `gh` subprocess, alongside `scripts/pm/recheck_milestone.py` (same convention — no special env; `gh` on PATH). Run once to drive the backfill; kept afterward as a board-hygiene drift-checker and as the tool the DoR/best-next blocked-check calls instead of ad-hoc grep.

## 4. Architecture — four thin layers

Each layer is independently testable; the parse layer holds the risk and gets the real test coverage.

| Layer | Job | How |
|-------|-----|-----|
| **fetch** | get open issues + bodies | `gh issue list --state open --limit 1000 --json number,title,body,state`. Uses the **issue list**, not the project board — so the board's Done-first pagination trap does not apply (that trap is `projectV2.items`-specific). `--limit 1000` per the default-limit convention. |
| **parse** | body → `(dependent N, blocker M)` pairs | regex over a blocker-phrase allowlist; narrative-phrase exclude. Detailed in §5. |
| **reconcile** | diff prose pairs vs the native graph | for each pair: is M open? does a native `blockedBy N→M` edge already exist? → classify into `needs-wiring` / `already-wired` / `closed-blocker (skip)` / `un-wireable (PR blocker)`. |
| **act** | report / wire / exit-code | mode-dependent (§6). |

## 5. Parse layer — false-positive control

The parse layer is the only risky part: a false positive that gets wired creates a **wrong `blockedBy` edge**, which marks real work un-pullable (the inverse failure of the bug we're fixing). Controls:

- **Blocker-phrase allowlist** (case-insensitive), each immediately followed by `#<digits>`:
  `depends on`, `blocked by`, `blocked on`, `blocked-by:`, `gated on`, `requires`.
- **Narrative-phrase exclude** (these are *why*-context, NOT *what*-blockers — wiring them is wrong):
  `informs`, `consumes`, `relates to`, `related to`, `follow-up to`, `supersedes`, `superseded by`, `see`, `cf.`.
- **Reconcile-stage filters** (not parse-stage — keep parse purely syntactic):
  - **Closed blocker M → drop** (convention: only edges for currently-open blockers).
  - **M is a PR, not an issue → report as `un-wireable (prose-only)`**, never auto-acted. Native dependencies are issue↔issue; a "blocked by PR #714" stays prose. Detect via the issue/PR type on lookup.
  - **Self-reference (N == M) → drop.**
  - **Duplicate pairs → dedupe.**
- **Never auto-wire blind.** Default mode only reports. Wiring requires a human review gate (§7).

**Heuristic, not perfect.** The allowlist will miss exotic phrasings and may over-catch a blocker phrase used in historical/negated context ("previously blocked by #X, now resolved"). Both failure modes are caught by the human review gate before `--apply`. The scanner's job is to *surface candidates*, not to decide autonomously.

## 6. Modes (one tool, four entry points)

- `--report` *(default)* — print the drift table: `dependent | blocker | blocker-state | native-edge? | action`. Read-only. Drives the backfill review and the hygiene sweep.
- `--apply <issue-list>` — wire native edges for the **PM-reviewed subset** (passed explicitly, drawn from the `--report` output — the report *is* the review gate), via the REST `dependencies/blocked_by` POST. Refuses to wire anything outside the `needs-wiring` classification (closed-blocker, PR-blocker, already-wired are no-ops). Prints each edge as created. Idempotent.
- `--check` — exit non-zero if any drift exists (open prose blocker with no edge). For the board-hygiene sweep / future CI.
- `--issue N` — single-issue scope; composes with the above. This is what the **DoR commitment-act check** and the **best-next backstop** call: "does candidate N have an unwired open prose blocker?"

## 7. Backfill execution flow (the AC)

1. `scan_prose_deps.py --report` → full drift table across open issues.
2. **PM reviews** — strike false positives, PR-blockers, and any narrative the allowlist over-caught. Human gate; non-negotiable per the blind-wire risk (§5).
3. `scan_prose_deps.py --apply` on the reviewed set → native edges wired.
4. **Verify:** `gh search issues "is:blocked"` returns the real blocked set (AC item 2). Spot-check a known chain (#707→#708→#709).

## 8. Integration — the "sustain" half

- **DoR / commitment-act** (`shared/feedback_board_hygiene.md`) — the blocked-check added 2026-06-16 calls `--issue N` instead of ad-hoc grep.
- **Best-next backstop** (`shared/feedback_best_next_issue.md` Step 1) — same `--issue N` call.
- **Board-hygiene sweep** — periodic `--check` surfaces drift (a prose dep authored without an edge) as a hygiene finding — the drift-detector that keeps "set edge at creation" honest.

These three call-sites are *documented* now; whether they invoke the script literally vs. remain a convention pointer is a plan-level choice (the script makes either cheap).

## 9. Testing

- **pytest on the parse layer** (`workflow/tests/.venv`) — fixture bodies → expected pairs: allowlist hits, narrative excludes, multi-blocker, closed-blocker filtering at reconcile, PR-blocker flagging, self-reference, dedupe. The false-positive risk lives here, so this is where coverage concentrates.
- **Light integration** on reconcile/act — mock the `gh` calls (or a dry-run against a couple of known live chains). No live mutations in tests.

## 10. Doc / convention updates

- Cross-ref the scanner into `shared/feedback_dependency_tracking.md` "How to apply" (the scanner is *how* you find/reconcile prose deps) and name the three call-sites. The convention ("stated → edge at creation") is unchanged; #722 gives it a drift-checker, not a new rule.

## 11. Acceptance criteria (mapped from the Issue)

- [ ] Existing open prose dependency chains reflected as native `blockedBy` edges (§7 steps 1-3).
- [ ] `is:blocked` / blockedBy audit returns the real blocked set (§7 step 4).
- [ ] Convention documented — scanner cross-ref'd into `feedback_dependency_tracking.md`; the three call-sites named (§10).
- [ ] `scripts/pm/scan_prose_deps.py` lands with the four modes (§6) and parse-layer pytest coverage (§9).

## 12. Open questions for the plan

1. **Single-issue `--issue N` perf** — does the DoR/best-next call need its own native-edge lookup per candidate, or can it reuse a cached full-scan? (Likely per-candidate; the call is rare and interactive.)
2. **PR-blocker representation** — leave un-wireable PR blockers as prose-only forever, or open a follow-up to track them another way? (Lean: prose-only is fine; PR blockers resolve fast.)
