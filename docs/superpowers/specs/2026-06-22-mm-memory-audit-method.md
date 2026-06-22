# MM Recurring Memory-Audit — Method Spec

**Status:** Operational · 2026-06-22 · Author: MM (this session)

**Context:** A recurring audit of the personas memory repo (`claude-personas-splice-neoepitope-pipeline`), migrated to the Memory Manager role from cerebrum's retired weekly `RemoteTrigger` ([Issue #20](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/20); cerebrum PR #33). Run 1 shipped 2026-06-22 via personas [PR #66](https://github.com/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/pull/66). This is the durable home for the method + the rolling baseline-watermark, so neither is lost when each run's tracking issue closes. Sibling to [2026-05-27-memory-manager-role-design.md](2026-05-27-memory-manager-role-design.md).

---

## 1. Purpose + scope

Keep the personas memory corpus internally consistent and lean: catch inline-vs-source **drift**, classify **new files**, scan for **stale path-refs**, **dedup** near-duplicate always-loaded rules, and apply the safe fixes — without silently deleting enforcement.

**Scope:** operative memory files — `{developer,pm,scientist,memory_manager}/MEMORY.md`, `shared/MEMORY.md`, `shared/feedback_*.md` + role `feedback_*.md`/`reference_*.md`, `CLAUDE.md`. **Out of scope:** `episodes/`, `_retired/`, `.silent/`, `shared/drafts/` (immutable / journal / archive).

## 2. Cadence + recurrence model

- **Queue-driven** (not scheduled). Each run is pulled from the `role:memory_manager` queue during a normal MM session. Rationale: an unattended scheduled run can't satisfy the personas push-gate (push only on the user's OK, or under the standing memory-push authorization with the diff surfaced) — so the audit runs attended. (A scheduled `JasonEtco/create-an-issue` cron is the alternative if that ever changes.)
- **One issue per run, not a standing tracker.** Each run = a fresh `role:memory_manager` + `audit` Issue, closed by that run's PR; its findings live in that PR. Do **not** reopen-and-reuse one issue across runs — that destroys per-run history. Full convention: personas `shared/feedback_recurring_issue_tracking.md` (GitHub's own scheduled-issue pattern; web-validated 2026-06-22).
- **The rolling baseline-watermark lives here** (§3), not in a closed issue.

## 3. Baseline rule

- **First run:** baseline at the merge-base commit `9173558` — audit everything landed since (a run-*date* would skip already-landed commits).
- **Later runs:** baseline at the **prior run's coverage commit** (the merge commit of the prior run's PR).

**Watermark log:**

| Run | Date | Coverage | Next baseline |
|---|---|---|---|
| 1 | 2026-06-22 | `9173558..57676bd` (PR #66, merge `420a2a6`) | **`420a2a6`** |

Update this table at the end of each run; the bottom "Next baseline" is where the next run starts.

## 4. The audit method

1. **Drift check.** For every `<!-- src: … -->` annotation in the "Always in effect" sections of the role `MEMORY.md` files, diff the inline rule against its cited source. Flag `source-newer` / `inline-newer` / `contradiction` / `broken-ref`. Faithful compression (terser but behaviorally identical) is **not** drift. This is the highest-value check — an always-loaded rule loads into every session.
2. **New-file classification.** Files added/changed since baseline → keep / dedup / move / stale.
3. **Stale path-ref scan.** Grep for `.claude/memory/…` or `memory/<role>/…` prefixes (the convention is bare file-relative). Distinguish real drift from legitimate mentions (the convention rule *describing* the forbidden form, branch-name examples, cerebrum's actual path, immutable drafts).
4. **Dedup.** Near-duplicate rules across role `MEMORY.md` + `shared/`. A role inline that merely restates a rule **already in `shared/MEMORY.md` Always-in-effect** is real per-session cliff cost (remove it; the role gets it from shared). A role inline that carries genuine role-specific nuance, or whose content is *not* in shared's always-loaded set, is **not** a safe dedup — keep it.
5. **Apply safe fixes** → the personas PR flow (§5).

**Method note (run-1 calibration):** "dup of shared L*N*" means the *directive* is the dup — not that every word is in shared's always-loaded text. An elaboration that lives only in a reference file is demoted, not deduped; weigh that before removing.

## 5. Applying fixes — safety + flow

- **MM-safe levers:** *compress* (strip incident-history tails, keep directive + Why + pointer) and *delete* (a genuinely dead rule). Pure personas-repo work, zero enforcement risk.
- **NOT MM-only:** *hook/skill demotion* moves enforcement into the code repo (`.claude/`). **Never strip a rule's `MEMORY.md` text before its hook/skill replacement is live in the code repo** — that silently deletes enforcement. Carve the hook/skill bundle to a Developer sub-issue; don't do it from an MM audit.
- **Bullet-removal hygiene:** after deleting an always-loaded bullet, re-check newlines (`grep -nP '\S-->- \*\*'` / read the surrounds) — a swallowed newline silently merges the next bullet (invisible in a unified diff).
- **Flow:** the audit Issue is board-tracked → branch (`gh issue develop`) in a worktree → PR → review → stop at the merge gate (the human merges). Verify each agent/tool finding with your own `Read` before editing. A `Workflow` fan-out (per-role drift + adversarial verify + survey) scales the read-heavy phases.

## 6. Worked example — run 1 (2026-06-22)

13-agent `Workflow` over a ~1-month backlog (143 changed operative files). Found **4 drift fixes** (1 contradiction, 2 broken-refs, 1 skip-list/carve-out) + **3 dedups** (each a verbatim restatement of a shared Always-in-effect rule), −18 net always-loaded lines. Adversarial verify killed 1 false positive (a "TodoWrite drift" — the source deliberately retained the term + carried an in-place translation note). Two recommended dedups were **held by judgment** (a shared-promotion that would load an MM-irrelevant rule; a demotion that wasn't actually shared-covered). Carryforwards D2 + NF1 verified moot. Path-ref scan clean. Shipped via [PR #66](https://github.com/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline/pull/66) (reviewed: approve + one calibration note, addressed).

## 7. Open items

- **Automation:** if attended cadence ever relaxes, a scheduled `create-an-issue` cron can open each run's issue (`shared/feedback_recurring_issue_tracking.md`).
- **Watermark durability:** this table is the watermark home today; if a tooling CLI for the personas store lands, the watermark could move into it.

## Implementation handoff

Operational now — run from the `role:memory_manager` queue. Trail: [Issue #527](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/527) (MM rollout) · [Issue #20](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/20) (audit migration) · personas `shared/feedback_recurring_issue_tracking.md` (recurrence convention).
