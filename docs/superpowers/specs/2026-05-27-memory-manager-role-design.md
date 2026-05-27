# Memory Manager — 4th Role Design

**Status:** Draft for review · 2026-05-27 · Author: PM (this session)

**Context:** Brainstormed in response to the asymmetry surfaced during [PR #520](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/520) (the rule-lift ship for [Issue #496](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/496)) — memory file edits live in the personas repo with no clear owner of the git lifecycle, and 5-7 open `role:pm` Issues are actually memory-curation work that dilutes PM's focus. Community research ([DevelopersIO global-memory pattern](https://dev.classmethod.jp/en/articles/claude-code-global-memory-with-git/), [Letta Context Repositories](https://www.letta.com/blog/context-repositories), [Anthropic subagent docs](https://code.claude.com/docs/en/sub-agents), [Anthropic Issues #42844 / #19903 / #1628 / #60151](https://github.com/anthropics/claude-code/issues/42844)) confirmed: cwd is locked at session start, so "skill invoked from PM session" is structurally non-viable; the right shape is a session whose cwd IS the personas repo.

---

## 1. Role profile

**Name:** Memory Manager (MM). Matches functional naming of PM/Sci/Dev.

**Purpose:** Own the meta-process workstream of the project — keep the personas-repo memory durable, consistent, and ergonomic for PM/Sci/Dev. Explicitly cross-cutting and meta-process; not a bioinformatics role.

**In scope:**

- **Personas-repo git lifecycle.** All `git add` / `commit` / `push` on `claude-personas-splice-neoepitope-pipeline/` flow through MM. Replaces the current "personas-repo is not your responsibility" rule.
- **Memory-curation Issues.** [Issue #248](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/248) (frontmatter audit metadata), [Issue #326](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/326) (periodic consolidation), [Issue #346](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/346) (audit script for Always-in-effect rules), [Issue #353](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/353) (Pattern Language re-description), plus the pending MEMORY.md slimming audit (`project_memory_md_slimming.md` in PM role memory).
- **Always-in-effect rule promotions / lifts.** The shape of [Issue #496](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/496) (Dev rule lifted to shared) — MM authors + commits.
- **Cross-rule consistency.** Orphan files, broken cross-refs, duplicate rules, stale `Caught YYYY-MM-DD` dates.
- **Memory tooling.** Audit scripts, hook fire-log threshold reviews, MEMORY.md slimming when shared/MEMORY.md exceeds the documented ~14-rule cliff.

**Out of scope:**

- Bioinformatics work (Sci/Dev domain).
- Project board / milestone / triage / closure rituals (PM keeps).
- Code, tests, manuscript, lab notebook entries for project work.
- Standup posting (no MM channel).
- Lab notebook entries (personas-repo commit log is MM's journal).

**Role-label semantics:** New label `role:memory_manager` on project repo + personas repo. Existing `role:pm` curation Issues get relabeled (PM stays bystander on the label optionally).

**Domain-bespoke justification (per `feedback_domain_bespoke_roles.md`):** MM's justification is **workload separation**, not domain analogy — measurably 5-7 queued Issues + a known git-lifecycle gap + ongoing rule-curation rhythm that conflates with real PM work when PM owns it. External framing: "Memory Manager is meta-process (cross-domain, would apply to any multi-role Claude Code setup with growing memory infrastructure); PM/Sci/Dev remain bioinformatics-bespoke."

---

## 2. Workspace + session model

**Cwd at session start.** `cd ~/dev/GitHub/Jin-HoMLee/claude-personas-splice-neoepitope-pipeline && claude --add-dir ~/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline`

The personas repo IS the working directory. Decisive constraint: Claude Code's cwd is locked at session start ([Anthropic Issues #42844, #19903, #1628, #60151](https://github.com/anthropics/claude-code/issues/42844)) — there is no `/cd` command, the Bash tool resets cwd between calls, so any approach involving "session in repo A doing primary work on repo B" requires permanent `git -C` + absolute paths. Opening the session in the personas repo dissolves this.

**Memory dir:** Auto-derived by Claude Code from personas-repo path's hash → `~/.claude/projects/<personas-hash>/memory/`. At init, a one-time symlink wires that hash dir to `memory_manager/` inside the personas repo. Same mechanism PM/Sci/Dev use.

**Shared memory access:** `memory_manager/shared/` symlinks to `../shared/`. MM reads cross-role rules transparently like other roles.

**Project-repo access via `--add-dir`:** v2.1.20+ flag ([ClaudeLog reference](https://claudelog.com/faqs/--add-dir/)) extends reachable filesystem + loads any CLAUDE.md from the added directory. MM gets read access to project-repo files without its own clone.

**Personas-repo `.claude/` setup (new):**

- `CLAUDE.md` **at personas-repo root** (per [Anthropic .claude docs](https://code.claude.com/docs/en/claude-directory) — root is standard convention; `.claude/CLAUDE.md` is also valid but non-standard). Documents repo structure, memory file conventions, role-mapping, `<role>/shared` symlink mechanism.
- `.claude/settings.json` — permissions for `git` + `gh` CLI; hooks (e.g. closing-keyword foot-gun detection like we hit on PR #520).
- `.claude/agents/` (optional) — sub-agents for parallel work (cross-rule-consistency-scan, etc.).

**Personas-repo git workflow:** Direct-to-main, no PR. MM stages + commits + surfaces diff for user confirmation + pushes. The shared "commit, push, merge — three separate steps" rule applies. PR overhead skipped because (a) single-author post-transition, (b) the GitHub Action bot isn't configured to subscribe to personas-repo events, so automated review triggers don't fire on personas-repo PRs — PR overhead is unrewarded, (c) direct commits keep cadence appropriate for bursty workload.

**Standup + lab notebook:** Skipped. Coordination is GitHub-Issue-mediated. Journaling is via commit messages.

**Session attribution:** `**Created by:** Memory Manager` on Issues / PR comments / commit trailers. Personas-repo commits use `Co-Authored-By: <session-model-id> <noreply@anthropic.com>` matching the session's actual model per project convention (e.g. `claude-opus-4-7`, `claude-sonnet-4-6` — whichever the running MM session is on).

**Cross-role visibility:** PM/Sci/Dev sessions read personas-repo memory files via their existing `<role>/shared` symlinks (no change). They also **write/edit** personas-repo files mid-session as naturally as today (the memory dir is symlinked into their session's filesystem). The boundary under MM is *git lifecycle*, not authoring — only MM commits + pushes.

---

## 3. Authority + handoff model

**The boundary is git lifecycle, not authoring.** Any role can edit personas-repo memory files mid-session (transparent through symlinks). Only MM commits + pushes.

| Activity | Who |
|---|---|
| Author/edit memory files | Any role, mid-session, as part of normal work |
| Commit + push personas-repo | Memory Manager exclusively |
| Queued curation work (slimming, dedup, audit scripts) | MM authors AND commits (no other role naturally engaged) |

### Flavor A — Session-incident edits (common case)

Active role hits a slip, edits personas-repo files inline (lift a rule to shared, fix a stale `Caught` date, add `<!-- src: -->` annotation). At session end:

- Records edit briefly in lab notebook entry under a "Memory edits — for MM to commit" bullet, OR
- Tags via project-repo Issue comment / quick note (low-ceremony — single sentence).

**Migration-window note:** Until Phase 3 lands the "Memory edits — for MM to commit" bullet convention in `shared/feedback_lab_notebook.md`, this Flavor A pattern is informal — sessions during migration use ad hoc handoff notes (any sentence in the lab notebook entry calling out uncommitted personas edits works).

Next MM session:

- Bare `git status` at morning-routine start (MM's cwd IS personas repo — no `-C` needed) surfaces uncommitted files.
- Reads active role's note for context.
- Reviews diff, commits as-is (clean) or surfaces inconsistencies for resolution (e.g. "PM added inline rule X but didn't update Sci's mirror — should we?").

### Flavor B — Queued curation work (MM-authored)

Sustained tasks where no role is naturally the author. Flow:

1. Any role files Issue in **personas repo** with `role:memory_manager` label + adds to project board #9 (cross-repo ProjectV2 supports this — `gh issue create --repo Jin-HoMLee/claude-personas-splice-neoepitope-pipeline --project "JH M Lee Lab" ...`). Skeletal body OK.
2. Project board surfaces it.
3. MM opens personas-repo session, authors + commits the change(s), closes the Issue with commit SHA reference.
4. If AC requires a project-repo lab notebook entry (e.g. cross-role memory promotion), MM files a project-repo PR with just the lab notebook entry — same pattern as PR #520 today.

### Updated `feedback_memory_escalation.md`

Current rule:
> **On "you forgot X":** Find the existing memory. If it's only behind a link, copy it inline here under Always in effect. If no memory exists, create one AND inline it here.

New form:
> **On "you forgot X":** Find the existing memory. If it's only behind a link, copy it inline here under Always in effect. If no memory exists, create one AND inline it here. The edit ships at the next Memory Manager session (active role edits; MM commits). For sustained curation needs, file a `role:memory_manager` Issue.

### Cross-repo coordination mechanics

- **Inbound to MM:** GitHub Issues with `role:memory_manager` label (in either repo, both on board #9) + uncommitted personas-repo state.
- **Outbound from MM:** Issue closures, personas-repo commit log, shared/MEMORY.md changes surface at next `/memory` load (existing self-documenting rule).

---

## 4. Migration plan

**Principle:** minimum viable rollout + 4-week validation gate. Fully reversible (text/config/relabels only).

### Phase 1 — Pre-flight (project-repo state)

The MM rollout epic is filed as a **native GitHub parent Issue** per `feedback_parent_sub_issues.md`: no branch / no PR / no Size on the parent itself; Status mirrors from subs; sub-issues inherit milestone + priority; size rolls up. Recent precedent: [pm-i5 parent Issue #480](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/480) (closed 2026-05-26 via summary comment after all 4 subs Done).

PM session (this one or follow-up) files:

- **Parent Issue** "feat(roles): introduce Memory Manager as 4th role" — pm-i6 milestone, P2, `role:pm` + `role:memory_manager`. **Skeletal body** — no AC of its own; just links to this design doc and lists the sub-issues once they exist (linked via REST API per `reference_github_sub_issues.md`).
- **Native sub-issues** (each carries own branch + PR + Size + AC; linked to parent via REST API):
  - Sub 1: **Design — write + commit Memory Manager design doc** (this current work; lands the spec to `docs/superpowers/specs/`)
  - Sub 2: Personas-repo CLAUDE.md + `.claude/` setup
  - Sub 3: `memory_manager/MEMORY.md` initial content
  - Sub 4: Shared memory updates (single PR; rules are interrelated)
  - Sub 5: Role MEMORY.md updates (PM/Sci/Dev acknowledgments + escalation rule rewrite)
  - Sub 6: Enable Issues + create `role:memory_manager` label on personas repo
  - Sub 7: Relabel existing curation Issues
  - Sub 8: Bootstrap first MM session + 4-week validation
  - Sub 9: Cross-repo bot blind-spot follow-up — "should `scripts/audit_and_merge.sh` gate project-repo merges on the existence of a corresponding personas-repo commit (MM-side completion check)?"
- Create label `role:memory_manager` on project repo (one command, no Issue).
- Parent closes after all subs Done, via summary comment per pm-i5 parent-close precedent.
- Decision: keep pm-i6 milestone name as-is for now. Defer rename to pm-i7 if validation passes.

### Phase 2 — Personas-repo structural setup (one-shot, via `git -C` from PM session)

Chicken-and-egg: MM can't bootstrap its own dir from within an MM session — personas-repo cwd doesn't have MM's setup yet. PM session creates structure:

- `mkdir memory_manager/`
- `ln -s ../shared memory_manager/shared`
- Initial `memory_manager/MEMORY.md` with frontmatter + Always-in-effect placeholder.
- Personas-repo root `CLAUDE.md` (new).
- Personas-repo `.claude/settings.json` (minimum permissions).

User runs one-time `~/.claude/projects/<personas-hash>/memory/` symlink wiring (same per-clone setup as PM/Sci/Dev needed). To obtain `<personas-hash>`: run `claude` once in the personas-repo cwd — the auto-created `~/.claude/projects/<hash>/` path will contain exactly one new hash directory; that's the personas-repo hash. Symlink that dir's `memory/` to the personas-repo's `memory_manager/`.

### Phase 3 — Shared memory updates (PM session, one PR)

- `shared/feedback_team_structure.md` — add MM subsection.
- `shared/feedback_domain_bespoke_roles.md` — add MM justification paragraph (workload separation, not domain analogy).
- `shared/feedback_memory_escalation.md` — rewrite per Section 3.
- `shared/feedback_lab_notebook.md` — (a) note MM exempt from lab notebook entries (commits ARE the journal); (b) add the "Memory edits — for MM to commit" bullet convention for active-role sessions that edited personas-repo files inline (per Section 3 Flavor A).
- `shared/MEMORY.md` — drop "personas-repo git state is not your responsibility" line; replace with "Personas-repo git lifecycle is Memory Manager's responsibility" pointing at the new `team_structure.md` section.

### Phase 4 — Role MEMORY.md + Issue relabels (PM session, second PR)

- `pm/MEMORY.md`, `scientist/MEMORY.md`, `developer/MEMORY.md` — minor updates: replace personas-repo no-touch reflex with "file `role:memory_manager` Issue for curation work."
- **Update each role's morning routine** (in `<role>/feedback_morning_routine.md`) to include a `git status` scan on the personas repo at session start — surfaces uncommitted personas-repo state for chat acknowledgment + flagging for MM. This is the structural mitigation for the "active roles edit memory mid-session, forget to flag for MM" risk in Section 5.
- Relabel [Issue #248](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/248), [Issue #326](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/326), [Issue #346](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/346), [Issue #353](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/353) + the to-be-filed MEMORY.md slimming audit Issue from `role:pm` → `role:memory_manager`.

### Phase 5 — First MM session bootstrap

User opens MM session: `cd ~/dev/.../claude-personas-splice-neoepitope-pipeline && claude --add-dir ~/dev/.../splice-neoepitope-pipeline`. MM:

- Validates setup (memory loads, shared symlink works, `--add-dir` reaches project repo).
- Fleshes out `memory_manager/MEMORY.md` initial Always-in-effect rules.
- Picks first concrete curation Issue (likely MEMORY.md slimming audit — pending longest + natural "introduce yourself by curating your own neighbors" task).
- First personas-repo direct commit + push.

### Phase 6 — Validation gate (4 weeks post-bootstrap)

PM morning routine at +4 weeks measures:

- `role:memory_manager` Issues opened vs closed; net rate.
- MM session count (user self-report).
- Personas-repo commit frequency — uncommitted state lingering longer than acceptable?
- Subjective: does PM session feel less burdened by memory-meta work?

**If validation passes:** keep MM. Optionally rename pm-i6 to drop the PM prefix at pm-i7 carve.

**If validation fails:** rollback — relabel Issues back to `role:pm`, retire MM session pattern, keep personas-repo structure as inert dead code. Lab notebook entry captures the trial + outcome.

**Estimated effort:** ~5-7 project-repo PRs over the rollout — Sub 1 (this design doc), Subs 2-3 (personas-repo structural setup, each carrying a lab notebook entry as the project-repo deliverable), Sub 4 (shared memory updates journal), Sub 5 (role MEMORY.md updates journal), Sub 7 (Issue relabels journal), Sub 9 (cross-repo bot blind-spot follow-up — Dev tier). The personas-repo file changes themselves land via direct commits from PM during Phases 2-4 (the chicken-and-egg structural prep before MM exists), then via direct commits from MM in Phase 5+. Sub 8 (bootstrap MM session) has no project-repo PR (MM is exempt from lab notebook entries per the new shared rule). Phase 6 is observation.

---

## 5. Open questions / risks / non-goals

### Open questions

1. **Project board membership for personas-repo Issues — resolved.** Personas repo's Issues join project board #9 directly. Requires (a) Issues enabled on `claude-personas-splice-neoepitope-pipeline` (one-time setting), (b) MM files Issues in personas repo for substantive curation work + adds to board #9 at creation. All cross-repo references use full URLs (already the convention per link+prefix+keyword rule) to disambiguate `#N` namespace collisions. Routine commits (single-rule edits flagged by an active role's lab notebook bullet) stay un-Issued.
2. **MM session morning routine.** Sketch in `memory_manager/feedback_morning_routine.md` during Phase 5 bootstrap. Likely shape: read role memory, check `gh issue list --label role:memory_manager --state open` (both repos), bare `git status` for uncommitted state.
3. **Concurrent MM + active-role edits to same shared/* file.** Mostly theoretical (MM is async). Mitigation: MM session begins with `git pull` before any edits. Active roles' edits are uncommitted local state — no git-level conflict until MM commits. Edge case: concurrent MM + PM sessions both editing same file → last-write-wins on disk. **Defer to bootstrap session for practical test.** Test outcomes: (a) **acceptable risk** — across N MM sessions (≥10) where active-role edits also occurred, no observed conflict → no further mitigation needed; (b) **mitigation required** — ≥1 conflict observed → propose a lock-file convention (touch `.mm_session_active` on session start, error on overlap) or a hard-rule "no MM session while another role session is active" added to MM Always-in-effect.
4. **Cross-repo bot blind-spot follow-up.** Re-scoped: not "make the bot see personas changes" (by-design boundary) but "should `scripts/audit_and_merge.sh` gate project-repo merges on the existence of a corresponding personas-repo commit (MM-side completion check)?" Folded into Phase 1 sub-issues (Sub 9).
5. **Directory naming `memory_manager/` vs `mm/`.** Default `memory_manager/` for consistency with `pm/`, `scientist/`, `developer/` long-form naming.

### Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| MM workload too sparse → role atrophies → uncommitted state accumulates | Medium | 4-week validation gate measures commit cadence; rollback to per-role commits if sparse |
| Active roles edit memory mid-session, forget to flag for MM, edits sit uncommitted | Medium-High | Every role's morning routine includes `git status` scan on personas repo; lab notebook entry checklist for non-routine sessions asks "any uncommitted personas edits?" |
| Latency for "you forgot X" rule lifts annoys active roles | Low-Medium | Section 3 hybrid: edits are immediate, only commits are deferred. Fix is immediate as today. |
| Personas-repo `.claude/settings.json` permissions drift from project-repo over time | Low | Keep minimal; document intentional divergence in personas-repo CLAUDE.md |
| Migration takes longer than 4 PRs | Medium | Phases independent; can land separately. MM can bootstrap with stale shared memory and update incrementally. |
| External readers interpret "4-role multi-agent" as overclaiming | Low | `feedback_multi_role_not_multi_agent.md` already addresses this; MM is a role in our multi-role workflow, not an autonomous agent |

### Non-goals

- MM does NOT become a 4th bioinformatics role. Justification is workload separation, not domain analogy.
- MM does NOT take over project-repo PR management — PM keeps.
- MM does NOT post to team_standup.md.
- MM does NOT replace per-role memory authoring — each role still writes own role-specific MEMORY.md content; MM curates + commits but doesn't dictate substance.
- No PRs on personas-repo. Direct-to-main; user-confirmed push.
- No automatic / cron-driven MM session — user-launched, same as other roles.

---

## Implementation handoff

Next step: invoke `superpowers:writing-plans` to translate this design into a step-by-step implementation plan with review checkpoints. Plan will sequence Phases 1-5 across discrete commits/PRs; Phase 6 is observation, not implementation.
