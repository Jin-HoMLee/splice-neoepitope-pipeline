# `.claude/` → `.agents/` canonical config-dir migration - design

- **Issue:** [#861](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/861)
- **Date:** 2026-06-25
- **Author:** PM (driving end-to-end per owner decision)
- **Status:** design - pending user review before plan
- **Predecessor:** [#857](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/857) (`AGENTS.md` canonical, `CLAUDE.md` symlink) - this is the directory-level sibling.

## 1. Goal & honest scope

Make `.agents/` the canonical project-local agent-config directory and turn `.claude/` into a committed symlink pointing at it (`.claude/ -> .agents/`).
This extends the vendor-neutral canonicalization already done for instructions down to the config directory: the source of truth stops carrying the `claude` label, while Claude Code keeps working unchanged through the symlink.

This is an **organizational / naming** change, **not** new functional portability.
Roughly 80% of the directory's contents (`settings.json`, `hooks/*.py`) are Claude Code's runtime/event model and stay Claude-specific regardless of the folder name - no other tool reads them.
The win is a vendor-neutral canonical path plus a spec-aligned home (`.agents/skills/`) for the one genuinely portable surface (skills), matching the project's agnostic-by-intent direction.

**Non-goals:**

- No change to hook *logic* - only their declared command paths.
- No skills are moved: none exist yet under `.claude/skills/` (the `authoring-research-decks` skill is still only a plan doc). We merely bless `.agents/skills/` as the future home.
- No attempt to make Codex/Cursor/etc. read our hooks or settings.

## 2. Current state (verified 2026-06-25, all three clones)

All three clones are separate clones (not worktrees), all on `main`, all clean.

**Tracked under `.claude/` (migrates via git, propagates on pull):**

- `commands/coordination.md`
- `hooks/check_at_claude.py`, `hooks/check_board_query_pagination.py`, `hooks/check_gh_issue_develop_parent.py`, `hooks/post_gh_pr_create.py`, `hooks/recheck_dispatch.py`
- `settings.json`

**Per-clone local, gitignored (does NOT migrate via git - handled per clone):**

| File | Nature | Handling on receiving clone |
|---|---|---|
| `memory` → `../../claude-personas-splice-neoepitope-pipeline/<role>` | role-specific symlink | **must be re-created** at `.agents/memory` with the same target |
| `settings.local.json` | per-user permissions allow-list | preserve (move aside → restore) |
| `hook_fires.jsonl` | guard fire log | droppable (regenerates) |
| `hooks/__pycache__/` | transient bytecode | droppable (regenerates) |
| `scheduled_tasks.{lock,json}` | transient scheduler state | droppable if present (regenerates) |

Per-clone `memory` symlink targets (the script reads the existing target rather than hardcoding a role):

| Clone | `memory` target |
|---|---|
| `splice-neoepitope-pipeline` (base) | `…/developer` |
| `splice-neoepitope-pipeline-pm` | `…/pm` |
| `splice-neoepitope-pipeline-scientist` | `…/scientist` |

**Hook command paths today** are relative (`.claude/hooks/X.py`), not `${CLAUDE_PROJECT_DIR}`-anchored.

## 3. Why two clones can't just `git pull` it cleanly

Git tracks **commits**, not working-tree appearance. The migration is one tracked commit (the `git mv` rename + the committed `.claude` symlink). When a receiving clone pulls that commit, git must remove the old tracked files from `.claude/` and lay a `.claude` **symlink** in their place - but git **refuses to replace a directory that still holds untracked files** with a symlink. Each receiving clone's `.claude/` still physically holds its own gitignored locals (the `memory` symlink, `settings.local.json`, the fire log, `__pycache__`), so the pull half-applies or aborts.

Pre-renaming locally before the pull does **not** help: a manually-`mv`'d `.agents/` is *untracked*, and the incoming commit wants to create those same paths as *tracked* files, so `git pull` aborts with `The following untracked working tree files would be overwritten by merge`. The only clean sequence is **clear the untracked locals out of `.claude/` → pull (git does the tracked swap) → restore the locals under `.agents/`**.

## 4. Architecture / approach

Two distinct execution paths, one shared idea (rename + symlink-back):

### 4a. Originating clone (PM) - authored by hand in the PR branch

1. `git mv .claude/commands .agents/commands`
2. `git mv .claude/hooks .agents/hooks` (carries the untracked `__pycache__` along physically; still ignored)
3. `git mv .claude/settings.json .agents/settings.json`
4. Relocate this clone's gitignored locals into `.agents/`: move `settings.local.json` and `hook_fires.jsonl`; record the `memory` target and re-create `.agents/memory -> <same target>`; drop `__pycache__`.
5. `rmdir .claude` (now empty) and `ln -s .agents .claude`.
6. Rewrite `.agents/settings.json` hook command paths (§5).
7. Rewrite `.gitignore` (§6).
8. `git add .agents .claude .gitignore` (+ staged deletions); commit; open PR.

The committed objects: the renamed tracked files under `.agents/`, the `.claude` symlink (git mode 120000), and the updated `.gitignore`.

### 4b. Receiving clones (base, scientist) - committed one-shot script

`scripts/migrate_claude_to_agents.sh` - idempotent + self-verifying. It reads the *existing* clone's role from the current `memory` symlink target, so the same script is correct in every clone.

Procedure:

1. **Preconditions:** assert inside a git repo with a clean tracked tree; assert `.claude/` is still a real directory (else exit 0 - already migrated, idempotent).
2. **Save state:** read and store the `memory` symlink target; move `settings.local.json` aside to a temp; remove the droppable locals (`hook_fires.jsonl`, `hooks/__pycache__/`, `scheduled_tasks.{lock,json}`) so nothing untracked remains under `.claude/`.
3. **Pull:** `git pull --ff-only origin main` - git removes the now-only-tracked files from `.claude/` and lays down `.agents/` + the `.claude` symlink.
4. **Restore:** move `settings.local.json` into `.agents/`; re-create `.agents/memory -> <saved target>`.
5. **Verify** (fail loudly on any miss):
   - `.claude` is a symlink resolving to `.agents`;
   - `.agents/memory` resolves to a directory whose basename matches the saved role;
   - `.agents/hooks/*.py` and `.agents/settings.json` exist;
   - `git status --porcelain` is empty (no stray diff, no accidentally-tracked local).

The script is disposable post-rollout (like prior one-shot migration helpers); it stays in-tree as the documented procedure and AC evidence.

## 5. `settings.json` hook-path rewrite

Rewrite each hook `command` from the relative `.claude/hooks/X.py` to the absolute canonical path:

```
${CLAUDE_PROJECT_DIR}/.agents/hooks/check_at_claude.py
${CLAUDE_PROJECT_DIR}/.agents/hooks/check_gh_issue_develop_parent.py
${CLAUDE_PROJECT_DIR}/.agents/hooks/check_board_query_pagination.py
${CLAUDE_PROJECT_DIR}/.agents/hooks/post_gh_pr_create.py
${CLAUDE_PROJECT_DIR}/.agents/hooks/recheck_dispatch.py --scope shared
```

Rationale (belt-and-suspenders): guard **execution** then resolves the real `.agents/` path directly and never depends on symlink resolution; only Claude Code's **discovery** of `.claude/settings.json` rides the symlink. The `--scope shared` arg on `recheck_dispatch.py` is preserved verbatim.

## 6. `.gitignore` rewrite

Repoint the per-clone artifact rules to the canonical `.agents/` paths (keep the `.claude/` forms too - harmless, both resolve to the same inode through the symlink, and they protect any clone mid-migration):

```gitignore
# Claude Code - transient scheduler state (per-session lock + durable cron jobs)
.agents/scheduled_tasks.lock
.agents/scheduled_tasks.json
.claude/scheduled_tasks.lock
.claude/scheduled_tasks.json

# Claude Code - per-user local settings (permissions allow-list, absolute paths)
.agents/settings.local.json
.claude/settings.local.json

# claude-personas role-memory symlink
/.agents/memory
/.claude/memory
.agents/hook_fires.jsonl
.claude/hook_fires.jsonl
```

`.claude` (the symlink itself) is **not** ignored, so it commits. The generic `__pycache__/` rule (line 5) already covers `.agents/hooks/__pycache__/`.

## 7. Verification & rollout sequence

1. PM clone: author 4a on a branch (`gh issue develop 861` via `scripts/new_branch.sh`), open PR, request `@-claude review`, merge via `scripts/audit_and_merge.sh`.
2. ⚠️ **Mid-session hook death:** rewriting `settings.json` silently deactivates all hooks for the *current* session (documented gotcha - no hot-reload). So the live-hook ACs are verified in a **fresh session after merge**, not in the authoring session.
3. Fresh-session verification (PM clone): `/hooks` shows all three PreToolUse guards; `/doctor` clean; a trace-probe confirms a guard actually fires; `/reload-plugins` then confirm any future `.agents/skills/` discovery path.
4. Receiving clones: run `scripts/migrate_claude_to_agents.sh` in base + scientist; confirm the script's self-verify passes; spot-check role memory loads and a guard fires in each.
5. **Docs sweep:** update lingering `.claude/` references in `docs/` / `AGENTS.md` where they should point at `.agents/`; confirm any that stay (symlinked) still resolve. This is the relative-link-breakage class that bit the #860 re-review - verify each, don't assume.

## 8. Acceptance criteria (from #861, mapped)

- [ ] `.agents/` is the real directory; `.claude/` is a committed symlink (git mode 120000).
- [ ] `/hooks` shows all three PreToolUse guards live and `/doctor` is clean in a fresh post-migration session.
- [ ] Hook command paths point at `${CLAUDE_PROJECT_DIR}/.agents/hooks/...`; a trace-probe confirms a guard fires.
- [ ] `.agents/skills/` is the blessed location and is discovered when populated (verify after `/reload-plugins`).
- [ ] `.gitignore` covers the per-clone artifacts under the new path; no gitignored file is accidentally tracked.
- [ ] Each clone's role-memory symlink is re-created at `.agents/memory` and resolves to that clone's own role dir; role memory still loads.
- [ ] `docs/` / `AGENTS.md` references to `.claude/` still resolve (via symlink or updated to `.agents/`).
- [ ] The per-clone migration is encoded in `scripts/migrate_claude_to_agents.sh`, and **all three clones** (base, `-pm`, `-scientist`) are migrated + re-verified.

## 9. Risks & mitigations

| Risk | Mitigation |
|---|---|
| Half-migrated clone → silently dead guards or broken role memory | Idempotent, self-verifying script that fails loudly; fresh-session `/hooks` + trace-probe check |
| `git pull` aborts on populated `.claude/` | Script clears untracked locals before pull (§3, §4b) |
| `${CLAUDE_PROJECT_DIR}` not set in some context | Discovery still works via the `.claude` symlink; execution path is absolute - both belt and suspenders |
| Doc relative-links break on path change | Explicit per-reference verify in the docs sweep (§7.5) |
| Local `settings.local.json` lost during swap | Preserved via move-aside → restore; never deleted |
