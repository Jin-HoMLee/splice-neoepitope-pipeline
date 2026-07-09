# Memory-path cwd-drift guard (Issue #1053)

Design spec, 2026-07-08. Developer.

## Problem

The four role agents each work from their own clone. The main-session Bash shell
persists its working directory across tool calls (the standard persistent-shell
model), so the cwd can silently drift out of the clone root - most commonly via a
legitimate `cd` into the scratchpad temp root that is never undone. A later
command that uses a **clone-root-relative** path (e.g. `grep .agents/memory/...`)
then resolves against the wrong directory.

## Premise audit (AC#1)

Auditing the three recorded drift incidents (esp. the 2026-07-04 PM episode)
resolves the failure shape the issue left open ("persisted cwd vs realpath
mis-resolution vs both"):

- **The recurring shape is persisted cwd drift, not wrong-clone memory writes.**
  On 2026-07-04 a `cd .../scratchpad` persisted and made a later bare-relative
  `grep .agents/memory/shared` resolve wrong. The PM's own clone-realpath check
  that session confirmed edits still landed in the *right* clone.
- **The "silent wrong-file memory I/O" class has never actually fired** - it is a
  latent risk warned about by the bare-path memory convention, with zero observed
  occurrences.
- **The existing `check_no_cd_outside_cwd` hook (PR #1029) does not cover this.**
  It denies `cd` into another repo / sibling clone / `$HOME`, but *deliberately
  allows* `cd` into `/private/tmp` (the scratchpad root) - exactly the door the
  drift walks through. It landed hours after the last incident.
- **PreToolUse JSON carries `cwd`**, which is the lever a residual guard needs.

## Decision

Build a narrow residual guard as cheap insurance against the latent-but-serious
wrong-file memory-write class, complementing (not duplicating) the existing
cd-guard. Grounded in agentic-coding best practice: persistent-shell harnesses
accept cwd drift as a known tradeoff and mitigate it by (a) preferring absolute
paths / `git -C` and (b) making drift *loud* rather than silent. This guard makes
the one high-harm drift moment loud.

Rejected: a blanket "warn on any command while cwd drifted" guard - it would fight
the existing hook's deliberate temp-root allowance and generate noise. Rejected:
closing the issue - the wrong-file class, though never fired, is a genuine
memory-corruption hazard and the guard is cheap.

## Design

New `PreToolUse` hook `.agents/hooks/check_memory_path_cwd_drift.py`.

**Registered on** the `Bash` matcher and the `Edit|MultiEdit|Write` matcher (two
settings entries, one script).

**Denies** iff *both* hold:
1. `cwd` from the hook JSON resolves **outside the clone root** (reuse the
   existing hook's root-derivation: `Path(__file__).resolve().parent.parent.parent`,
   symlink-safe, no env var). A cwd inside the clone root or unresolvable → allowed.
2. The tool references a **relative** path under `.agents/memory/` or
   `.claude/memory/`:
   - `Bash`: scan the command string for those segments as a relative token.
   - `Edit`/`MultiEdit`/`Write`: inspect `file_path` (and for `Bash`, the whole
     command). An **absolute** memory path is always safe → allowed.

**Fails open** (allow) on: unresolvable cwd, dynamic paths (`$VAR`, backticks,
command substitution), tokenization failure, unexpected JSON shape, missing
fields. A guard must never break a legitimate command on ambiguity.

**Escape hatch**: `CLAUDE_ALLOW_CWD_DRIFT=1`.

**Fire log**: one line per deny to `.agents/hook_fires.jsonl` (Issue #453 infra),
`hook: "check_memory_path_cwd_drift"`.

**Deny message**: names the drift, shows the drifted cwd vs clone root, and nudges
to the real fix - `cd` back to the clone root or use an absolute path / `git -C`.

## Testing

Unit tests in `tools/ci/test_check_memory_path_cwd_drift.py`, mirroring
`test_check_no_cd_outside_cwd.py`: pure helpers (path classification, memory-path
detection, cwd-inside-clone) plus `main()` end-to-end via stdin JSON. Cover:
drift + relative memory path (deny), drift + absolute memory path (allow), drift +
non-memory path (allow), no-drift + relative memory path (allow), dynamic path
(allow / fail-open), escape hatch, each tool type, fire-log line shape, and an
`os.access(HOOK, os.X_OK)` regression assertion (the exec-bit lesson from #1032).

## Acceptance criteria mapping

- AC1 (confirm premise) → this spec's Premise audit section + a comment on #1053.
- AC2 (guard fires on a probe) → live-fire probe in the wiring step.
- AC3 (fail-open + fire-log) → design + tests above.
- AC4 (strip stopgap bullet) → same PR, cross-repo edit to personas `pm/MEMORY.md`.
