# Hook-liveness contract tests - design

Issue [#1140](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1140) (Sub D of epic [#1135](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1135)).
Date: 2026-07-15.
Author: Developer.

## The class this closes

Four governance mechanisms shipped **completely inert** this month; each passed its unit tests, none was ever invoked through its real path.
The common shape: the test exercised everything downstream of the trigger and nothing at it, so it could only ever come back green.
`post_gh_pr_create`'s `matches_pr_create()` returned `False` for a heredoc-then-`gh pr create` (the way every PR with a real body is opened) and was dead for months while every unit test passed.

Mutation testing (the sibling, Issue #1141) structurally cannot see this class: a hook that is never invoked has no surviving mutant, because no test calls it through its real path.
This mechanism is the half that catches an inert hook.

## What we build

A **hook-liveness contract suite**: for every hook registered in `.agents/settings.json`, a test that

1. constructs the hook's **real trigger** as the harness would deliver it (the real command *string*, including the heredoc-then-command shape),
2. drives it through the hook's **real entry path** - the hook script as a subprocess reading the harness JSON envelope on stdin, exactly as the harness invokes it, not by calling an inner function or `main()` with a hand-built payload,
3. asserts the **observable artifact** appears (a `permissionDecision: deny` on stdout, or an appended `.agents/hook_fires.jsonl` line),
4. and, as a matched-pair control, asserts a should-not-fire counterexample produces **no** artifact.

**Absence of the artifact is the failure signal.** The suite must be able to come back red when a hook is inert.

## Why "real entry path" = subprocess-with-real-command, not the harness `if:` gate

The harness applies the `if:` matcher (e.g. `Bash(gh *)`) before invoking the hook command.
That gate is vendor-documented **best-effort** and, per project convention, every hook self-matches internally and fails open - so the gate is a performance optimization, not a correctness boundary.
We therefore do **not** re-implement the harness gate (it is harness code we cannot see and would only be guessing at).
We drive the hook script directly, piping the true command string in the real envelope - which is exactly what the harness does once the gate passes, and it is the layer where every one of the four dead hooks actually failed (their own internal matcher, `matches_pr_create` et al.).

This is the concrete engineering form of the standing AC: **drive the real trigger, never pipe a synthetic payload into `main()`.**

## Architecture

Two files, one CI job.

### `tools/ci/hook_contract.py` - the shared harness (single-authored)

- `registered_hooks()` - parse `.agents/settings.json`, return one record per `(event, matcher, if, command)` registration. This is the enumeration AC and the drift guard: the suite is parametrized off *this*, so a newly registered hook with no contract entry fails the completeness test rather than being silently uncovered.
- Envelope builders per event: `pretooluse_bash(command)`, `pretooluse_edit(tool, file_path, old, new, content)`, `posttooluse_bash(command, tool_response)`, `stop_envelope(...)` - each returns the JSON dict the harness pipes on stdin for that event/matcher.
- `drive(hook_path, envelope, *, gh_stub=None, env=None)` - run the hook as a subprocess (`sys.executable hook_path`), stdin = `json.dumps(envelope)`, capture stdout/stderr/exit. When `gh_stub` is set, prepend a temp dir holding a stub `gh` to `PATH` (the `test_set_status.py` pattern) so board hooks reach their handler without a network call.
- Observable detectors: `deny_decision(result)` reads stdout JSON for `permissionDecision == "deny"`; `fired(before, after)` diffs `.agents/hook_fires.jsonl` line count around a `drive()` call. A `fire_log_guard()` context manager snapshots the real (gitignored) log so the suite leaves it untouched.

### `tools/ci/test_hook_liveness_contract.py` - the suite

- A per-hook **registry**: `{hook_basename: HookContract(...)}` where each `HookContract` carries the real fire command, the real no-fire counterexample, the observable kind (`DENY` | `FIRE_LOG`), and whether a `gh` stub is needed. Populated per hook (the fan-out output below).
- `test_every_registered_hook_has_a_contract()` - the completeness AC: `set(registered_hooks basenames) == set(registry)`; an unregistered-but-present hook and a registered-but-uncovered hook both fail here.
- Parametrized `test_hook_fires_on_real_trigger(contract)` and `test_hook_silent_on_counterexample(contract)` - the matched pair, per hook.
- `test_heredoc_gh_pr_create_shape()` - explicit coverage of the exact string that killed two hooks.
- `test_breaking_a_matcher_turns_the_contract_red()` - the red-on-break **demonstration**, not assertion: copy one hook to a temp file, break its matcher, drive the real trigger through the copy, assert the contract detector now reports no fire. Demonstrated, not claimed.

### CI

Extend the existing hook-test CI job (or add a path-filtered one) to run the suite on any change under `.agents/hooks/` or `.agents/settings.json`.

## Per-hook fan-out

Ten distinct hook scripts, each with a different trigger shape, observable, and counterexample.
Determining each hook's *real* fire trigger and a genuine no-fire counterexample requires reading that hook's own matcher - independent per-hook source analysis, which is the natural unit of parallel work.
Each hook is analyzed by one agent (produce the `HookContract` spec) and adversarially verified by a second (confirm the fire command truly fires and the counterexample truly does not, by reasoning about the hook's real matcher - the check against a wrong trigger, the `matches_pr_create` failure mode).
The main session assembles the verified specs into the single registry and runs the suite.

## The live-`gh` wrinkle

`post_gh_pr_create`, `post_gh_pr_review_request`, and `recheck_dispatch` write their fire-log line only *after* a live `gh` mutation.
Their contract drives the real trigger with a stub `gh` on `PATH` (returning canned board/PR JSON), so the matcher-plus-handler path runs end-to-end and the fire-log append is the observable - proving liveness without a network call.
`write_session_watermark` (Stop event) writes a watermark file, not a deny/fire-log line; its observable is that file appearing, driven with a temp `CLAUDE_PROJECT_DIR`.

## Out of scope (surfaced, per AC 1)

The `@claude` mention guard (`check_at_claude.py`) is registered at **user** level (`~/.claude/settings.json`), not in `.agents/settings.json`, so it is outside this suite's registry-driven scope.
That it is unreachable from the project registry is itself worth noting (AC 1: "a hook absent from the registry is a liveness bug worth surfacing"); we record it as a known gap rather than silently omitting it.

## Testing

The suite *is* the test. Its own credibility rests on `test_breaking_a_matcher_turns_the_contract_red` - if that demonstration cannot make the suite red, the suite is itself a hollow check and the whole exercise has failed.
