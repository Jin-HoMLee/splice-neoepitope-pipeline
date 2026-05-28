# Hook Fire-Log Infrastructure + Threshold-Crossing Promotion Prompt — Design

**Issue:** [#453](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/453)
**Branch:** `feat/pm/issue-453-hook-fire-log-infra`
**Audience:** Implementer of #453 (PM or whoever picks it up).

## 1. Goal

Replace the "vibes-check" gate for hook promotion (`.claude/settings.local.json` → committed `.claude/settings.json`) with deterministic instrumentation. Every warning fire writes one line to `.claude/hook_fires.jsonl`; an aggregator script summarizes; a threshold-crossing helper surfaces the promotion-decision prompt automatically when a hook hits K=3 fires for the first time.

## 2. Background / why

Three sibling PostToolUse rechecks live in `.claude/settings.local.json` today (PM-clone only, gitignored):

| Hook | Logic location | Filed via |
|------|---------------|-----------|
| `target_sync_check` | Inline function in [`.claude/hooks/recheck_dispatch.py:153-185`](../../../.claude/hooks/recheck_dispatch.py#L153) | [PR #450](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/450) / [Issue #448](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/448) |
| `recheck_milestone` | [`scripts/pm/recheck_milestone.py`](../../../scripts/pm/recheck_milestone.py), subprocess-invoked | [PR #397](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/397) / [Issue #247](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/247) |
| `recheck_parent_status` | [`scripts/pm/recheck_parent_status.py`](../../../scripts/pm/recheck_parent_status.py), subprocess-invoked | [PR #407](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/407) / [PR #415](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/415) |

Each is a candidate for eventual promotion to committed `.claude/settings.json` (so Scientist + Developer sessions also benefit). The "should we promote this?" decision today is vibes-check. This infra replaces that with counted-fires data + a threshold-crossing nudge.

The structural asymmetry (target_sync_check inline vs sibling scripts) is intentional, tracked in [Issue #533](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/533) as separate cleanup follow-up. This Issue does NOT depend on it being resolved.

## 3. Architecture overview

3 helper functions + 1 wiring change inside `recheck_dispatch.py`, plus one new bash aggregator:

```
.claude/hooks/recheck_dispatch.py
├── _log_fire(hook_name, issue, **metadata) → None
│   └── appends one JSONL line to .claude/hook_fires.jsonl (atomic)
├── _count_fires(hook_name) → int
│   └── reads jsonl, substring-counts lines matching hook_name
├── _threshold_prompt(hook_name) → str | None
│   └── if count == K (exactly): returns "🎯 …" line; else None
├── _is_fire(hook_name, output) → bool
│   └── per-hook predicate distinguishing fire vs clean output
├── _wrap_warning(hook_name, issue, warning, outputs) → None
│   └── always appends warning to outputs (so user sees context);
│       if _is_fire: logs fire + appends threshold prompt to outputs
└── dispatch() — call sites refactored to use _wrap_warning

scripts/check_hook_health.sh  (NEW)
└── reads .claude/hook_fires.jsonl, summarizes per-hook fires + dock state

.gitignore
└── .claude/hook_fires.jsonl  (NEW LINE)
```

## 4. Components

### 4.1 Data shapes

**JSONL line schema** (single line per fire, append-only, never edited):
```json
{"ts": "2026-05-28T15:42:11Z", "hook": "target_sync_check", "issue": 533, "ms_due_on": "2026-06-11", "target_date": null}
```

Required keys: `ts` (ISO 8601 UTC), `hook` (string), `issue` (int or null). Optional hook-specific kwargs land as extra keys via kwargs passthrough. No enforced per-hook schema — keeps the writer dead-simple.

**Per-hook config** — module-level constants in `recheck_dispatch.py`:
```python
HOOK_CONFIG = {
    "target_sync_check":     {"threshold": 3, "dock": 454},
    "recheck_milestone":     {"threshold": 3, "dock": None},  # to be filed when promotion is sought
    "recheck_parent_status": {"threshold": 3, "dock": None},  # to be filed when promotion is sought
}
```

K defaults to 3 across all hooks but is per-hook configurable. `dock=None` → threshold prompt says `"(no dock Issue filed yet)"`; `dock=N` → `"Issue #N"`.

**Decision: constants, not a JSON config file.** YAGNI — 3 entries that change rarely. Move to a separate config file only if a 4th+ hook joins.

### 4.2 Concurrency + counting

**Append strategy — POSIX line-atomic write:**
```python
LOG_PATH = Path(__file__).resolve().parent.parent.parent / ".claude" / "hook_fires.jsonl"

def _log_fire(hook_name: str, issue: int | None = None, **metadata) -> None:
    payload = {"ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
               "hook": hook_name, "issue": issue, **metadata}
    line = json.dumps(payload, separators=(",", ":")) + "\n"
    if len(line.encode("utf-8")) >= 4096:
        print(f"hook_fires: oversized line (>=4KB) for {hook_name}, skipped", file=sys.stderr)
        return
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(line)
```

POSIX guarantees each `write()` syscall on an `O_APPEND` regular file is atomic — no interleaved bytes from concurrent writers, regardless of length. **No `fcntl.flock` needed.** Our lines are ~150-200 bytes; the 4 KB guard is sanity hygiene to keep the log structure parseable if a future schema change pushes a line above ~4 KB (the write is skipped rather than persisted; stderr warning surfaces the issue).

**Counting — naive substring match:**
```python
def _count_fires(hook_name: str) -> int:
    if not LOG_PATH.exists():
        return 0
    needle = f'"hook":"{hook_name}"'
    return sum(1 for line in LOG_PATH.read_text().splitlines() if needle in line)
```

K=3 threshold + low fire trajectory → log stays tiny for weeks. JSON-parse-per-line wins nothing over substring match at this scale. **No sidecar count file.** (Note: `json.dumps` with `separators=(",", ":")` ensures the needle pattern matches the serialized output exactly — no whitespace drift.)

### 4.3 Fire predicate + wiring

**Per-hook fire predicate:**
```python
def _is_fire(hook_name: str, output: str) -> bool:
    if hook_name == "target_sync_check":
        return bool(output)  # function returns None when no fire
    if hook_name in ("recheck_milestone", "recheck_parent_status"):
        return "Status: [No change]" not in output  # both scripts share this sentinel
    return False
```

Both subprocess scripts emit a `Status: [No change]` line when no drift/staleness is detected (verified in [`scripts/pm/recheck_parent_status.py:158`](../../../scripts/pm/recheck_parent_status.py#L158); same convention in `recheck_milestone.py`). Predicate-matched directly in dispatch.py — scripts stay untouched.

**Threshold prompt:**
```python
def _threshold_prompt(hook_name: str) -> str | None:
    cfg = HOOK_CONFIG.get(hook_name)
    if cfg is None:
        return None
    if _count_fires(hook_name) != cfg["threshold"]:  # == K, not >=, so fires exactly once
        return None
    dock = f"Issue #{cfg['dock']}" if cfg["dock"] else "(no dock Issue filed yet)"
    return f"🎯 {hook_name} has now fired {cfg['threshold']} times — review for promotion: bash scripts/check_hook_health.sh ; {dock}"
```

Equals-once (`==`, not `>=`) means the prompt fires on exactly the K-th fire and never again. No nagging.

**Wiring helper:**
```python
def _wrap_warning(hook_name: str, issue: int | None, warning: str, outputs: list[str]) -> None:
    outputs.append(warning)  # always emit so user sees recheck context
    if not _is_fire(hook_name, warning):
        return
    _log_fire(hook_name, issue=issue)
    prompt = _threshold_prompt(hook_name)
    if prompt:
        outputs.append(prompt)
```

**`dispatch()` call-site refactor:** every existing `outputs.append(...)` for a check's output becomes `_wrap_warning(hook_name, issue, output, outputs)`. The `dispatch()` function shrinks slightly (no behavioral change for clean / non-fire paths).

For `target_sync_check` (inline function returning `str | None`): the call sits at the existing call-site (line 235 area); the function's return value flows directly into `_wrap_warning`.

### 4.4 Aggregator — `scripts/check_hook_health.sh`

Operator-facing inspection tool. Pure bash + `jq`.

**Output shape:**
```
Hook health summary (.claude/hook_fires.jsonl)

target_sync_check
  Fires:  5
  First:  2026-05-21T10:14:33Z
  Last:   2026-05-28T14:25:47Z
  Dock:   Issue #454

recheck_milestone
  Fires:  12
  First:  2026-05-13T09:02:11Z
  Last:   2026-05-28T15:42:11Z
  Dock:   (no dock Issue filed yet)

recheck_parent_status
  Fires:  0
  Dock:   (no dock Issue filed yet)

Total: 17 fires across 3 hooks
```

**No-fires-yet case:** when `.claude/hook_fires.jsonl` doesn't exist, prints `"No fires logged yet."` + per-hook config snapshot. Exit code 0 (informational).

**Dock display:** hardcoded 3-entry mapping at top of script (acceptable duplication with HOOK_CONFIG in dispatch.py at 3 entries; refactor only if 4th hook joins).

### 4.5 Testing

Test file: `workflow/tests/test_recheck_dispatch.py` (creates the file if not present). Run via `workflow/tests/.venv/bin/python -m pytest workflow/tests/test_recheck_dispatch.py -q`.

**Three test cases** (covering the 3 ACs that mention tests):

1. **`test_log_fire_jsonl_append`** — call `_log_fire("hook_A", issue=1)` 3× into a tmp_path log → file has 3 lines, each valid JSON, each contains `"hook":"hook_A"` and `"issue":1` and a valid ISO ts.

2. **`test_count_fires_filters_by_hook`** — log 2 × `hook_A` + 1 × `hook_B` → `_count_fires("hook_A") == 2`, `_count_fires("hook_B") == 1`, `_count_fires("hook_C") == 0`.

3. **`test_threshold_prompt_equals_once`** — log 2 × fires → `_threshold_prompt` returns None. Log 3rd → returns the 🎯 line (and contains "fired 3 times"). Log 4th → returns None.

Tests use `monkeypatch` to redirect `LOG_PATH` and `HOOK_CONFIG` to tmp_path-scoped values, so they don't touch the real log.

**No tests for the bash aggregator** (glue script; manual inspection sufficient; if it grows logic, retrofit).

### 4.6 Cross-cutting

**`.gitignore` append:**
```
.claude/hook_fires.jsonl
```

**Shared memory cross-ref — `shared/feedback_hook_proves_out.md`:**
Update to cite Issue #453 as the landed infra + add a one-line "How to apply" pointing operators at `bash scripts/check_hook_health.sh` before any future hook-promotion decision. (Personas-repo edit; lands via direct commit per shared-rule "personas-repo git state is not your responsibility".)

## 5. Open questions / risks / non-goals

### Open questions

None blocking. Implementer should verify on first run:
1. `Status: [No change]` sentinel string is preserved on any future `recheck_milestone.py` / `recheck_parent_status.py` refactor — the fire predicate depends on it. (Mitigation: a passing unit test plus the test_threshold_prompt_equals_once test exercises the predicate via mock outputs.)

### Risks

1. **Schema growth past 4KB** — if a future hook logs verbose metadata pushing the line over PIPE_BUF, the oversized-line guard skips the write (no log corruption, but loses that fire). Mitigation: keep metadata terse; if it grows, switch to `fcntl.flock`.
2. **Test-environment LOG_PATH leak** — if `monkeypatch` doesn't fully isolate `LOG_PATH`, tests could pollute `.claude/hook_fires.jsonl`. Mitigation: use tmp_path fixture; assert log_path is under tmp before writing.

### Non-goals (explicit)

- Promoting any specific hook to committed `.claude/settings.json` — per-hook dock Issues (e.g. [Issue #454](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/454)) handle that decision separately.
- Wiring fire-log into non-recheck hooks (e.g. `check_at_claude.py`) — different lifecycle and rationale; retrofit if useful.
- Resolving the inline-vs-script asymmetry of `target_sync_check` — tracked separately at [Issue #533](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/533); not blocking.
- Reading `HOOK_CONFIG` at runtime from the bash aggregator (no Python-introspection round-trip; tiny manual duplication is acceptable at 3 entries).

## 6. Estimated effort

M-sized. Single implementation pass: ~3-4 hours focused work.

**Project-repo PR** (closes [Issue #453](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/453)) — 4 file changes:
- `.claude/hooks/recheck_dispatch.py` (5 helpers added + `dispatch()` call-sites refactored)
- `scripts/check_hook_health.sh` (new file)
- `workflow/tests/test_recheck_dispatch.py` (new file)
- `.gitignore` (1-line append)

**Personas-repo direct commit** (separate from project-repo PR; per shared-rule "personas-repo git state is not your responsibility"):
- `shared/feedback_hook_proves_out.md` cross-references [Issue #453](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/453) + adds "How to apply" pointing to `check_hook_health.sh`.

**Created by:** PM
