# Visual Morning Routine — Design Spec

**Date:** 2026-06-04
**Status:** Draft for review (brainstorming output → next: writing-plans) — **revised after 7-lens adversarial review**
**Issues** (split by work-type 2026-06-04): Design — this spec + the Phase-1 plan — is tracked by [Issue #656](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/656) (`role:pm`); the implementation epic is [Issue #655](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/655) (`role:developer`), whose build sub-issues are carved at execution.
**Author:** PM session (Claude) + Jin-Ho Lee

## 1. Goal

Make the daily morning routine **visual, interactive, and coupled to the chat session** — delivered as a reusable Claude Code **skill** that all three roles (PM, Scientist, Developer) can use. Today the routine is Markdown rendered in the terminal, paced one phase per message with a live TodoWrite list ([`feedback_morning_routine.md`](../../../.claude/memory/feedback_morning_routine.md), [`shared/feedback_morning_routine.md`](../../../.claude/memory/shared/feedback_morning_routine.md), [`shared/feedback_morning_routine_pacing.md`](../../../.claude/memory/shared/feedback_morning_routine_pacing.md)). This adds a browser cockpit beside the chat that mirrors the session live and lets the user click to drive it — **without replacing** the chat transcript, which stays the durable record and the place real work (gh, board edits, governance) happens. A user who ignores the browser still gets the full paced text routine.

This spec is the *what/why*. The step-by-step build plan is produced next by the writing-plans skill.

## 2. Decisions locked (from brainstorming, 2026-06-04)

| # | Decision | Choice |
|---|----------|--------|
| D1 | Coupling direction | **Bidirectional** (session ↔ page) |
| D2 | Center of gravity | **Co-equal** — chat and page both first-class |
| D3 | Form factor | **Cockpit** — pinned summary bar + left progress rail + focus panel + footer |
| D4 | Liveness | **Hybrid** — page-local nav acts live; **any state-mutating action confirms in chat first** |
| D5 | Launch | Built **as a skill** (`/morning`); "good morning" routes to it **via an instruction rule, not harness auto-dispatch** (see §2.2) |
| D6 | Generalization | **PM-first, config-driven** — engine is role-agnostic; News/Status need per-role formatters (see §11 for the honest scope) |
| D7 | Build substrate | **Purpose-built local stdlib server** in the repo (not the superpowers companion infra; not static HTML) |
| D8 | Concurrency | **One server per role/clone** (not one shared server) |
| D9 | Item references | Issues / PRs / milestones / board filters render as **hyperlinks to GitHub** (third interaction type) |

### 2.1 Spike status — live bridge PROVEN on the purpose-built substrate (Phase-0 run 2026-06-04)

Two spikes, run during design:

**Spike 1 (browser → session, n=1):** a `run_in_background` watcher exited on one real browser click and the task-notification re-invoked the idle session with zero user input (armed 09:33:18, fired 09:35:40). Establishes the browser → server → events → watcher → wake path end-to-end, once.

**Spike 2 — hardened Phase-0 (the re-arming loop, on the purpose-built stdlib server):** a minimal `morning_server.py` (token-gated `POST /event`, server-stamped `seq`, `/health` nonce, redacted logs, dynamic port) received 9 auto-driven HTTP POSTs; a re-arming watcher (exit on `seq > cursor`, re-armed each wake) drove the session through **6 consecutive auto-wakes.** Measured:

| Metric | Result |
|---|---|
| Consecutive auto-wakes delivered | **6 / 6 (100%)** — no degradation, no dropped/coalesced notifications across cycles |
| Events consumed (cursor dedup) | **9 / 9, zero lost** |
| Miss-rate | **0%** (go/no-go abort threshold was > 5%) |
| Coalescing | Confirmed — events arriving faster than the cycle collapse into one wake (wakes 1–3 consumed 2 events each); spaced wider than latency → clean 1:1 (wakes 4–6) |
| End-to-end latency (event → session acts) | **~7–35 s, variable** (fastest clean cycles 7–8 s) |
| Server security gates | token → 403, bad content-type → 415, **neither appended** ✓ |

**Verdict: GO on the full live bridge** (not the §8.4 batch-reconcile fallback). Two honest caveats remain: (a) the loop was driven by `curl` POSTs, not real browser clicks — but Spike 1 proved the browser leg and this server token-gated + stamped 9 real HTTP POSTs, so the composition is sound (real-browser *cadence* is the only untested seam, re-verified in build Phase 0); (b) **token cost is real** — each wake is a full model turn on a growing transcript, so the §8.2 wake-minimisation levers (nav is page-local with **no** wake; coalescing; builder caching) are **mandatory, not optional** — coalescing alone already turned 9 events into 6 wakes. The ~7–35 s latency means the cockpit is "live within seconds," not snappy-instant; the UX is framed accordingly.

### 2.2 How "good morning" routing actually works (D5 correction)

There is **no harness-level dispatcher** that intercepts a plain-text "good morning" and auto-invokes a skill. Skills are invoked via the **Skill tool by the already-running session**. Today "good morning" works because an Always-in-effect memory rule instructs the model to load the routine. So routing is **instruction-driven and best-effort**: the existing memory rule must be amended (personas repo, MM-committed — §17) to say "invoke the `/morning` skill," and that edit **must land before the skill is usable**. The reliable primary entry is the user typing **`/morning`**; "good morning" is a soft alias. This is the same class of instruction-nudge that CLAUDE.md documents as occasionally slipping — so we do not treat auto-invoke as guaranteed.

## 3. Architecture

Five components. All code lives in the repo; all runtime state lives in a gitignored per-clone dir.

| Component | What | Role |
|---|---|---|
| **`morning` skill** | `.claude/skills/morning/SKILL.md` (entry) + Python under `tools/morning/` | Orchestration; `/morning` (and the "good morning" alias) routes here |
| **`morning_server.py`** | ~200-line stdlib `http.server`, no third-party deps | Serves cockpit page + `state.json`; token-gated `POST /event` → appends `events.jsonl`; `GET /health` (instance nonce); self-owned idle watchdog |
| **`cockpit.html` + `cockpit.js`** | static shell, theme-matched | Polls `state.json` (~1 s); diff-renders; POSTs token-authenticated clicks |
| **`state.json`** | session **writes atomically**, page **reads** | Single source of truth for what is on screen |
| **`events.jsonl`** | server **appends** (stamping a server-side `seq`), session **reads** from a persisted cursor on wake | The click/command queue |
| **background watcher** | `run_in_background` Bash loop, **single live instance** (pid-tracked) | Blocks until `events.jsonl` has lines beyond the session's cursor → exits → auto-wakes session → re-armed **before** processing |

### Two independent half-loops

- **Push (session → page):** session writes `state.json` **atomically** (temp file + `os.replace`) → page poll picks it up → UI updates. Works everywhere; no special mechanism.
- **Pull (page → session):** click → token-authenticated `POST /event` → server stamps `seq` + appends `events.jsonl` → watcher (exit condition: `max(seq) > cursor`) exits → session wakes → re-arms watcher → processes new events → writes new `state.json`. **Mechanism still to be hardened (§2.1).**

**Hard rule: the page is pure I/O, but that is NOT the security boundary.** The page holds no logic and never runs `gh`. However, the auto-wake makes the page (or any forger of `events.jsonl`) a **control channel into a privileged agent** that holds a `gh` token + shell. Security therefore rests on §14's controls (token auth, server-authoritative action lookup, confirm-gating every mutation), not on "the client has no code." See §14 threat model.

## 4. Protocol (the two schemas that decouple page from session)

Both schemas carry `schema_version` (integer). `cockpit.js` is built for a version; on a major mismatch it shows a visible "cockpit out of date — reload" banner rather than mis-rendering (the page and session drift independently behind a browser tab left open across a skill update).

### 4.1 `state.json` (session writes atomically)

```jsonc
{
  "schema_version": 1,
  "state_seq": 42,                 // monotonic; cockpit.js ignores any poll whose state_seq <= last rendered (stale/out-of-order guard)
  "boot_nonce": "a1b2c3",          // per-run nonce; page suppresses ALL action buttons until it sees the current run's nonce (kills stale-prior-run clicks)
  "generated_at": "2026-06-04T09:35:00Z",
  "role": "pm",
  "bridge": "armed",               // armed | processing | resume_needed | server_down | torn_down  (measured, not self-asserted — see §8.3)
  "summary": {
    "greeting": "Good morning",
    "line": "Thu 4 Jun · memory ✓ · 17 open · 2 new overnight · 1 to flag · Ready 4 · 0 overdue",
    "alerts": ["personas repo: 1 stranded edit — flag for MM"],   // Step −1 anti-stranding scan surfaced here, not collapsed into 'memory ✓' (§17)
    "links": [ {"label": "17 open", "url": "https://github.com/.../issues"} ]
  },
  "phases": [
    {"id":"sdr","emoji":"✅","title":"Service Delivery Review","status":"done",
       "sub_phases":[{"id":"recap","title":"Board recap","status":"done"}, {"id":"closure","title":"Closure audit","status":"skipped"}, ...]},
    {"id":"standup","emoji":"🚶","title":"Daily Stand-up","status":"done"},
    {"id":"replenishment","emoji":"🎯","title":"Replenishment","status":"active"},
    {"id":"signals","emoji":"📡","title":"Signals","status":"pending"},
    {"id":"warmup","emoji":"☕","title":"Warm-up","status":"pending"}
    // Friday adds a {"id":"friday_cleanup", ...} phase — phases[] is DAY- and ROLE-driven, not a fixed 5 (§6, §11)
  ],
  "focus": {
    "phase_id": "replenishment",
    "blocks": [ /* see §4.3 block enum */ ],
    "actions": [   // the ACTIVE phase's actions live HERE, not in a flat top-level array (prevents stale-phase clicks)
      {"id":"triage_651","label":"Triage #651 → Backlog","kind":"github_write"},
      {"id":"commit_644","label":"Commit #644 → Ready","kind":"github_write"}
    ]
  }
}
```

- `phases[].status` ∈ `done | active | pending | skipped | error`. A `skipped` phase renders a **dimmed `– skipped` rail row** (NOT removed — removing rows mid-routine shrinks the rail and breaks jump affordance) but **no focus panel**. An `error` phase (builder exited non-zero — e.g. `check_ready_queue.sh` exits 1 on gh-auth/network failure) renders an `error` block with a retry action.
- `focus` is **nullable** (the initial all-pending state before any builder runs has `focus: null`).
- `sub_phases[]` carry SDR's four independently-skippable sub-sections as nested rail rows (matching the agenda's sub-section skip-empty granularity).
- The page renders **only `focus.actions`** for the active phase; the session **rejects any event `id` not in the current `focus.actions`** (defence-in-depth).

### 4.2 `events.jsonl` (server appends, one JSON object per line)

```jsonc
{"seq": 17, "recv_at": 1780558540464, "event_kind": "action", "id": "commit_644"}
{"seq": 18, "recv_at": 1780558550123, "event_kind": "nav",    "id": "next"}
{"seq": 19, "recv_at": 1780558560777, "event_kind": "jump",   "id": "signals"}
```

- **`seq` is server-assigned and strictly increasing** (the server mediates every append). It is the **only** dedup/progress key. The session persists a **cursor** (`last_seq`) in `.morning/cursor` so resume after a session restart/compaction reads from the cursor, never re-processing history. The browser's wall clock is **never** used for dedup (it's non-monotonic); a client `ts`, if present, is advisory display metadata only.
- **`event_kind`** ∈ `action | nav | jump` (the *click type*). This is distinct from an action's **safety class** (`safe | github_write`), which is **never** carried on the wire and **always** re-derived server-side from the action registry by `id` (§14). Renaming avoids the two-`kind` overload.
- **Anchor/link clicks fire NO event** (navigation only). Link navigation is intentionally not recorded; the GitHub view is its own record (documented asymmetry, not a silent gap).
- The watcher exit condition and the session dedup filter are both `seq > cursor`. Readers ignore a trailing partial line (read only up to the last newline).
- **Lifecycle:** `events.jsonl` is append-only within a routine; **truncated and the cursor reset to 0 atomically at run start and at teardown.** A size cap + max-new-lines-per-wake cap bound DoS/wake-amplification (§14).

### 4.3 Block-type enum (closed, versioned) — must cover all five PM phases

`{text, table}` is **not** sufficient. The closed enum, with a coverage check against real phase outputs:

| `type` | Shape | Phase(s) it serves |
|---|---|---|
| `text` | `{text}` | queue line, prose |
| `bullet_list` | `{items:[{text, url?, change?, indent?}]}` | Stand-up blockers/WIP, closure-audit findings list |
| `kv_list` | `{rows:[{label, value, url?, status?}]}` | milestone-health, roadmap-overdue (link + status verb per item) |
| `table` | `{columns, rows}`; a cell is a string **or** `{text?, url?, before?, after?, change?, emphasis?}`; a row carries optional `indent`/`parent_ref` for `↳` sub-issue nesting | the 7-column parent-nested triage diff-table |
| `error` | `{text, detail?}` (paired with a retry action) | any builder that exits non-zero |

The page renders a visible "unsupported block" placeholder for an unknown `type` (never silently drops). A worked example encoding the canonical triage diff-table (parent #171 + sub-issues #177/#178, `before → after` cells, `↳` nesting) is part of the build's protocol test fixture.

### 4.4 Three interaction types

| Type | Render | Behaviour |
|---|---|---|
| **Link** ↗ | underlined accent + `↗`, `rel="noopener noreferrer"`, **https-scheme-allowlisted** (a `javascript:`/`data:` url renders as inert text) | Opens the GitHub source (`target="_blank"`). No event, no session involvement. |
| **Safe action** (`nav`/`jump`/focus) | green button | **Page-local only** — changes which phase is focused. **No state mutation, no `gh`.** May be handled **page-side without waking the session** (§8.2) to slash wake count/cost. |
| **State-mutating action** (`github_write`) | amber button + `⌨` | Fires event → session **surfaces a one-line confirm in chat** naming the exact action + target → on explicit yes, acts once. **Every action that touches GitHub/the board is this tier** — including triage→Backlog and commit→Ready. |

## 5. The cockpit page

Single page, four regions: **summary bar** (pinned; greeting + memory ✓ + the at-a-glance line; counts are links; Step −1 anti-stranding alerts surfaced here), **left rail** (live phase progress incl. `sub_phases` and a Friday-cleanup row when present; doubles as the TodoWrite list; clicking a phase emits a page-local `jump`), **focus panel** (renders `focus.blocks` + `focus.actions`, colour-coded by tier; a legend documents the three interaction types), **footer** (`‹ prev` / `next phase ›` + a **measured** bridge-state indicator).

**Diff-render (specified, not hand-waved):** `cockpit.js` assigns stable `data-id`s to phases, rows, and buttons keyed by their schema `id`; on each poll it compares `state_seq` and per-region content and **re-renders only changed regions**; it **never re-renders the region containing the focused element or an open anchor**; on a JSON parse failure it **retains the last good state and retries next tick** (belt-and-suspenders for the atomic-rename window). "No flicker / stable focus+scroll" is a manual smoke assertion.

## 6. Phase-config model & generalization

**One role-agnostic engine, one per-role config.** A config supplies: the **day-aware ordered phase list** (PM gains a Friday-cleanup phase on Fridays), each phase's **builder** (gathers data → returns structured blocks + rail status + actions, wrapping the routine scripts that already exist — `check_ready_queue.sh`, `check_milestone_health.sh`, the No-Status triage query, the closure-audit scan), the per-phase **actions** (each `safe` or `github_write`), a **per-role block-formatter**, and the **summary-line** composer.

**Role detection:** resolve `os.path.basename(os.path.realpath('.claude/memory'))` against an **explicit allowlist `{pm, scientist, developer}`**. The symlink target (not cwd basename) is used precisely because the developer clone is bare-named `splice-neoepitope-pipeline` — a cwd-basename check would misdetect it. A 4th target `memory_manager` exists but has no morning routine → **out of scope**. An unrecognized basename, dangling/missing symlink, or non-clone cwd ⇒ **abort with a clear "cannot detect role" error**, never default to PM.

See §11 for the **honest** shared/per-role seam (the earlier "News/Status reuse generic builders" claim was overstated).

**Pacing preserved.** Message 1 (greeting + memory ✓ + agenda) is pushed **without arming a blocking watcher and without a stop** (the memory-check line is explicitly *not* a stop point — `feedback_morning_routine_pacing.md`). Phase advance is **non-blocking by default**: the chat narration flows phase→phase as today; the watcher only intercepts an explicit click. An ignore-the-browser user gets the full paced text routine with no hang. Skip-empty applies at the **sub-section** level (SDR's four sub-sections skip independently).

## 7. Concurrency & lifecycle identity

**One server per (clone, role).** Runtime files live in `<clone>/.morning/` (gitignored, `0700`; files `0600`). The reuse/identity decision is keyed on **(clone-path, role, owner-session-token)** jointly — never "a server exists" alone.

**Idempotent-singleton launch — three-part liveness probe (all required to reuse):**
1. `server-info` exists AND records this clone path + **this role**;
2. the recorded PID is alive AND its process command matches `morning_server.py` for this clone (defends against PID reuse — a bare `kill -0` is insufficient);
3. `GET /health` on the recorded port returns the **instance nonce** stored in `server-info` (defends against port reuse by a foreign listener).

Any failure ⇒ treat as stale: atomically remove `server-info`/pidfile, relaunch. The whole detect-or-spawn sequence is serialized under an OS lock (`flock` on `.morning/launch.lock`) so a concurrent same-role launch can't double-spawn; `server-info` is written atomically (temp + rename).

**Two same-role sessions in one clone → refuse-and-warn** (the chosen, fully-specified behaviour; not "attach read-only"). A second launch finds a live server whose **owner token ≠ its own** and aborts: *"a morning routine is already active in this clone; close it or use the other window."* Owner-token equality is what distinguishes *my-own-prior-run* (reuse) from *a live sibling* (contend).

**Stale prior-run state** is never shown: on launch the server resets `state.json` to a `booting` placeholder (with a fresh `boot_nonce`) and truncates `events.jsonl` + resets the cursor, atomically, **before serving**. The page suppresses all action buttons until it sees the current `boot_nonce`.

The shared-server model's only benefit (a unified team view) is deferred to an **optional read-only aggregator** that reads the three per-role `state.json` files. YAGNI until requested.

## 8. The live bidirectional mechanism

### 8.1 Watcher loop — re-arm BEFORE processing (fixes the lost-click race)

Arm-lifetime is **decoupled from event-processing**. The watcher's exit condition is `max(seq in events.jsonl) > cursor`. On wake: the session reads the new lines, **advances the cursor and immediately re-arms a fresh watcher seeded with the new cursor, THEN processes**. Because the cursor (not the watcher's liveness) is the source of truth, a click landing during the multi-second processing window leaves `events.jsonl` longer than the cursor, so the freshly-armed watcher exits immediately and re-wakes — no lost click. **Exactly one watcher is alive at any time** (pid in `.morning/watcher.pid`; a new arm replaces/refuses if a prior pid is still alive). Per-watcher timeout is generous (10–15 min, matching realistic phase dwell) AND re-arming, so a slow human never races the clock.

### 8.2 Wake economy — nav is page-local, builders cache (controls token cost)

- **Pure `nav`/`jump` are handled page-side and do NOT wake the session.** The page changes focus locally and the session reconciles lazily on the next *action* wake. Wakes are reserved for actions that genuinely need session work. This is the primary lever against the §2.1(4) token blow-up.
- **Builders run once when a phase becomes active**; their structured blocks are cached in `state`. A `github_write` action mutates the cached state in-place (or targeted-refetches only the affected row) rather than re-running the whole gh-heavy builder. Live `gh` calls are capped per routine (secondary-rate-limit exposure noted).
- A **token budget** is part of the Phase-0 spike: estimate tokens/wake = (system prompt + current transcript) and model the worst case (N action-wakes × growing transcript); set a hard click-budget and a go/no-go threshold.

### 8.3 Hybrid dispatch & the confirm state machine

On wake the session reads all new events (`seq > cursor`) and processes them in order:
- **`safe`** (re-derived from registry) → execute page-local reconcile, push new state.
- **`github_write`** → **do not execute**; surface a confirm in chat naming the exact `id` + target (#NNN). A pending confirm **freezes further event processing** (subsequent clicks queue at the cursor and process after yes/no) so a human-blocked action can't be skipped or interleaved unsafely. A confirm authorizes **exactly one** registry handler invocation against one named target — never a batch. Cap confirms surfaced per wake at 1; queue the rest.

**`bridge` is measured, not asserted.** The watcher (or server) heartbeats a timestamp into `.morning/`; the server reflects current armed-ness at `GET /health`. `cockpit.js` shows `armed` only if the heartbeat is fresh (< ~2 poll intervals); otherwise the footer flips to **`resume_needed` → "click may have been missed — send any message in chat to resume."** A dropped wake becomes **visible**, not silent — the difference between this being usable and being worse than plain chat.

### 8.4 Graceful fallback if the live loop fails its Phase-0 gate

If the hardened spike (§18 Phase 0) shows the re-arming loop is unreliable (miss-rate over threshold) or too token-expensive, the design **falls back to turn-gated / batch-reconcile**: the page stays a live *view* (push works everywhere) and queues clicks; the session reconciles them on the user's next natural chat turn ("click, then go"). The cockpit is still valuable as a glanceable dashboard; only the *per-click auto-wake* is dropped. This fallback is a first-class outcome, not a failure.

## 9. Skill orchestration flow

```
/morning  (or "good morning" → instruction rule invokes the skill — §2.2)
  1. detect role (memory symlink allowlist; abort if unknown)        → load role config
  2. launch-or-reuse server (flock; 3-part liveness; reset state+events; fresh boot_nonce)
       → print URL; open VSCode Simple Browser (fallback: system browser)
  3. Step −1 checks: memory check + personas-repo anti-stranding scan → surface results in summary.alerts (NOT a stop point)
  4. write Message-1 state (summary + rail, all phases pending, focus:null) — NO watcher armed, NO stop
  5. for each phase (day/role-driven; skip-empty at sub-section level):
       a. mark active; run builder (cache result) → focus blocks + actions
       b. push state; echo the phase in chat (durable record, existing formatting rules)
       c. arm watcher (re-arm model §8.1)
       d. on wake: re-arm → safe: reconcile+repush;  github_write: confirm-in-chat → act+repush;  (nav/jump never wake)
       e. mark done; advance (non-blocking — chat flows even with no click)
  6. Warm-up: suggest one XS/S role-queue item
  7. teardown in a finally-wrapper (runs even on mid-routine exception): stop server, disarm watcher, truncate events
```

The chat narration per phase is unchanged from today's routine (same formatting, diff tables, link+keyword convention). The cockpit is **additive**.

## 10. Lifecycle & edge cases

- **Server self-owns its lifecycle** (does not depend on the skill being alive): an **idle watchdog thread inside `morning_server.py`** exits the process after N minutes with no poll/POST/state-write — so an abandoned session is cleaned up even if the skill never reaches teardown. **Idle resets on cockpit POLLS** (a polling page = a watching user), so a legitimate >30-min paced pause does not tear down a live server; abandonment = closed browser = no polls.
- **Server died mid-routine** is detected actively, not push-only: the watcher periodically pings `GET /health`; on death it wakes the session with a synthetic `server_down` event → relaunch (preferring the **same** recorded port so existing tabs keep working) + repush. Independently, `cockpit.js` flips the footer to `server_down` / `reconnecting` on POST/poll failure and **queues the click locally for replay**, so the user is never clicking into a silent void.
- **Browser closed** → session continues; reopening the URL re-renders from current `state.json`. If a relaunch changed the port, the session prints the new URL.
- **Teardown** is idempotent and crash-safe (finally-wrapper + the server's own watchdog backstop).
- **Two same-role sessions / same clone** → refuse-and-warn (§7).
- **Watcher died / wake missed** → `bridge:resume_needed` surfaced on the page; "send any message to resume" re-reads from the persisted cursor (no double-processing) and relaunches the server if dead.

## 11. Generalization to Dev/Sci — the honest scope of D6

**Genuinely role-agnostic (write once, never change):** the engine, server, protocol, dispatch, lifecycle, the watcher loop, and the **Warm-up** builder (one XS/S from `role:<role>` — a true scope swap).

**Shared *plumbing*, per-role *formatting/side-effects*:** News/Status share data-gathering (Zotero+Issues dedup, `poll_releases.py`, the issue-creation cap, the standup-hygiene two-halves, the open-`role:<role>`-work survey) — but their **output formats diverge** and the renderer must too:
- PM Replenishment → `before → after` **diff-table**; PM Stand-up adds a **WIP-awareness** flag + **cross-role** walk; PM Signals carries a **landscape-doc maintenance** side-effect (`research/multi_agent_landscape.md`).
- Scientist News → **numbered-footnote citation** prose + sources list.
- Developer News → one-line **signal-tagged bullets** (4-value taxonomy).

So a "second role" is **config + a per-role block-formatter**, and Sci/Dev will require **new renderer block types** in `cockpit.js` (`citation_list`, `signal_bullet`) — an engine touch the original "config-only" framing denied. **D6 is corrected to: the orchestration core is role-agnostic; presentation is per-role.** PM-specific logic (WIP flag, cross-role walk, landscape hook) lives in `pm.py` as config-supplied extension blocks the generic builder appends — not in the engine.

**D6 acceptance test (strengthened):** exercise a **real** Sci or Dev builder end-to-end through `cockpit.js` rendering (or a fixture emitting the citation-list / signal-bullet shapes). The earlier "stubbed builders emitting `text` blocks" dry-run **cannot** catch this gap — it stubs out exactly the non-generic part. If the test requires new block types, that surfacing **is** the corrected-D6 boundary, by design.

## 12. Testing strategy

Python lives under **`tools/morning/`** (mirroring `tools/ci/`, `tools/news/`), tests run via **`pytest tools/morning/ -v` on Python 3.11**, wired into the existing **`ci-tools-pytest`** CI job (`.github/workflows/tests.yml`). No new Python env is invented (CLAUDE.md enumerates exactly 4; this reuses the CI-tools convention). Only `SKILL.md` + cockpit assets live under `.claude/skills/morning/`.

- **Unit — `morning_server.py` (headless):** serves `state.json`, token-gates + validates `POST /event`, stamps `seq`, appends well-formed `events.jsonl`, atomic state writes, `GET /health` nonce, idle watchdog, free-port selection, `0600` perms, log redaction.
- **Protocol round-trip:** write state → fetch → POST event → assert appended with monotonic `seq`; cursor-based re-read idempotence; a forged event carrying `kind:"safe"` for a registered `github_write` is **still confirm-gated**; a `javascript:` url renders inert; an unknown/wrong-phase `id` no-ops.
- **Engine `--dry-run`:** fixture role-config, stubbed builders, no `gh` → asserts the phase loop, state transitions, skip-empty (incl. sub-phase), error-phase rendering, and dispatch. **CI-safe.**
- **D6 acceptance:** real Sci/Dev block shapes through the renderer (§11).
- **Phase-0 live-loop spike (manual, gated):** ≥5 consecutive re-arm wakes on the **purpose-built** server; measure **miss-rate** and **tokens/wake**; single-live-watcher + no-orphan + clean-teardown assertions; **go/no-go thresholds** (e.g. abort live-bridge if miss-rate > 5%). Documented as a smoke step (needs harness + browser), not CI.

## 13. File layout

```
.claude/skills/morning/
  SKILL.md                 # entry: triggers, role detection, flow
  cockpit.html             # static shell
  cockpit.js               # poller + diff-renderer + token POST
  cockpit.css
tools/morning/
  engine.py                # role-agnostic orchestration + dispatch
  server.py                # morning_server.py (stdlib http.server)
  roles/{pm,developer,scientist}.py   # day-aware phase lists + builders + actions + formatter
  builders/                # shared plumbing + per-role builders (wrap existing routine scripts)
  tests/                   # unit + protocol + dry-run + D6 acceptance
  README.md                # interpreter (CI Python 3.11), smoke-test steps
<clone>/.morning/          # gitignored runtime (0700; files 0600): state.json, events.jsonl, cursor,
                           # server-info, server.pid, watcher.pid, launch.lock, server.log, token
```

`.gitignore`: add `.morning/` under its **own** comment (not bundled with the `.superpowers/` brainstorming-mockups rationale). Ensure the non-clone fallback runtime location (§15) is always under an ignored path or a session tmp dir.

## 14. Security & safety — threat model

**Actor:** a malicious local web page (CSRF), another same-host process, or a forged `events.jsonl`. **Asset:** the privileged agent (holds a `gh` token + project-board write + a shell) reachable via the auto-wake. **The page-is-pure-I/O property is NOT a control** — the auto-wake is precisely a control channel into that agent.

Controls (all required):

1. **Per-launch secret token** generated at server start, written only into `cockpit.html` + `.morning/server-info` (never logged, passed via **header** never URL). Every `POST /event` must present it → `403` otherwise. Loopback binding is necessary but **not sufficient** against CSRF (a cross-origin simple POST is *sent* even if the response can't be read — and the side effect, not the response, is the attack).
2. **Origin/Referer check** on `POST /event` (reject cross-origin); reject any body whose content-type isn't `application/json`. Defence-in-depth with the token.
3. **`safe` redefined narrowly = page-local navigation only** (no state mutation, no `gh`). **Every GitHub/board mutation is `github_write` (confirm-gated)** — including triage→Backlog. A dispatcher assertion: a handler tagged `safe` that attempts any subprocess/`gh` call is a **hard error**, not silent execution.
4. **Server-authoritative action kind.** The dispatcher **ignores any `kind` on the wire** and re-derives safety class solely by `id` → action registry. Unknown/absent/wrong-phase `id` → no-op. (Safety class is not carried on the wire at all.)
5. **`events.jsonl` is UNTRUSTED, typed command input** — never prose for the LLM. Parse each line, match `id` against the registry, act only on the registry handler, and **never surface raw event-supplied strings into reasoning or chat** (re-derive every label from trusted `state.json`/registry). A **per-wake event-count cap + rate limit** bounds wake-amplification (token-cost DoS). A malformed line advances the cursor without crashing.
6. **File perms + log redaction.** `.morning/` is `0700`, files `0600` (correcting the earlier "world-readable is fine" note). `morning_server.py` overrides `log_message` to record method + path-without-query + status only — **no bodies, no query strings, no token**. The "no secrets" invariant covers `state.json`, `events.jsonl`, **and** `server.log`.
7. **Link hardening.** `cockpit.js` allowlists `https` schemes (rejects `javascript:`/`data:`/`file:` → inert text) and adds `rel="noopener noreferrer"` to every `target="_blank"`. Builders only emit URLs constructed from known GitHub patterns, never echo a user-supplied URL field.
8. **Idle-wake window.** The auto-wake channel is open while the server lives; the server's idle watchdog (§10) bounds an abandoned-but-live window. On idempotent reuse, the token is **regenerated** so a stale/leaked token can't be replayed.

Boundary after mitigations: **same-user processes that can read `server-info`/`cockpit.html` are trusted; cross-user and browser-CSRF are not.**

## 15. Open questions (resolve in writing-plans or Phase 0)

- VSCode Simple Browser `target="_blank"` behaviour (system browser vs new tab) — confirm during build; spec assumes system browser (desirable for GitHub).
- Watcher implementation (`fswatch` if present vs portable poll loop) — poll loop is the portable default.
- Realistic per-click wake latency + miss-rate + tokens/wake — **measured in Phase 0**; these decide §8.4 (live vs fallback).
- Non-clone-cwd runtime dir: default `$PWD/.morning/`; must always be under an ignored path or a session tmp dir.
- Supported-environment matrix for the auto-wake (harness/IDE versions); the skill runs a **startup self-check** (arm watcher → POST synthetic event → confirm round-trip within N s) before telling the user the bridge is live; on failure it falls back to read-only/chat-only mode rather than presenting dead buttons.

## 16. Out of scope (v1)

Unified multi-role aggregator view; SSE/WebSocket push (polling suffices); PDF/handout export; mobile/responsive beyond desktop-beside-IDE; editing issue *content* from the page (only status/triage/commit actions + links).

## 17. Memory / governance follow-ups

- **Tracking Issue (do first):** file a role-labelled, sized Issue with a Priority rationale; add its `[Issue #N](url)` to this spec's header. This is a **propose-and-confirm** item (cross-role framework + ongoing maintenance) per `feedback_ask_for_help.md`, not an auto-committed build.
- **"good morning" routing edit (personas repo, MM-committed):** the Always-in-effect rule must be amended to "invoke the `/morning` skill," and **must land before the skill is usable** (§2.2). Flag for MM; this role does not commit the personas repo.
- **Lab-notebook merge gate:** the build PR ships via `scripts/audit_and_merge.sh` (not bare `gh pr merge`); a `research/lab_notebook/pm.md` entry (or `<!-- skip-lab-notebook: routine -->`) is required at merge.
- The skill code itself lives in the **project repo** `.claude/` + `tools/morning/` and is committed normally.

## 18. Suggested build phasing (detail → writing-plans)

0. **Hardened live-loop spike (GATE) — ✅ DONE 2026-06-04, PASSED.** On the purpose-built server: **6/6 consecutive re-arm wakes, 0% miss, 9/9 events consumed, ~7–35 s latency, coalescing confirmed** (see §2.1). Go/no-go ⇒ **GO on full live bridge**. Build Phase 0 remaining work: re-verify with **real browser** clicks at routine cadence (vs the spike's `curl` driver) and add the single-live-watcher / no-orphan / clean-teardown assertions as the productionised smoke test.
1. **Server + protocol** — `server.py` + schemas (`seq`, cursor, versioning, atomic writes, token, `/health`, idle watchdog, log redaction) + unit/round-trip tests, wired into `ci-tools-pytest`.
2. **Cockpit page** — `cockpit.html/js/css` against a fixture `state.json`; diff-render, link allowlist, bridge-state footer, token POST.
3. **Engine + watcher loop** — phase loop, re-arm-before-process, dispatch, confirm state machine, `--dry-run`.
4. **PM role config** — wrap real routine scripts as cached builders; wire actions (all mutations `github_write`); day-aware phases incl. Friday cleanup + sub-phases; productionise the smoke loop.
5. **SKILL.md + routing** — role-detection allowlist, idempotent launch, crash-safe teardown, startup self-check.
6. **Dev/Sci configs** — prove **corrected** D6 (per-role formatter + any new block types) via the strengthened acceptance test.

---

### Appendix — revision provenance

This spec was revised after a 7-lens adversarial review (architecture, protocol, concurrency, generalization, consistency, security, feasibility; 2026-06-04). Material corrections folded in: spike demoted to partially-de-risked with a Phase-0 go/no-go gate (feasibility); CSRF token + narrow-`safe` + registry-authoritative kind + untrusted-events model (security); re-arm-before-process lost-click fix + atomic writes + diff-render model (architecture); `seq`-based dedup + cursor persistence + versioning + block-type enum + error/skipped states (protocol); health-nonce liveness + crash-safe teardown + idle-on-polls + refuse-and-warn (concurrency); honest D6 scope + role-detection allowlist (generalization); `tools/morning/` + CI wiring + tracking Issue + day-aware phases (consistency).
