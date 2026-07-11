Answer one question: **is anyone waiting on me?** Scan the two team-coordination channels, report what each one covered, and surface anything addressed to or relevant to your role. Confirm with the user before acting on any request.

This is the single definition of the coordination scan. The morning routine (Beat 2) and the resume routine both delegate here - do not restate these mechanics in memory files, or the copies drift (they already did once: see the last section).

## 1. Board pings - board-wide, NOT role-scoped

A comment addressed `**To:** <your role>` can land on **any** Issue regardless of its `role:*` label (a Developer replies `To: PM` on a `role:developer` Issue). Scanning only your own `role:*` Issues therefore **structurally misses every cross-role ping**.

```bash
python3 scripts/pm/scan_addressed_comments.py --role <pm|scientist|developer|memory_manager>
```

It reads the session watermark itself (1-day overlap, 7-day floor when the marker is absent), scans board-wide, and prints pings grouped by Issue with author, timestamp, and snippet. A named `To:` addressee **acks on arrival**.

Do **not** substitute `board_open_items.py --role <role>` here. That is the role-scoped scan this command used until [Issue #1114](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1114); it silently drops cross-role pings.

## 2. Open Discussions - Team Coordination category

Open means it still needs attention.

```bash
gh api graphql -f query='{ repository(owner:"Jin-HoMLee", name:"splice-neoepitope-pipeline"){ discussions(first:30, categoryId:"DIC_kwDORwn9EM4C-Jo6", states:[OPEN]){ totalCount nodes{ number title author{login} } } } }'
```

<!-- first:30 is unpaginated - fine at 4-role scale. Check totalCount; add a hasNextPage/after cursor loop if open threads ever exceed 30. -->

## 3. Report coverage on every run - including when empty

State what was actually scanned, so each run re-teaches its own scope:

```
Inbox: board pings (board-wide, since <floor>): N - open Discussions: N
```

An empty inbox is a result worth printing. Silence is indistinguishable from a scan that never ran.

## Then

Bring any pending item to the user, summarize what it asks, and act only after confirmation. Close a Discussion / board thread once its ask is satisfied (open = live). Full protocol: `shared/feedback_team_coordination.md`.

## Why the board-wide scan is load-bearing

`scan_addressed_comments.py` ([Issue #901](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/901)) landed board-wide precisely because the role-scoped scan had already lost a real ping: PM missed a Developer's `To: PM` reply on `role:developer` [Issue #887](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/887) (2026-06-29). The morning routine adopted the helper; this command was left on the old scan for three weeks and kept losing cross-role pings until [Issue #1114](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1114). One scan, one definition - that is what keeps this from happening a third time.
