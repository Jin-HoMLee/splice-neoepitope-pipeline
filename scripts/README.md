# `scripts/`

Operational helpers for the pipeline and the PM board routine. Pipeline setup /
run scripts (`setup_*.sh`, `run_cloud_gpu.sh`, `prepare_*_data.sh`, etc.) are
documented at their call sites (`docs/`, `CLAUDE.md`); this README covers the
**board-health check family** and the convention that governs it.

## The `check_*_health` family

Three sibling checks back the PM morning-routine Service Delivery Review. Each
follows the same exit-code contract:

| Script | Watches | Routine phase |
|--------|---------|---------------|
| `check_milestone_health.sh` | open milestone `due_on` clock (overdue / imminent) | milestone health |
| `check_ready_queue.sh`      | Ready-queue depth (demand-aware, per-role) | replenishment trigger |
| `check_roadmap_health.py`   | per-issue board **Target date** (overdue) | roadmap-overdue sweep |

**Exit-code contract (shared):** `0` all-clear · `2` at least one finding · `1`
usage / runtime error. Each takes `-h/--help` and a header docstring.

`check_roadmap_health.py` is **not** redundant with `check_milestone_health.sh`:
the milestone clock watches `due_on`, while the Target date is a per-issue field
that can **desync** from `due_on` (a milestone change not paired with a Target
re-sync, or an un-milestoned issue still carrying a stale Target). The roadmap
check is the detection backstop for that desync (Issue #704).

## Convention: language-by-data-source

The `check_*` family splits cleanly by **data source**, and the split is
load-bearing — pick the language to match:

- **Plain GitHub REST objects** (milestones, issues, counts) → **shell + `jq`**.
  Date math via `fromdateiso8601 - now`; no pagination. → the `.sh` checks.
- **ProjectV2 item fields** (Status / Priority / Size / **Target date**),
  paginated across the whole board (~700+ items, sorted Done-first so a
  single-page query silently truncates open work) → **Python**, reusing
  `board_open_items.py`'s paginator + `fieldValues`-union unpacking
  (`... on ProjectV2ItemFieldDateValue`) rather than re-hand-rolling the cursor
  loop in `jq`, where the union unpacking gets ugly. → `board_open_items.py`,
  `check_roadmap_health.py`.

New board-health checks: if they read a ProjectV2 field, add them on the Python
side and import from `board_open_items` (`fetch_all_items` + `normalize`); if
they read only REST objects, a `.sh` + `jq` script is the lighter fit.
