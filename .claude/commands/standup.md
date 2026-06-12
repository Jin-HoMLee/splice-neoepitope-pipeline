Scan the two team-coordination channels and surface anything addressed to or relevant to your role; confirm with the user before acting on any request.

1. **Board pings** — comments `**To:** <your role>` on your `role:*` Issues + recent comment activity: `python3 scripts/board_open_items.py --role <role> --status "In progress"` (then `"Ready for review"`, `"In review"`).
2. **Open Discussions** — Team Coordination category (open = needs attention):
   `gh api graphql -f query='{ repository(owner:"Jin-HoMLee", name:"splice-neoepitope-pipeline"){ discussions(first:30, categoryId:"DIC_kwDORwn9EM4C-Jo6", states:[OPEN]){ totalCount nodes{ number title author{login} } } } }'`
   <!-- NOTE: first:30 is unpaginated — fine at 4-role scale, but add a hasNextPage/after cursor loop if open threads ever exceed 30 (totalCount surfaces the count to check). -->

Bring any pending item to the user, summarize what it asks, and act only after confirmation. Close a Discussion / board thread when its ask is satisfied (open = live). Full protocol: `shared/feedback_team_coordination.md`.
