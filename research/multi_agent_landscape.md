# Multi-agent / agentic-coding orchestration — landscape

Curated reference of the multi-agent / agentic-coding field with one-line summaries, our rationale per item, and links to internal Issues and framing memories. Doubles as a portfolio artifact (we surveyed the field; we're not cargo-culting) and a project memory aid (where prior reasoning on Symphony, Fugu, Managed Agents, etc. lives without re-grepping `news_log.md`).

Maintained by PM. Updated when `news_log.md` entries get tagged `methodology-signal` — see [Maintenance](#maintenance) at the bottom.

For our position relative to this landscape — why PM/Sci/Dev, why "multi-role" instead of "multi-agent" — see [Our position](#our-position).

---

## Frameworks

Concrete products / tooling that could plausibly displace, augment, or be borrowed from for our PM/Sci/Dev workflow.

### Anthropic Managed Agents — Dreaming + Multiagent Orchestration

- **Source:** [Anthropic blog, 2026-05-06](https://claude.com/blog/new-in-claude-managed-agents)
- **Summary:** Agents review past sessions to self-improve ("dreaming"); a lead agent delegates to specialists rather than the human dispatching each role.
- **Our rationale:** **Observe.** Dreaming is the closest production analog to our `/cerebrum` memory-rehydration pattern. Orchestrator-delegation is what PM does for Sci/Dev today, *but in our setup the human opens the next session* — we don't have agent-to-agent triggers (per [multi-role, not multi-agent](#our-position)). Tracking via [Issue #305](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/305) (PM eval + memory tightening borrow).

### OpenAI Symphony

- **Source:** Reported late April 2026 (no public docs yet — pattern-only)
- **Summary:** Enterprise multi-agent orchestration; uses PM-style project boards as the control plane for coding agents.
- **Our rationale:** **Pattern-confirmation, no action.** The "boards-as-control-plane" framing matches what we do — GitHub Projects v2 IS our PM/Sci/Dev coordination surface (Status / Priority / Size / Target Date fields drive triage). Product itself targets corporate teams (not solo/research), so no adoption path. Confirms the convergent pattern.

### Sakana Fugu

- **Source:** [Sakana 2026 (beta)](https://sakana.ai/fugu-beta/)
- **Summary:** Multi-agent orchestration as a foundation model — dynamically assembles agent teams from a pool, no fixed roles.
- **Our rationale:** **Explicit counter-position.** Opposite of our fixed PM/Sci/Dev split — per [Role decomposition is domain-bespoke](#our-position), our roster is intentionally bespoke to bioinformatics research (Scientist exists because literature/biology reasoning has no analog in pure coding). Fugu optimises for generality; we optimise for domain fit. Still borrowed three secondary mechanics: [Issue #324](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/324) (model routing per task type), [Issue #325](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325) (post-merge critic), [Issue #326](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/326) (memory consolidation).

### GitHub Copilot desktop app

- **Source:** [GitHub changelog 2026-05-14](https://github.blog/changelog/2026-05-14-github-copilot-app-is-now-available-in-technical-preview/) (technical preview)
- **Summary:** GitHub-native desktop client to start agentic dev sessions from an issue/PR/prompt; isolates work in-flight, supports mid-session steering, lands via PR review.
- **Our rationale:** **Pattern-confirmation, no action.** Same parallel-session UI shape as Claude Code Agent View (news_log 2026-05-13) and convergent with our PM/Sci/Dev orchestration — but GitHub-anchored rather than terminal/IDE-anchored. Fits the "boards-as-control-plane" framing (cf. [OpenAI Symphony](#openai-symphony)). Tech preview only; too early for adoption decision but worth watching for GA + agent-callable surface (per [UI features ≠ agent capabilities](#our-position)). No action.

---

## Methodology / framing

Essays and reports that shape how the field talks about itself. Distinct from frameworks above — these change vocabulary and mental model, not tooling.

### Anthropic 2026 Agentic Coding Trends Report

- **Source:** Anthropic 2026 (PDF; tracked via [Issue #235](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/235))
- **Summary:** Industry data on multi-agent coding workflows. Proposes canonical 4-specialist split: Architecture & Design / Implementation & Coding / Testing & Validation / Review & Docs.
- **Our rationale:** **Pattern borrow + roster reject.** Convergent pattern (one orchestrator + multiple specialists, each with own context) applies. Roster doesn't — that report is scoped to *coding*; we're scoped to bioinformatics research, where literature integration and biological interpretation have no analog in pure coding (so Scientist instead of Designer). Primary source for both [Role decomposition is domain-bespoke](#our-position) and [multi-role, not multi-agent](#our-position). Skim + 5 PM-practices cross-checks in [Issue #235](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/235).

### Mason — "Coherence Through Orchestration, Not Autonomy"

- **Source:** [Mike Mason, Jan 2026](https://mikemason.ca/writing/ai-coding-agents-jan-2026/)
- **Summary:** Argues "multi-agent" coding works because of orchestration discipline, not agent autonomy.
- **Our rationale:** **Direct reinforcement.** Primary external source for [multi-role, not multi-agent](#our-position) — we don't claim autonomy, we claim disciplined orchestration. The user is the message bus; agents communicate via human-readable artifacts (standup, sub-issues). No action — already internalised.

### AddyOsmani — "The Code Agent Orchestra"

- **Source:** [AddyOsmani blog, 2026-05-12](https://addyosmani.com/blog/code-agent-orchestra/)
- **Summary:** Essay on what makes multi-agent coding work — section roles, conductor metaphor.
- **Our rationale:** **Convergent vocabulary.** Maps cleanly onto our 3-layer orchestration: human (conductor) → PM (section leader) → Sci/Dev (instrumentalists). No new action items — this landscape doc itself ([Issue #337](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/337)) is the response to the essay.

### digitalapplied — "Pattern Language 2026"

- **Source:** [digitalapplied, 2026-05](https://www.digitalapplied.com/blog/multi-agent-orchestration-patterns-producer-consumer)
- **Summary:** Formal taxonomy reducing every reliable multi-agent system to 5 archetypes — producer (decomposes ambiguity into work items), consumer (executes them), coordinator (routes + bounded fan-out), critic (suggestions, no gate authority), judge (binary go/no-go). 12 composition rules (acyclicity, idempotent consumers, bounded fan-out, etc.); 8 predictable failure modes (cycle formation, critic-judge deadlock, coordinator overload, silent drift, others). Maps cleanly to LangGraph, CrewAI, OpenAI Agents SDK, Claude Agent SDK without changing the underlying shape.
- **Our rationale:** **Vocabulary candidate — open question.** Distinct from AddyOsmani above (essay vs. taxonomy) and from the Trends Report (coding-scoped vs. domain-agnostic). Possible re-description of our shape: `human = meta-coordinator → PM = coordinator + judge → Scientist = producer + critic + Developer = consumer + producer (infra) + critic`. The 8-failure-mode list doubles as a portfolio audit lens (have we hit critic-judge deadlock? silent drift?). Decision tracked via [Issue #353](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/353) (PM vocabulary adoption).

---

## Our position

We describe this project as a **human-orchestrated multi-role workflow**, not a "multi-agent system." Three internal framing notes (PM-shared memory files; reproduced inline here so this doc is self-contained) pin this down:

1. **"multi-role, not multi-agent."** Industry "multi-agent" usage implies agent-to-agent autonomy — fire-and-forget orchestration, self-organising task assignment. What we have is role-bound Claude Code sessions (PM / Sci / Dev) where *the human is the message bus*: when one role flags work for another, the human opens the other session and surfaces the flag. Agents never trigger each other. Calling this "multi-agent" overclaims; "multi-role workflow" is honest.

2. **Role decomposition is domain-bespoke.** The *pattern* (one orchestrator + multiple specialists, each with own context window) is convergent across the field. The *roster* (PM/Sci/Dev) is bespoke to bioinformatics research — literature integration, biological interpretation, hypothesis design, and results analysis have no analog in pure coding work, so we factor that out as a Scientist role. Other domains would warrant different specialists: UI product → Designer; deployed service → SRE; ML research → Modeller/Data specialist. Future role additions (Cloud Eng, Writer) get justified by domain need, not by analogy to other projects' rosters. *External validation:* Gartner ([arXiv 2601.13671](https://arxiv.org/html/2601.13671v1)) forecasts >40% agentic AI projects cancelled by 2027 at the orchestration ↔ domain-logic junction — supports domain-bespoke roster over generic 4-specialist borrow.

3. **Cerebrum vs project scope.** Cerebrum is the multi-agent meta-framework that lives *above* all projects and roles — it defines the role pattern, memory format conventions, role-binding mechanics. This project (splice-neoepitope-pipeline) is one *instance* using Cerebrum. Frameworks evaluated in this landscape doc are weighed against *this project's* needs (bioinformatics-research workflow), not Cerebrum's meta-level concerns. Cross-project generality is a separate question.

**3-layer orchestration shape:**

```
human  (conductor / high-level orchestrator — decides what to build, who to dispatch)
  └─ PM  (section leader / low-level orchestrator — triages, scopes Issues, runs morning routine)
       ├─ Scientist  (literature, biology, hypothesis, results interpretation)
       └─ Developer  (implementation, testing, CI, infrastructure)
```

The standup file (`pm/shared/team_standup.md`), the GitHub Projects v2 board, and the per-role Cerebrum memory directories are the durable substrates that make role-bound sessions coherent across time and across context-window resets.

---

## Maintenance

PM-owned. Update triggers:

- News_log entry tagged `methodology-signal` lands → PM scans during morning routine, decides whether to add an entry here (or whether the item belongs only in news_log).
- New framing memory in `pm/shared/` cross-references this doc → update [Our position](#our-position) to mention it.
- New Issue concretely borrows from or extends a framework listed here → backfill the framework's rationale block with the Issue link.
- Framework deprecated or pivoted → don't delete the entry; update the rationale with the deprecation note (the doc is a journal of *which* options we considered, not just current options).

Last sweep: 2026-05-15 (GitHub Copilot desktop app backfilled to Frameworks as pattern-confirmation).
Prior: 2026-05-14 (Gartner cancellation forecast backfilled to `our_position` #2 as external validation).
Prior: 2026-05-12 (initial scaffold; 6 backfilled items: Symphony, Trends Report, Managed Agents, Mason, Fugu, Code Agent Orchestra).
