# Multi-agent / agentic-coding orchestration — landscape

Curated reference of the multi-agent / agentic-coding field with one-line summaries, our rationale per item, and links to internal Issues and framing memories. Doubles as a portfolio artifact (we surveyed the field; we're not cargo-culting) and a project memory aid (where prior reasoning on Symphony, Fugu, Managed Agents, etc. lives without re-grepping conversation history).

Maintained by PM. Updated when a morning briefing surfaces a `methodology-signal` item — see [Maintenance](#maintenance) at the bottom.

For our position relative to this landscape — why PM/Sci/Dev, why "multi-role" instead of "multi-agent" — see [Our position](#our-position).

**Companion doc:** this file surveys the *frameworks/framing* (what exists in the field). For the *governance-evidence* survey — what the 2024-2026 literature says about **how to run** an agent/role team, which human-borrowed Kanban dials are miscalibrated, and the MAST failure taxonomy — see [`agent_team_governance_research_2026-07.md`](agent_team_governance_research_2026-07.md) (PM deep-research report, 2026-07).

---

## Entry status: the funnel

Every entry carries a **`status:`** ring: a maturity signal borrowed from the ThoughtWorks Technology Radar, adapted so this curated landscape doubles as a staging -> commit funnel.

| Ring | Meaning | Terminal? |
|------|---------|-----------|
| `Assess` | Watching ("too early" / "Observe"); carries a one-line **DoR** naming the event that would promote it to `Trial`. | no |
| `Trial` | Eval candidate; has (or gets, on promotion) an open board #9 eval Issue. | no (-> `#N`) |
| `Adopt` | Adopted into the pipeline or our practice. | yes (-> `#N`) |
| `Hold` | Counter-position: evaluated, actively avoid, never start. | yes |
| `Reference` | Evaluated; no product adoption, but the idea or pattern informed our thinking or confirmed our approach. | yes |
| `Rejected` | Evaluated a plausible candidate and walked away; records *why*. | yes |

We borrow the Radar's ring taxonomy but drop its calendar cadence.
Entries re-ring on a *named event* - a promotion trigger firing, or an eval Issue closing - never on a bi-annual Radar date.

**Promotion is the `Assess -> Trial` flip.**
The same edit that flips the ring opens the board #9 eval Issue (role label + priority rationale + per-tool Quarto deck, per the eval-Issue convention), cross-linked both ways.
Terminal flips (`Adopt` / `Hold` / `Reference` / `Rejected`) are recorded back here when the eval Issue closes, via the close-time outcome-routing ritual.

> **Ring-set note (2026-07-02, [#553](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/553)):** `Reference` is a 6th ring added during the retrofit, beyond #553's original 5-ring design, and **confirmed by PM**.
> About half the existing entries are "pattern-confirmation / reinforcement" - evaluated, no product adopted, but the idea confirmed our approach or lent vocabulary - which none of `Assess/Trial/Adopt/Hold/Rejected` captured honestly (not watched, not to-avoid, and not "Rejected" when we endorse the idea).
> The distinction it draws: `Rejected` = walked away, not useful; `Reference` = walked away from the product, kept the idea.

---

## Frameworks

Concrete products / tooling that could plausibly displace, augment, or be borrowed from for our PM/Sci/Dev workflow.

### Anthropic Managed Agents — Multiagent Orchestration + Dreaming + Outcomes

- **status:** Reference
- **Source:** [Anthropic blog, 2026-05-06](https://claude.com/blog/new-in-claude-managed-agents); [9to5Mac coverage, 2026-05-07](https://9to5mac.com/2026/05/07/anthropic-updates-claude-managed-agents-with-three-new-features/).
- **Summary:** Three features. **Multiagent Orchestration:** a lead agent delegates to autonomous specialist agents on a shared filesystem; specialists run in parallel, lead checks back mid-workflow, persistent events so every agent remembers what it has done. **Dreaming:** scheduled review of past sessions to extract patterns and curate memory across agents — auto-update by default, optional human review gate. **Outcomes:** users write success rubrics; agents iterate autonomously toward them, with an independent grader that "isn't influenced by the agent's reasoning"; webhook on completion.
- **Our rationale:** **Observe — per-feature mapping.** **Multiagent Orchestration** maps directly onto the autonomous pattern that [multi-role, not multi-agent](#our-position) explicitly does *not* claim — in our setup the human opens the next session and there is no agent-to-agent trigger. The model vendor itself using "multi-agent" in this autonomous sense validates the framing distinction. **Dreaming** is the closest production analog to our `/memory` memory-rehydration pattern, but inverted: their default is auto-update with optional review; ours is always human-initiated and reviewed at write time. **Outcomes** has no direct analog in our setup (we don't yet codify success rubrics for autonomous iteration; closest precursor is the GitHub-Actions `Claude Code` review trigger, which has a human-set criterion but no agent-driven re-iteration). Closed [Issue #305](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/305) with landscape-tightening (the concrete contrast belongs here; the framing rule in `shared/feedback_multi_role_not_multi_agent.md` stays unchanged).

### OpenAI Symphony

- **status:** Reference
- **Source:** Reported late April 2026 (no public docs yet — pattern-only)
- **Summary:** Enterprise multi-agent orchestration; uses PM-style project boards as the control plane for coding agents.
- **Our rationale:** **Pattern-confirmation, no action.** The "boards-as-control-plane" framing matches what we do — GitHub Projects v2 IS our PM/Sci/Dev coordination surface (Status / Priority / Size / Target Date fields drive triage). Product itself targets corporate teams (not solo/research), so no adoption path. Confirms the convergent pattern.

### Sakana Fugu

- **status:** Hold (roster counter-position; the 3 borrowed secondary mechanics are tracked separately as Adopt at #324 / #325 / #326)
- **Source:** [Sakana 2026 (beta)](https://sakana.ai/fugu-beta/)
- **Summary:** Multi-agent orchestration as a foundation model — dynamically assembles agent teams from a pool, no fixed roles.
- **Our rationale:** **Explicit counter-position.** Opposite of our fixed PM/Sci/Dev split — per [Role decomposition is domain-bespoke](#our-position), our roster is intentionally bespoke to bioinformatics research (Scientist exists because literature/biology reasoning has no analog in pure coding). Fugu optimises for generality; we optimise for domain fit. Still borrowed three secondary mechanics: [Issue #324](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/324) (model routing per task type), [Issue #325](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325) (post-merge critic), [Issue #326](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/326) (memory consolidation).

### GitHub Copilot desktop app

- **status:** Assess
- **DoR (-> Trial):** GitHub ships GA with an *agent-callable* surface (API / CLI), not just the desktop UI (per [UI features ≠ agent capabilities](#our-position)).
- **Source:** [GitHub changelog 2026-05-14](https://github.blog/changelog/2026-05-14-github-copilot-app-is-now-available-in-technical-preview/) (technical preview)
- **Summary:** GitHub-native desktop client to start agentic dev sessions from an issue/PR/prompt; isolates work in-flight, supports mid-session steering, lands via PR review.
- **Our rationale:** **Pattern-confirmation, no action.** Same parallel-session UI shape as Claude Code Agent View (news_log 2026-05-13) and convergent with our PM/Sci/Dev orchestration — but GitHub-anchored rather than terminal/IDE-anchored. Fits the "boards-as-control-plane" framing (cf. [OpenAI Symphony](#openai-symphony)). Tech preview only; too early for adoption decision but worth watching for GA + agent-callable surface (per [UI features ≠ agent capabilities](#our-position)). No action.

### Microsoft Conductor

- **status:** Hold
- **Source:** [Microsoft Open Source blog, 2026-05-14](https://opensource.microsoft.com/blog/2026/05/14/conductor-deterministic-orchestration-for-multi-agent-ai-workflows/); [microsoft/conductor](https://github.com/microsoft/conductor) (MIT, Python 3.12+)
- **Summary:** YAML-first CLI for multi-agent workflows. Routing between agents is **deterministic** — Jinja2 templates + expression evaluation handle conditions and branching; the orchestration layer itself spends zero LLM tokens. Static parallel groups (`fail_fast` / `continue_on_error` / `all_or_nothing`); dynamic for-each over variable-length arrays with batched concurrency. Backends: GitHub Copilot SDK or Anthropic Claude.
- **Our rationale:** **Explicit counter-position to LLM-as-orchestrator, partial echo of our shape.** Most multi-agent frameworks make the orchestrator itself an LLM that plans which agents to call — that pays token cost, latency, and unpredictability on every routing decision. Conductor argues: when workflow structure is *known*, declarative YAML + deterministic routing wins (diffable like CI/CD, zero token spend on routing). Our equivalent of "the routing layer" is the human + the GitHub Projects board + standup file — also zero LLM tokens spent on routing (the routing artifacts live in non-LLM substrates), also declarative-ish (Status / Priority / Size / Target Date fields are the YAML-analog). Difference: ours stays human-discretionary at every hop (PM judges triage, doesn't follow a fixed graph), Conductor's is fixed-at-definition-time. Useful as a thought-experiment: *which subset of PM's morning routine is workflow-known-enough to encode in a Conductor-style YAML?* Likely candidates: closure audit (mechanical), milestone health check (already script-driven). Unlikely candidates: triage scoping, news interpretation. No adoption action — too early, and Python-3.12+ pin doesn't fit our 3.11 baseline — but the deterministic-vs-LLM-orchestrator axis is now a first-class evaluation criterion when surveying future frameworks (per the user's *deterministic-before-semantic* preference, established 2026-05-21).

---

## Methodology / framing

Essays and reports that shape how the field talks about itself. Distinct from frameworks above — these change vocabulary and mental model, not tooling.

### Anthropic 2026 Agentic Coding Trends Report

- **status:** Reference
- **Source:** Anthropic 2026 (PDF; tracked via [Issue #235](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/235))
- **Summary:** Industry data on multi-agent coding workflows. Proposes canonical 4-specialist split: Architecture & Design / Implementation & Coding / Testing & Validation / Review & Docs.
- **Our rationale:** **Pattern borrow + roster reject.** Convergent pattern (one orchestrator + multiple specialists, each with own context) applies. Roster doesn't — that report is scoped to *coding*; we're scoped to bioinformatics research, where literature integration and biological interpretation have no analog in pure coding (so Scientist instead of Designer). Primary source for both [Role decomposition is domain-bespoke](#our-position) and [multi-role, not multi-agent](#our-position). Skim + 5 PM-practices cross-checks in [Issue #235](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/235).

### Mason — "Coherence Through Orchestration, Not Autonomy"

- **status:** Reference
- **Source:** [Mike Mason, Jan 2026](https://mikemason.ca/writing/ai-coding-agents-jan-2026/)
- **Summary:** Argues "multi-agent" coding works because of orchestration discipline, not agent autonomy.
- **Our rationale:** **Direct reinforcement.** Primary external source for [multi-role, not multi-agent](#our-position) — we don't claim autonomy, we claim disciplined orchestration. The user is the message bus; agents communicate via human-readable artifacts (standup, sub-issues). No action — already internalised.

### AddyOsmani — "The Code Agent Orchestra"

- **status:** Reference
- **Source:** [AddyOsmani blog, 2026-05-12](https://addyosmani.com/blog/code-agent-orchestra/)
- **Summary:** Essay on what makes multi-agent coding work — section roles, conductor metaphor.
- **Our rationale:** **Convergent vocabulary.** Maps cleanly onto our 3-layer orchestration: human (conductor) → PM (section leader) → Sci/Dev (instrumentalists). No new action items — this landscape doc itself ([Issue #337](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/337)) is the response to the essay.

### digitalapplied — "Pattern Language 2026"

- **status:** Trial
- **Next (-> Adopt / Rejected):** [Issue #353](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/353) decides whether to re-describe our shape in producer / consumer / coordinator / critic / judge terms.
- **Source:** [digitalapplied, 2026-05](https://www.digitalapplied.com/blog/multi-agent-orchestration-patterns-producer-consumer)
- **Summary:** Formal taxonomy reducing every reliable multi-agent system to 5 archetypes — producer (decomposes ambiguity into work items), consumer (executes them), coordinator (routes + bounded fan-out), critic (suggestions, no gate authority), judge (binary go/no-go). 12 composition rules (acyclicity, idempotent consumers, bounded fan-out, etc.); 8 predictable failure modes (cycle formation, critic-judge deadlock, coordinator overload, silent drift, others). Maps cleanly to LangGraph, CrewAI, OpenAI Agents SDK, Claude Agent SDK without changing the underlying shape.
- **Our rationale:** **Vocabulary candidate — open question.** Distinct from AddyOsmani above (essay vs. taxonomy) and from the Trends Report (coding-scoped vs. domain-agnostic). Possible re-description of our shape: `human = meta-coordinator → PM = coordinator + judge → Scientist = producer + critic + Developer = consumer + producer (infra) + critic`. The 8-failure-mode list doubles as a portfolio audit lens (have we hit critic-judge deadlock? silent drift?). Decision tracked via [Issue #353](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/353) (PM vocabulary adoption).

### "Toward Agentic Software Project Management: A Vision and Roadmap" (SPM 3.0)

- **status:** Reference
- **Source:** [arXiv 2601.16392](https://arxiv.org/html/2601.16392v1) (Jan 2026). Surfaced in the PM Signals scan 2026-07-07; Zotero methodology collection `DA3EWEJ9`.
- **Summary:** A vision/roadmap (explicitly "preliminary, without empirical evidence") for an **Agentic PM**, framed as a maturity model: SPM 1.0 (manual) -> 1.5 (tool-supported) -> 2.0 (GenAI copilot) -> **3.0 (a coordinating agent autonomously orchestrating task-specialized sub-agents under human oversight)**. Four capabilities (perceive / decide / act / learn-and-adapt via a central store) and **four working modes** scaling autonomy by task complexity x risk (Guided-AI-Autonomy -> Supervised-AI -> Human-AI-Collaborative -> AI-Assisted); it deliberately **excludes a fully-autonomous "Level 5"** and keeps accountability with the human PM (evolved into an "ethical / strategic leader / coach").
- **Our rationale:** **Direct reinforcement + advancement lens - our strongest external validation to date.** The paper theorizes as *future* empirical work almost exactly the system we already run: human-in-control, no Level-5 autonomy, a central learning store, task-specialized roles. Concept map: their working modes ~ our reversibility / blast-radius autonomy gate plus the autonomy-merge-gate cadence; their learn-and-adapt (envisioned as RL over a DB) ~ our file memory + episodic recall + the **mechanism-over-memory ladder** (an implemented, auditable learning loop the paper lacks). **The one deliberate divergence** is that theirs is the SPM-3.0 target and ours is not (yet): their *coordinating agent is an agent*; ours is the **human + the board** ([multi-role, not multi-agent](#our-position)). Where we lead the paper: we named the **binding constraint** (human review bandwidth, per [`agent_team_governance_research_2026-07.md`](agent_team_governance_research_2026-07.md)) that it never identifies, and we are the running instance it says is missing. **Two gaps it exposes in us:** (1) an explicit **ethics / data-integrity** governance dimension (they stress privacy + fabricated references; ours is engineering-process correctness, only partially covered by the Scientist `scientific-rigor` rule); (2) a **formalized autonomy / working-mode taxonomy** (we currently gate implicitly on reversibility). **Actionable spawn:** [epic #1072](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1072) scopes the advancement toward an agentic coordinating layer - closing the *decide-and-dispatch* loop the human currently owns - **review-axis first** (widen the binding constraint before automating dispatch), with the merge + novel-decision human gates preserved. Rungs: review-axis auto-request ([#1073](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1073)), autonomy-envelope spike ([#1074](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1074)), then Dev-led dispatch loops.

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

PM-owned.

**Funnel convention.** Every entry carries a `status:` ring (see [Entry status: the funnel](#entry-status-the-funnel)); active entries (`Assess` / `Trial`) carry a one-line DoR naming the event that promotes them.
The board-hygiene sweep greps for active entries whose promotion trigger has fired or whose tracked Issue went stale (most sweeps answer "nothing ready").

Update triggers:

- Morning-briefing item tagged `methodology-signal` surfaces → PM decides whether to add an entry here (ring it `Assess` with a DoR, or a terminal ring - `Reference` / `Hold` / `Rejected` - if no promotion path), or whether the item is chat-only and not landscape-doc-worthy.
- Board-hygiene sweep: `grep -nE "status: (Assess|Trial)" research/multi_agent_landscape.md` → for each active entry, has its DoR trigger fired (promote to `Trial` + open the eval Issue in the same edit), or has its tracked Issue gone stale?
- New framing memory in `pm/shared/` cross-references this doc → update [Our position](#our-position) to mention it.
- New Issue concretely borrows from or extends a framework listed here → backfill the framework's rationale block with the Issue link.
- Framework deprecated or pivoted → don't delete the entry; update the rationale with the deprecation note (the doc is a journal of *which* options we considered, not just current options).

Last sweep: 2026-07-07 (backfilled the SPM-3.0 vision paper [arXiv 2601.16392] to Methodology / framing as `Reference` - our strongest external validation of the multi-role model; captured the concept map + the two gaps it exposes (ethics/data-integrity dimension; formalized autonomy taxonomy) and seeded [epic #1072](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1072) (advance the coordinating layer toward agentic dispatch, review-axis first) with subs [#1073](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1073) / [#1074](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1074)).
Prior: 2026-07-02 (adopted the Technology Radar-style `status:` ring funnel [#553](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/553): retrofitted all entries with a ring + per-entry DoR on active ones, added the funnel convention + board-hygiene grep step; introduced a 6th `Reference` ring, confirmed by PM).
Prior: 2026-07-01 (companion governance-evidence report [`agent_team_governance_research_2026-07.md`](agent_team_governance_research_2026-07.md) cross-linked; adds MAST failure taxonomy + review-debt/review-column-WIP finding, neither previously tracked here).
Prior: 2026-05-21 (Microsoft Conductor backfilled to Frameworks as counter-position to LLM-as-orchestrator).
Prior: 2026-05-15 (GitHub Copilot desktop app backfilled to Frameworks as pattern-confirmation).
Prior: 2026-05-14 (Gartner cancellation forecast backfilled to `our_position` #2 as external validation).
Prior: 2026-05-12 (initial scaffold; 6 backfilled items: Symphony, Trends Report, Managed Agents, Mason, Fugu, Code Agent Orchestra).
