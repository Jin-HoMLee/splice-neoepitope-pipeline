# Governing an AI agent/role team - deep-research findings (2026-07)

PM-owned deep-research report.
Companion to [`multi_agent_landscape.md`](multi_agent_landscape.md) (the *frameworks/framing* survey); this doc is the *governance-evidence* survey - what the 2024-2026 literature says about **how to run** a team of AI personas, and which of our human-borrowed Kanban parameters are miscalibrated.

**Research question.**
How should PM/workflow governance be designed for a team of AI agents rather than humans, and which standard human-team Kanban/Agile assumptions break down when the "team" is AI personas?

**Method.**
Deep-research harness (fan-out web search across 6 angles -> 23 sources -> 109 extracted claims -> 3-vote adversarial verification -> synthesis).
Rate-limiting knocked out roughly half the verification votes on the final run, so the four decision-relevant claims that landed in the errored bucket were **re-verified by hand** against their primary sources (inline `WebFetch` fact-checks, 3/4 confirmed verbatim; the 4th was a non-actionable sub-statistic whose parent claim is independently confirmed).
Every claim below is either 3-0 confirmed by the harness or hand-verified; confidence and source-tier are noted per finding.

**Raw results / provenance.**
The complete machine output - all synthesized findings, every refuted and unverified claim verbatim, source-by-source breakdown, per-run logs across all three runs - is preserved in [`agent_team_governance_2026-07_deepresearch_raw.json`](agent_team_governance_2026-07_deepresearch_raw.json).
This doc is the *curated* read of that raw trail.

---

## The one-line reframe

For an AI-persona team the load-bearing parts of Kanban/Agile shift **from managing human fatigue and one-way flow to managing coordination cost, state persistence, and verification.**
Our *structure* already matches best practice; the *dials we borrowed from human teams* are what's miscalibrated.

## Reconciliation with our "multi-role, not multi-agent" position

The literature I pulled is the *multi-agent-systems* literature, but our [stated position](multi_agent_landscape.md#our-position) is that we are a **human-orchestrated multi-role workflow**, not an autonomous multi-agent system - the human is the message bus and personas never trigger each other.
That distinction is not a dodge here; it's load-bearing, because **our design already immunizes us against a whole class of the failure modes the literature describes:**

- **Inter-agent misalignment / infinite handoff loops / silent hallucination propagation** are largely *autonomy* failures - they need agent-to-agent triggering to occur.
  With a human at every hop, these degrade from "silent runaway" to "a human reads the handoff artifact and catches it."
  Our shared board is a *stronger* coordination substrate than the isolated-subagent pattern Anthropic describes (their subagents "can't coordinate mid-task"; ours coordinate through a durable, human-read board).
- What our design does **not** immunize against, and what therefore deserves the governance attention: **specification quality at handoff, verification of claimed completion, and state/memory staleness across context resets.**
  Those are exactly the findings that survived verification below.

So: read every "multi-agent failure mode" finding as *"the risk if we drifted toward autonomy"* - and note that most of our governance value is in staying on the disciplined-orchestration side of that line.

---

## Confirmed findings

### Structure - our shape is validated

1. **Orchestrator-worker is the canonical coordination architecture** and maps onto our human -> PM -> Sci/Dev dispatch. *(3-0; primary: Anthropic engineering.)*
   Caveat from verification: Anthropic's subagents can't coordinate mid-task, so the mapping is an analogy, not an identity - and our shared board is the *stronger* substrate.
2. **Structured, explicit task specs at handoff are essential** - "an objective, an output format, guidance on the tools and sources to use, and clear task boundaries"; without them "agents duplicate work, leave gaps, or fail to find necessary information." *(3-0; primary: Anthropic.)*
   This maps 1:1 onto our **Definition of Ready** and the `Backlog -> Ready` commitment gate. DoR *is* our handoff contract - harden it.
3. **External memory is load-bearing infrastructure, not overhead.** Long-running agents "summarize completed work phases and store essential information in external memory before proceeding"; past ~200K tokens context truncates. *(3-0; primary: Anthropic.)*
   This directly justifies the **Memory Manager role + on-disk memory stores**. Corroborating (hand-verified): "context rot" degrades accuracy *before* the nominal limit - a Chroma study of 18 frontier models found degradation at every increment, a 30%+ lost-in-the-middle drop, and significant degradation at 50K tokens despite a 200K window. *(hand-verified verbatim; secondary.)*

### Where the ROI is - governance, not model

4. **Most multi-agent failure is architectural/governance, not raw model capability.** The MAST empirical taxonomy (150+ traces, 7 frameworks, kappa=0.88, 1600+ annotated traces) clusters 14 failure modes into system-design (~44%), inter-agent misalignment (~32%), and task verification (~24%), and argues failures need "structural MAS redesigns," not a stronger base model. *(3-0; primary: arXiv 2503.13657.)*
   This is the strongest single justification that **investing in board governance itself is the high-ROI lever** - the payoff frontier is process design, not waiting for better models.
5. **Verification gaps are a quantified, high-prevalence failure class** ("Incorrect Verification" ~9.1%; the Task Verification category ~21-24% of failures across sources), flagged as "a significant challenge regardless of the LLM used." *(3-0; primary MAST + secondary Augment Code, both verbatim.)*
   Validates our **closure-ritual gate, `audit_and_merge` checks, and the review columns** as the mechanisms that catch hallucinated/premature "done."
6. **Human oversight at key decision points is load-bearing, not optional.** *(3-0; arXiv 2601.13671 - already in our landscape doc as the Gartner-cancellation source.)*
   For a 1-human/N-persona team the actionable read is: **concentrate the human's scarce attention at irreversible gates (merge, commitment), not spread across every status move.**

### The miscalibrated dials - what we borrowed from humans

7. **The binding constraint is human *review* bandwidth, not agent throughput ("review debt").**
   "The more you parallelize your code generation, the more review debt you create" (corroborated: PR volume +29% YoY post-AI; AI code ~12x costlier to review). *(3-0; secondary: InfoWorld + O'Reilly "comprehension debt".)*
   Practitioner guidance, **hand-verified verbatim**: apply WIP limits to the *human-review* columns separately from agent columns - "limit 'Awaiting Human Review' to 5-10 cards per reviewer"; "Human attention is also a constrained resource. Cap it." *(hand-verified; blog-tier, but the parameter is now confirmed as stated.)*
   -> Our per-role in-progress cap (~3) is defensible (it curbs duplicate-PR swarms - practitioner-supported), but it is aimed at the *wrong* constraint. The higher-value guard is a **WIP limit on the review/In-review columns keyed to the single human's review capacity.**
8. **Output/velocity metrics measure the wrong thing.** "Celebrating commit velocity ... you are not measuring productivity [but] how quickly your team can manufacture liability"; code generation "was never the bottleneck ... validation is." *(3-0; secondary: InfoWorld. Note: the *stronger* claim that agent metrics are "trivially gameable" was REFUTED 1-2 - the surviving, weaker claim is "measures the wrong thing.")*
   -> Our current "no flow data yet" state is *fine*; wall-clock throughput was never the right instrument. If we instrument anything, measure **verification lead time and rework rate**, not generation volume.

### New failure modes to design against (autonomy-conditional for us)

9. **Multi-agent-specific failure modes**: hallucination *amplification* (a hallucinated output in shared memory treated downstream as verified fact, "no error signal"), context loss at handoffs, cascading failures. *(3-0; primary: arXiv 2601.13671.)*
   For us the propagation vector is the shared board/memory - so **verify before one persona treats another's board/memory output as ground truth.** (Mitigated but not eliminated by the human message bus.)
   Hand-verified extension worth noting: an **independent judge agent** with isolated context/scoring is a recommended verification mechanism (PwC reported 7x accuracy, 10%->70%, via structured validation loops). *(hand-verified verbatim; secondary.)* This is the *semantic* sibling of our already-shipped *mechanical* post-merge critic ([Issue #325](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325)).

---

## What this changes - operationalization candidates

Proposed as Backlog options under `arc:board-governance` (not silent retuning; #902-style). Detail in the "Issues to file" section of the originating session.

1. **Review-column WIP limit** (start 5-10 cards/reviewer) - the review-debt fix; re-anchors WIP to the real (human-review) constraint. *Confirmed evidence.*
2. **Retire the "advisory until we have flow data" TODOs**; reframe the metric target as verification lead time + rework, not throughput. *Confirmed evidence.*
3. **(Stretch) semantic critic/judge pass at the merge gate** for hallucinated-completion detection - extends the shipped mechanical critic (#325), connects to the critic/judge vocabulary question (#353). *Confirmed as a pattern; needs a fit/cost check at our scale.*

## Caveats

- Field is fast-moving; all sources 2025-2026. MAST percentages and Anthropic's token multipliers are point-in-time and system-specific (self-reported approximations for their own product).
- Source-tier split: structure/failure-taxonomy/memory/human-oversight rest on **primary** sources (unanimous 3-0); the two *dial* findings (review debt, output-metric failure) rest on **secondary/blog** sources - directionally strong and corroborated, but conceptual, not controlled measurement. Treat the review-column parameter (5-10) as a defensible *starting* value to tune, not a law.
- One sub-statistic (individual MAST mode prevalences) remains unverified (arXiv PDF wouldn't parse in the fetch tool); its parent claim (aggregate MAST categories) is confirmed. Not decision-relevant.
- Framings like "orchestrator-worker maps to our PM board" and "DoR is the handoff contract" are the researcher's interpretive mappings, not source claims.

## Sources (verified subset)

- Anthropic, "How we built our multi-agent research system" (primary) - https://www.anthropic.com/engineering/multi-agent-research-system
- Cemri et al., "Why Do Multi-Agent LLM Systems Fail?" MAST taxonomy, arXiv 2503.13657 (primary)
- arXiv 2601.13671 (primary; also our landscape doc's Gartner-forecast source)
- InfoWorld, "AI agents and bad productivity metrics" (secondary) - https://www.infoworld.com/article/4135492/ai-agents-and-bad-productivity-metrics.html
- MindStudio, "Iterative Kanban pattern for AI agents" (blog; review-column WIP parameter) - https://www.mindstudio.ai/blog/iterative-kanban-pattern-ai-agents-feedback-loop
- Zylos, "Context window management / session lifecycle" (secondary; context rot) - https://zylos.ai/research/2026-03-31-context-window-management-session-lifecycle-long-running-agents/
- Augment Code, "Why multi-agent LLM systems fail" (secondary; judge agent, MAST categories) - https://www.augmentcode.com/guides/why-multi-agent-llm-systems-fail-and-how-to-fix-them
