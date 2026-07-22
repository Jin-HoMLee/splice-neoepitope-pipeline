# Developer-facing agent frameworks + multi-agent reliability evals - deep-research findings (2026-07)

PM-owned deep-research report.
Second companion to [`agent_team_governance_research_2026-07.md`](agent_team_governance_research_2026-07.md) (the *governance-evidence* survey) and [`multi_agent_landscape.md`](multi_agent_landscape.md) (the *frameworks/framing* survey).
This doc closes the **residual** of [Issue #265](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/265): the two sub-questions the governance report and landscape doc did not cover - (1) the developer-facing orchestration frameworks (LangGraph/CrewAI/AutoGen/OpenAI Agents SDK), and (2) multi-agent reliability benchmarks and evals.

**Research question (scoped to the residual).**
For a small human-orchestrated persona team (1 human message bus + 3-4 personas), which developer-facing multi-agent orchestration frameworks add anything over our shared-board + human-bus + on-disk-memory design, and what reliability benchmarks/evals actually transfer to measuring our own setup?
The prior report already settled failure modes (MAST), human-in-the-loop placement, the review-debt/WIP dials, and external-memory importance; those were explicitly out of scope here.

**Method.**
Deep-research harness (fan-out web search across 6 angles -> 25 sources -> 119 extracted claims -> 3-vote adversarial verification -> synthesis).
24 of 25 verified claims confirmed, 1 refuted (0-3), 0 unverified.
Every finding below is either 3-0 confirmed or explicitly flagged where a vote split 2-1.

**Raw results / provenance.**
The complete machine output (all synthesized findings, the refuted claim, per-source breakdown, per-agent logs) is preserved in [`agent_frameworks_reliability_2026-07_deepresearch_raw.json`](agent_frameworks_reliability_2026-07_deepresearch_raw.json).
This doc is the curated read of that trail.

---

## The one-line reframe

The frameworks do not upgrade our design - they invert it, and the reliability benchmarks do not transfer - but their metrics do.
Both sub-questions resolve the same way: our architecture is validated by contrast, and the actionable value is a handful of cheap methods, not a new tool.

---

## Sub-question 1: developer-facing orchestration frameworks

### The frameworks converge on two primitives, and both replace the human coordinator

The developer-facing frameworks (OpenAI Swarm / Agents SDK, Microsoft AutoGen v0.4+, and by extension CrewAI and LangGraph) converge on two coordination primitives: graph/pub-sub routing and handoff-via-tool-call delegation.

- **OpenAI Swarm is experimental/educational and superseded** by the production-oriented OpenAI Agents SDK (March 2025). *(3-0; primary: OpenAI Swarm README + Agents SDK repo.)*
  Swarm itself is not an adoption candidate; the Agents SDK is its production evolution.
- **Swarm's coordination model is handoff-based**: an agent transfers control by having a called function return another Agent object, optionally bundling return values plus context updates in one Result object - dynamic decentralized routing with no central router. *(3-0; primary: Swarm README.)*
- **AutoGen v0.4+ reimplements the same handoff mechanism**: agents delegate via a special delegate/handoff tool call rather than through a central router. *(3-0; primary: AutoGen core docs, self-attributed to Swarm.)*
- **The routing is decentralized and machine-mediated**: each agent decides whether to delegate; routing occurs through topic-based pub-sub subscriptions managed by the runtime; there is no central router *agent*. *(2-1; the precise reading is "no central router agent," not "no central message bus" - the pub-sub bus is itself the runtime.)*

### The load-bearing adoption finding

**AutoGen treats a human as a symmetric pub-sub endpoint (a `HumanAgent`), not a coordinator - the direct inverse of our design.** *(3-0; primary: AutoGen docs.)*
The `HumanAgent` "subscribes to `agent_topic_type` to receive messages and publishes to `user_topic_type`... a message endpoint/participant, not a central coordinator... operating through the same event-driven pub-sub infrastructure as AI agents, rather than as a special central authority."

For a 1-human/N-persona team whose entire design is human-as-bus plus human-read board, these frameworks' primitives are architecturally orthogonal: they *automate away* the human coordinator our design deliberately keeps in the seat.
Adopting one wholesale would require inverting our model, not extending it.

### What transfers

The transferable increment is a **named handoff/delegate convention on the board** (persona A hands to persona B with an explicit context bundle), not a framework runtime.
That is the one primitive worth borrowing, and it maps onto the coordination beat the governance report flagged as our weakest (handoff/conflict-resolution).

### Honest coverage gap

CrewAI and LangGraph were named in scope but produced **no surviving verified claims** in this run.
The framework evidence that survived is entirely Swarm/Agents-SDK and AutoGen; conclusions about CrewAI's role/crew delegation and LangGraph's explicit durable-state graph are inferred by analogy, not directly evidenced.
If a coordination primitive from either (LangGraph's durable shared-state graph is the most plausible candidate) is ever a real adoption question, it needs its own verified pass.

---

## Sub-question 2: multi-agent reliability benchmarks + evals

### Headline: single-run benchmark scores systematically overestimate real reliability

**Pass@1 is structurally blind to reliability.** *(3-0; primary: arXiv 2603.29231, arXiv 2601.06112.)*
GPT-4o scores 61% pass@1 on retail agent tasks but only ~25% pass^8 - the probability of at least one failure across 8 runs approaches 75%.
Single-run benchmarks (ToolBench, AgentBench, API-Bank) systematically overestimate production reliability.

### The benchmarks are gameable as-is

**Agent benchmarks are frequently and materially unreliable.** *(3-0; primary: arXiv 2507.02825, 2510.11977, 2605.12673, tau2-bench repo.)*
Design flaws under/overestimate performance by up to 100% in relative terms.
Named mainstream benchmarks have concrete validity bugs (SWE-bench Verified uses insufficient test cases; TAU-bench counts empty responses as success).
Agents shortcut tasks (searching HuggingFace for the benchmark instead of solving it).
tau2-bench itself shipped 75+ task-quality fixes.
Automated red-teaming (BenchJack) reward-hacked **9 of 10 popular benchmarks to near-perfect scores without solving a single task** (e.g. a 9-line `conftest.py` exploit -> 100% on SWE-bench Verified).

**Bottom line: do not treat headline agent-benchmark numbers as predictive of our setup's reliability.**
The benchmark *domains* (retail/airline/telecom/SWE/web) are large-software-org-specific and do not match research-bioinformatics work; only the metrics and methods transfer.

### LLM-simulated users/judges are miscalibrated

On tau-Bench retail with the agent fixed (GPT-4o), swapping only the user-simulator LLM shifts measured agent success by ~9 percentage points. *(3-0; primary: arXiv 2601.17087.)*
Simulated users are miscalibrated against real humans (ECE 15.1), underestimating success on the hardest tasks and overestimating on moderate ones.
Any harness we build that uses an LLM as user or as judge inherits the judge's biases.

### What transfers to us - the methods, not the benchmarks

1. **pass^k / repeat-run consistency** is the cheapest reliability instrument.
   Re-run the same persona beat N times and count *all-succeed*, not once.
   This directly extends the prior report's "rework rate" recommendation with a board-native measurement.
2. **Full-execution-trace failure attribution** (TraceElephant, ACL 2026).
   Providing full execution traces (inputs + context + handoff chain, not just outputs) improves failure-attribution accuracy by **up to 76%** over output-only observation. *(3-0; primary: arXiv 2604.22708.)*
   When a persona beat fails, capturing the full trace lets us attribute *which* persona and *which* step caused the rework - the prior report's rework-rate becomes actionable.
3. **Partial-credit / graceful-degradation scoring** (GDS) for partially-completed long tasks. *(3-0; primary: arXiv 2603.29231.)*
   A concept a small team can apply qualitatively without heavy infrastructure.
4. **A three-axis stress model** (consistency / robustness / fault-tolerance) from ReliabilityBench. *(3-0 vote, single-author preprint, medium source confidence; arXiv 2601.06112.)*
   Rerun a beat with reworded instructions (robustness) and under a simulated tool/network failure (fault tolerance).

### The negative result worth heeding

**A naive episodic-memory scaffold universally hurt long-horizon reliability** across all 10 models tested (4 neutral, 6 hurt, none positive). *(3-0; single preprint: arXiv 2603.29231.)*
Our curated, structured on-disk role memory is *not* a naive episodic scratchpad, so this is a caution, not a condemnation.
But it directly warns against assuming "more agent memory = more reliable," and it makes a memory-on/memory-off A/B on a few representative beats a defensible thing to run before we treat our memory layer as self-evidently net-positive for reliability.

---

## What this changes - operationalization candidates

Proposed as Backlog options under `arc:board-governance` (the follow-up Issues below).

1. **pass^k self-test** on a representative persona beat - re-run the same handoff/dispatch task 3-5x and count all-succeed; the cheapest reliability instrument, extends rework-rate.
2. **Full-execution-trace capture per failed/reworked beat** (TraceElephant-style) - capture inputs/context/handoff chain, not just the final output, so a rework is attributable to a persona and step.
3. **(Stretch) memory-on/memory-off A/B** - given the naive-episodic-memory negative result, test whether our curated memory layer measurably helps a few representative beats before assuming it does.

## Caveats

- Several reliability findings rest on single non-peer-reviewed 2025-2026 arXiv preprints (2603.29231, 2601.06112 single-author, 2605.12673); treat their specific numbers as directional, not settled.
  The robust, multiply-corroborated core is the qualitative finding that single-run scores overestimate real reliability (traces to Sierra's tau-bench pass^k, echoed across HAL, 2507.02825, and BenchJack).
- The framework landscape moves fast (Swarm -> Agents SDK, AutoGen v0.4); adoption/star signals shift monthly.
  The architectural *conclusion* (frameworks automate the human coordinator away) is stable; specific maturity claims should be re-checked before any adoption decision.
- Two votes were 2-1 (AutoGen "no central message bus" -> precise reading "no central router agent"; BenchJack -> could not confirm which specific benchmarks were among the 9 hacked). Both flagged inline.
- CrewAI and LangGraph produced no surviving verified claims despite being in scope (see the coverage gap above).
- The lightweight CI-native eval-tooling claim (Promptfoo) was refuted 0-3, so no verified evidence on that tooling survived; the pass^k / trace-capture increments above stand on their own primary sources.
- Nothing here directly benchmarks a research-bioinformatics workflow; all named benchmarks are large-software-org domains, so their scores do not transfer - only their metrics and methods do.

## Sources (verified subset)

- OpenAI, Swarm README (primary) - https://github.com/openai/swarm/blob/main/README.md
- OpenAI Agents SDK (primary) - https://github.com/openai/openai-agents-python
- Microsoft AutoGen, Handoffs design pattern (primary) - https://microsoft.github.io/autogen/dev//user-guide/core-user-guide/design-patterns/handoffs.html
- arXiv 2603.29231 - reliability metrics (RDC/VAF/GDS/MOP) + episodic-memory negative result (primary; preprint)
- arXiv 2601.06112 - ReliabilityBench, R(k,epsilon,lambda) surface (primary; single-author preprint)
- arXiv 2406.12045 - tau-bench (Yao et al., ICLR 2025), the peer-reviewed origin of the 61% pass@1 -> ~25% pass^8 headline number (primary; the two preprints above only corroborate it)
- arXiv 2506.07982 / sierra-research/tau2-bench - dual-control agent benchmark (primary)
- arXiv 2510.11977 - Holistic Agent Leaderboard (HAL) (primary; ICLR 2026)
- arXiv 2604.22708 - TraceElephant, failure attribution with full traces (primary; ACL 2026)
- arXiv 2601.17087 - "Lost in Simulation," LLM-user miscalibration (primary)
- arXiv 2507.02825 - Agentic Benchmark Checklist (ABC), benchmark validity audit (primary)
- arXiv 2605.12673 - BenchJack, automated benchmark reward-hacking (primary; preprint)
