# News Log

Shared log of news surfaced in morning briefings. Prevents repeating the same item across sessions and roles (PM, Scientist, Developer all log here).

**Format:**
- `## YYYY-MM-DD` — date-level section (one per day; new dates at the top)
- `### HH:MM UTC — Editor: <Role>` — time + editor attribution (one per session per day)
- `- **Item** (date) — keywords. → action. *Role. Signal.*` — one bullet per item, nested under the time section

Historical entries (before 2026-05-02) pre-date the time/editor sub-heading convention and are kept as-is.

---

## 2026-05-19

### 09:25 UTC — Editor: Scientist

- **Fisher et al.** ([JCI Insight 2025](https://insight.jci.org/articles/view/190287)) — Corin/CoREST inhibitor → drug-induced splice-NA + ICI synergy in melanoma; preclinical. → DISCUSSION ref; glossary + Zotero (`BIAEVMQ9`) same PR. *Scientist. methodology-signal.*
- **Shapiro & Bassani-Sternberg** ([Trends in Cancer 2026](https://doi.org/10.1016/j.trecan.2026.02.003)) — Immunopeptidomics-for-NA-discovery review. → Backfill; already in Zotero `Q6GIXN8N`. *Scientist. methodology-signal.*

### 08:20 UTC — Editor: PM

- **Claude Code 2.1.143** ([2026-05-15 changelog](https://code.claude.com/docs/en/changelog)) — `PostToolUse` hooks now break after 8 consecutive blocks (safety against runaway deny loops); new `worktree.bgIsolation: "none"` setting lets background sessions edit working copy directly. → No action; loop-break is free safety for our [PR #397](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/397) capacity-recheck hook + `check_at_claude.py` PreToolUse guard. *PM. tooling-relevant.*
- **Copilot cloud agent — REST API for repo config audit** ([GitHub changelog 2026-05-18](https://github.blog/changelog/)) — agent configuration on a repo (allowed actions, model defaults, secrets exposure) now queryable via REST without UI scraping. → Concrete hook for [Issue #295](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/295) (Copilot eval for PM micro-chores) — adds an auditable surface before any delegation trial. *PM. tooling-relevant.*

## 2026-05-18

### 08:15 UTC — Editor: Scientist

- **Garfinkle et al. — GNOS-PV01 DNA neoantigen vaccine** ([Nat Cancer 2026](https://www.nature.com/articles/s43018-026-01163-w)) — Phase 1 MGMT-unmethylated GBM, n=9, up to 40 NAs/patient; no DLTs, 2/3 PFS@6mo. → DISCUSSION ref; glossary bundle same PR. *Scientist. clinical-validation.*

## 2026-05-17

### 18:49 UTC — Editor: Developer

- **Snakemake 9.21.0** ([2026-05-14](https://github.com/snakemake/snakemake/releases/tag/v9.21.0)) — sqlite PRAGMA tuning + auto-detect network fstype (perf), logger plugin handler fix, filename arg-prepend helper. → No action; 8→9 migration target keeps moving while [Issue #200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200) is pending. *Developer. pipeline-relevant.*
- **Astral + prefix-dev toolchain cadence** — [ruff 0.15.x](https://github.com/astral-sh/ruff/releases) (now formats to 2026 style guide), [pixi 0.66.x](https://pixi.sh/dev/CHANGELOG/) (registered workspaces), [uv 0.11.14](https://github.com/astral-sh/uv/releases). → No action. *Developer. industry-standard.*

### 13:13 UTC — Editor: PM

- **GitHub Projects — Timestamp fields GA** ([changelog 2026-05-15](https://github.blog/changelog/2026-05-15-timestamp-fields-in-github-projects/)) — three built-in fields (Created, Updated, Closed) addable to any project view; agent-callable via GraphQL for "recently completed" views. → No action; closure-audit query (`--search "closed:>YYYY-MM-DD"`) already works. *PM. tooling-signal.*
- **Claude Code — Agent `subagent_type` matching + `/goal` evaluator fix** ([changelog](https://code.claude.com/docs/en/changelog)) — case-/separator-insensitive resolution ("Code Reviewer" → `code-reviewer`); `/goal` no longer fires while background shells or delegated subagents are still running. → No action; quality-of-life only. *PM. tooling-relevant.*

## 2026-05-15

### 16:42 UTC — Editor: PM

- **GitHub Copilot desktop app — technical preview** ([changelog 2026-05-14](https://github.blog/changelog/2026-05-14-github-copilot-app-is-now-available-in-technical-preview/)) — GitHub-native desktop for agentic dev from issue/PR/prompt; convergent with Claude Code Agent View. → Landscape doc Frameworks backfill (same PR). *PM. methodology-signal.*
- **`claude agents` flag expansion** ([Claude Code changelog](https://code.claude.com/docs/en/changelog)) — `--add-dir`/`--settings`/`--mcp-config`/`--plugin-dir`/`--permission-mode`/`--model`/`--effort` on dispatched background sessions. → No action; revisit if PM ever dispatches sub-agents. *PM. tooling-relevant.*

### 14:08 UTC — Editor: Scientist

- **Manakongtreecheep et al. — Pepyrus** ([Nat Biotech 2026](https://www.nature.com/articles/s41587-026-03003-9)) — *E. coli* user-defined peptide library + DIA-MS; 0.1 fmol rare-neoantigen detection. → Zotero add; DISCUSSION MS-validation hook; glossary DDA + DIA-MS same PR. *Scientist. portfolio-differentiator.*
- **Prélot et al.** ([bioRxiv 2025-09-16](https://www.biorxiv.org/content/10.1101/2025.09.10.674685v1)) — splice-neoepitope 69-86% cross-pipeline divergence benchmark; already in Zotero `ZXAUQAJL` + [Issue #304](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/304). → Audit-trail backfill. *Scientist. audit-trail.*

## 2026-05-14

### 10:17 UTC — Editor: PM

- **Claude Code May 2026 hook updates** ([changelog](https://code.claude.com/docs/en/changelog)) — `updatedToolOutput` now all-tools (was MCP-only); `terminalSequence` for desktop notifications. → comment on [Issue #247](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/247) (design refinement options). *PM. tooling-relevant.*
- **Multi-agent framework leaderboard 2026** ([Gurusup](https://gurusup.com/blog/best-multi-agent-frameworks-2026)) — AutoGen 32k★ / CrewAI 22k★ / LangGraph 8.2k★ on GitHub. → No action; pattern-only. *PM. methodology-signal.*
- **Gartner: >40% agentic AI projects cancelled by 2027** ([arXiv 2601.13671](https://arxiv.org/html/2601.13671v1)) — primary failure = orchestration ↔ domain-logic junction. → Backfill `our_position` #2 rationale in [multi_agent_landscape.md](research/multi_agent_landscape.md). *PM. methodology-signal.*

## 2026-05-13

### 09:00 UTC — Editor: Developer

- **PyTorch 2.12 cu126 keeps Pascal** ([CUDA 13.2 thread](https://dev-discuss.pytorch.org/t/introducing-cuda-13-2-and-deprecating-cuda-12-8-release-2-12/3337), [RFC #178665](https://github.com/pytorch/pytorch/issues/178665)) — corrects yesterday's "no revival": cu126 wheels list SM 6.0 for torch 2.8–2.12; `torch<2.5` permanence loose. → [Issue #352](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/352) (P100 spike). *Developer. pipeline-relevant.*

### 08:25 UTC — Editor: PM

- **Pattern Language 2026 — 5-archetype multi-agent taxonomy** ([digitalapplied 2026-05](https://www.digitalapplied.com/blog/multi-agent-orchestration-patterns-producer-consumer)) — producer/consumer/coordinator/critic/judge; 12 composition rules, 8 failure modes; maps to LangGraph/CrewAI/Claude Agent SDK. → Landscape doc Methodology entry via [PR #344](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/344) (parent landscape PR); vocabulary decision via [Issue #353](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/353). *PM. methodology-signal.*
- **Claude Code Agent View + `/goal` command** ([Anthropic 2026-05-11](https://claude.com/blog/agent-view-in-claude-code), v2.1.139) — Agent View = CLI dashboard + VS Code Agent Sessions sidebar listing all parallel sessions (running/blocked/done); `/goal <verifiable condition>` = autonomous multi-turn primitive iterating until met or constraints hit. → No action; revisit `/goal` for [Issue #200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200) (Snakemake 9.x migration — verifiable via dry-run). *PM. tooling-relevant.*

### 07:35 UTC — Editor: Scientist

- **Sahin et al.** ([Nature 2026](https://www.nature.com/articles/s41586-025-10004-2)) — Phase 1 TNBC personalized mRNA vaccine; 11/14 relapse-free at 6y; BioNTech (cf. Sahin 2017). → [Issue #351](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/351) (Sci, DISCUSSIONS citation). *Scientist. clinical-validation.*

## 2026-05-12

### 09:31 UTC — Editor: PM

- **GitHub Projects "Hierarchy view" GA** (2026-03-19, [changelog](https://github.blog/changelog/2026-03-19-hierarchy-view-in-github-projects-is-now-generally-available/)) — nested sub-issues now render on the board (8 levels), inline create + drag-to-reparent. → User-facing only; agent-side `gh api .../sub_issues` unchanged (new rule `feedback_ui_vs_agent.md`). *PM. tooling-signal.*
- **AddyOsmani — "The Code Agent Orchestra"** ([blog post](https://addyosmani.com/blog/code-agent-orchestra/)) — essay on what makes multi-agent coding work; convergent with PM/Sci/Dev orchestrator pattern. → [Issue #337](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/337) (PM, landscape doc + backfill). *PM. methodology-signal.*

### 09:06 UTC — Editor: Scientist

- **Onkar et al. — Bhardwaj-lab cancer vaccine field synthesis** ([Cell Rep Med 2026](https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(25)00648-2)) — neoantigen + ICI convergence; off-the-shelf shared-NA vaccines named as emerging strategy. → [Issue #334](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/334) (Sci, DISCUSSION clinical-translation framing). *Scientist. portfolio-differentiator.*

### 08:57 UTC — Editor: Developer

- **Snakemake 9.20.0 + 8→9 migration is small** ([changelog](https://snakemake.readthedocs.io/en/stable/project_info/history.html), [migration guide](https://snakemake.readthedocs.io/en/stable/getting_started/migration.html)) — only 1 breaking change (custom logger API via `--logger` / `OutputSettings.log_handler_settings`); much lighter than 7→8. → Re-scope [Issue #200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200) (Snakemake 9.x migration) lower. *Developer. pipeline-relevant.*
- **PyTorch 2.12 ships 2026-05-13** ([release key dates](https://dev-discuss.pytorch.org/t/pytorch-release-2-12-key-dates/3329)) — CUDA 13.2 experimental + Blackwell expansion; no Pascal revival. → No action; `torch<2.5` pin stays permanent on P100. *Developer. pipeline-relevant.*
- **libdeflate 1.26 still missing from bioconda samtools** (deps still pinned `<1.26.0a0`) — `hisat2.yaml` system-samtools workaround remains permanent. → No action; revisit [Issue #237](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/237) when bioconda lifts the cap. *Developer. pipeline-relevant.*

## 2026-05-11

### 10:04 UTC — Editor: Developer

- **Snakemake Hackathon 2026** ([BioHackrXiv 2026-04-27](https://index.biohackrxiv.org/2026/04/27/h6zqj.html)) — plugin / HPC / cloud-backend focus. → [Issue #66](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/66) (Google Batch) + [Issue #200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200) (Snakemake 9.x). *Developer. pipeline-relevant.*

### 09:11 UTC — Editor: PM

- **GitHub Issues from Slack — `@GitHub` natural-language mentions** ([May 2026 Changelog](https://github.blog/changelog/)) — NL mention in any channel creates structured issue in repo. → No action; we don't use Slack. *PM. tooling-signal.*
- **Sakana Fugu — multi-agent orchestration as a foundation model** ([Sakana 2026](https://sakana.ai/fugu-beta/)) — beta dynamically assembles agent teams from a pool, no fixed roles. Opposite of our PM/Sci/Dev split. → 3 borrow-Issues filed: [Issue #324](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/324) (model routing), [Issue #325](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/325) (post-merge critic), [Issue #326](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/326) (memory consolidation). *PM. methodology-signal.*

## 2026-05-10

### 15:17 UTC — Editor: PM

- **GitHub repo rulesets — per-user bypass actors** ([May 2026 Changelog](https://github.blog/changelog/)) — individuals as bypass actors via UI/REST/GraphQL. → No action; surgical alt to "branches up to date" rule removed via [PR #313](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/313) (rule removal). *PM. tooling-signal.*

### 08:43 UTC — Editor: Developer

- **OpenAI acquired Astral (uv/ruff/ty)** ([2026-03-19, dev.to](https://dev.to/max_quimby/openai-just-acquired-astral-what-it-means-for-uv-ruff-and-every-python-developer-41ah)) — Astral team joins Codex; uv/ruff/ty stay open-source under active development. → No action; revisit if we adopt uv/ruff for CI linting. *Developer. industry-standard.*
- **uv 0.11.12** ([2026-05-08 release](https://github.com/astral-sh/uv/releases)) — minor patch bump from 0.11.8 (logged 2026-05-01). → No action. *Developer. industry-standard.*
- **PyTorch 2.11.0** ([March 2026 release](https://github.com/pytorch/pytorch/releases)) — Pascal/Maxwell removal still stands per 2.8 cut; no change to our `torch<2.5` pin permanence on P100. → No action; pin permanent. *Developer. pipeline-relevant.*

### 08:30 UTC — Editor: Scientist (back-filled 2026-05-11)

- **Lu et al. — Benchmarking TCR-pMHC structure prediction** ([bioRxiv 2025-11-30](https://www.biorxiv.org/content/10.64898/2025.11.30.691400v1.full)) — n=70 unseen complexes, 10 methods; AF3 leads on accuracy + docking; CDR3-pLDDT rerank gives +4.3% Top-1 success on medium-quality predictions. → [Issue #316](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/316) (Sci, AF3 + CDR3-pLDDT eval). *Scientist. portfolio-differentiator.*

## 2026-05-09

### 09:21 UTC — Editor: PM

- **Copilot CLI cross-family Rubber Duck** ([release](https://releasebot.io/updates/anthropic/claude-code)) — Claude critic dispatched when GPT is orchestrator; formalizes our manual dual-LLM review pattern. → No action. *PM. methodology-signal.*
- **Anthropic doubles Pro/Max/Team limits, peak-hour reduction removed** ([release](https://releasebot.io/updates/anthropic)) — caps lift; relevant to Max→Pro downgrade after evening-compact usage drop. → Track usage 1–2 weeks. *PM. ops-signal.*

### 09:12 UTC — Editor: Scientist

- **Kim et al.** ([Cell 2025](https://www.cell.com/cell/fulltext/S0092-8674(25)00399-X)) — SF3B1/SRSF2 mis-splicing → public neoantigens + cognate TCRs vs SRSF2-mutant leukemia. → [Issue #311](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/311) (Sci, Kwok DISCUSSION one-liner). *Scientist. portfolio-differentiator.*
- **Zhang et al. — Safety-first TCR-T** ([Front Immunol 2026](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2026.1754735/full)) — AI specificity + base-editing anti-mispairing + public-NA QC <0.5 TPM normal-tissue. → DISCUSSION ref for [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299) threshold tradeoff. *Scientist. methodology-signal.*

## 2026-05-08

### 13:18 UTC — Editor: PM

- **Anthropic Managed Agents — Dreaming + Multiagent Orchestration** ([2026-05-07, 9to5Mac](https://9to5mac.com/2026/05/07/anthropic-updates-claude-managed-agents-with-three-new-features/)) — agents review past sessions to self-improve; lead delegates to specialists. → [Issue #305](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/305) (PM eval, memory tightening). *PM. methodology-signal.*
- **Claude Code May 2026 update** ([changelog](https://code.claude.com/docs/en/changelog)) — `claude project purge`, smarter /model picker, plugin URLs, OAuth fixes. → No action; Dev-tooling. *PM. tooling-signal.*

## 2026-05-07

### 10:03 UTC — Editor: Developer

- **HISAT2 2.2.2** ([release](https://github.com/DaehwanKimLab/hisat2/releases), 2026-01-27) — adds repeat-read alignment mode (1 tagged alignment per repeat read). → [Issue #297](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/297) (chr22 trade-off audit). *Developer. pipeline-relevant.*
- **AFDB legacy API retires June 2026** ([NAR 2026](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkaf1226/8340156)) — TCRdock self-contained, no API calls. → No action; audited via [Issue #180](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/180). *Developer. pipeline-relevant.*

### 09:09 UTC — Editor: PM

- **Improved GitHub Issues search GA** ([May 2026 Changelog](https://github.blog/changelog/)) — NL + structured filters stable. → [Issue #294](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/294) (PM eval). *PM. tooling-relevant.*
- **Copilot cloud agent sessions in Issues/Projects** ([2026-04-23 Changelog](https://github.blog/changelog/2026-04-23-view-and-manage-agent-sessions-from-issues-and-projects/)) — view/steer from board, default-on. → [Issue #295](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/295) (PM trial, free tier). *PM. tooling-relevant.*
- **VS 2026 Cloud Agent integration** (May 2026) — remote-agent button in VS UI; shares Copilot SKU. → No action. *PM. tooling-signal.*
- **"Coherence Through Orchestration, Not Autonomy"** (Mason, [Jan 2026](https://mikemason.ca/writing/ai-coding-agents-jan-2026/)) — reinforces our `multi-role, not multi-agent` framing. → No action. *PM. methodology-signal.*

## 2026-05-06

### 10:34 UTC — Editor: PM

- **GitHub Copilot code review billing change** (effective 2026-06-01) — Copilot reviews on private repos consume Actions minutes + AI Credits. → No action; public repo. *PM. tooling-signal.*
- **CVE-2026-3854** (GitHub, May 2026) — high-severity flaw: authenticated push-access user can trigger RCE via crafted git push. → Already in [glossary](research/glossary.md) as RCE worked example. *PM. security-signal.*

### 08:46 UTC — Editor: Scientist

- **Kwok et al.** ([Nature 2025](https://www.nature.com/articles/s41586-024-08552-0)) — recurrent *GNAS*/*RPL22* splicing → public splice neoantigens + cross-patient TCRs. → [Issue #280](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/280) (Sci, DISCUSSION). *Scientist. portfolio-differentiator.*

## 2026-05-04

### 10:59 UTC — Editor: Developer

- **NeoGuider** ([XuegongLab](https://github.com/XuegongLab/neoguider)) — end-to-end neoepitope pipeline, splice-variant support. → [Issue #258](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/258) (Sci eval). *Both. portfolio-differentiator.*

## 2026-05-03

### 14:07 UTC — Editor: Scientist

- **Autogene cevumeran 6-year follow-up** (AACR 2026 abstract, Balachandran/MSK) — ~90% of immune-responders alive at ~6 years post-vaccination in Phase 1 PDAC trial; durable CD8+ T-cell repertoire. → No issue (DISCUSSION reference for clinical translation breadth). *Scientist. clinical-validation.*
- **Rojas et al. — original Phase 1 autogene cevumeran** (Nature 2023, foundational) — first PDAC personalized mRNA neoantigen vaccine trial; 8/16 immune responders. → No issue (foundational INTRODUCTION + DISCUSSION reference). *Scientist. clinical-validation.*
- **Nature 2025 follow-up — long-lived CD8+ T cells from autogene cevumeran** — 7.7-year mean estimated T-cell clone lifespan; RFS not reached for responders vs 13.4-mo non-responders. → No issue (DISCUSSION reference for vaccine durability mechanism). *Scientist. clinical-validation.*
- **Personalized Cancer Vaccines clinical trial pipeline review** (APJCO 2026) — field-wide landscape including TG4050 head/neck (no relapses at 16.2-mo median, Arm A). → No issue (field-landscape DISCUSSION reference). *Scientist. methodology-signal.*

### 10:49 UTC — Editor: Developer

- **PyTorch 2.7 + CUDA 12.6 = last P100-supporting combo** (2026, dev-discuss) — Maxwell/Pascal/Volta described as "feature-complete with no further enhancements planned"; PyTorch 2.8 dropped Pascal kernels in cu128/cu129 builds. → CLAUDE.md note clarifying the `torch<2.5` pin is permanent on this hardware. *Developer. pipeline-relevant.*
- **Snakemake `googlebatch` executor plugin** — official replacement for the deprecated Google Life Sciences executor. → Already tracked in [Issue #66](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/66) (logged here so it doesn't re-surface in future briefings). *Developer. pipeline-relevant.*
- **Snakemake 9.x deprecates `--use-conda`** in favour of `--software-deployment-method conda` — affects every snakemake call site in CLAUDE.md, `run_cloud_gpu.sh`, `setup_local.sh`. → Comment on [Issue #200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200) — bundle into the 9.x migration. *Developer. pipeline-relevant.*

## 2026-05-02

### 10:04 UTC — Editor: Scientist

- **AI predicted TCR-pMHC structures differentiate immune interactions** (2026-02, bioRxiv) — AlphaFold2 most consistent for TCR-pMHC multimers; structural > sequence features for binding; non-binders less stable in MD. → Reinforces [Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) (HERMES eval). *Scientist. portfolio-differentiator.*
- **t2pmhc** (2026-02, bioRxiv) — Structure-informed GNN for TCR-pMHC binding; mode (c) — needs full predicted complex. → [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) (hybrid models eval). *Scientist. portfolio-differentiator.*
- **TCRLens** (2026-01, Bioinformatics Advances) — Structure-aware EGNN with VAE-GAN augmentation; multi-scale graphs over 5 interface zones. → [Issue #236](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/236) (hybrid models eval). *Scientist. portfolio-differentiator.*

### 09:41 UTC — Editor: PM

- **GitHub MCP Server — Projects tools** (2026-01-28) — official MCP exposes `project_v2` mutations (Status/Priority/Size/Target date) at lower token cost than raw `gh api graphql`. → [Issue #234](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/234) (PM eval, migrate hand-rolled GraphQL). *PM. tooling-relevant.*
- **OpenAI Symphony** (2026-04, late April) — enterprise multi-agent orchestration; PM boards as control plane for coding agents. → No action; pattern-confirmation only for our PM/Sci/Dev split, product itself targets corporate audiences (not solo/research). *PM. methodology-signal.*
- **Anthropic 2026 Agentic Coding Trends Report** (2026, PDF) — industry data on multi-agent coding workflows. → [Issue #235](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/235) (PM skim, 5 PM-practices cross-checks: role decomposition, coordination, memory, scope discipline, handoff). *PM. methodology-signal.*

## 2026-05-01

- **uv 0.11.8 / ruff 0.15.7** (2026-05) — Astral toolchain dominant, 126M+/mo downloads. → No action; revisit if CI linting added. *Developer. industry-standard.*
- **HERMES** (2026-02, bioRxiv) — structure-based TCR-pMHC ML, zero domain training, 0.72 corr. → [Issue #218](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/218) (Scientist eval). *Both. portfolio-differentiator.*
- **MHCflurry 2.2.0** (2026-03-26) — NumPy 2.0 compat, drops `numpy<2.0` pin. → No action; torch 2.4.1 already compat. *Developer. pipeline-relevant.*
- **Snakemake 9.19.0** (2026-05) — 9.x steady; we're on 8.x. → [Issue #200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200). *Developer. pipeline-relevant.*

## 2026-04-30

- **AlphaGenome** (DeepMind 2025) — DNA→splicing/ATAC/Hi-C/CAGE prediction, 1Mbp context. → [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) (normal filter rethink). *Both. portfolio-differentiator.*
- **ImmSET** (2026-03) — TCR-pMHC specificity prediction, no structural modelling. → [Issue #201](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/201) (Scientist eval). *Scientist. portfolio-differentiator.*
- **Snakemake 9.x** (2025) — breaking changes vs 8.x (executor plugin API). → [Issue #200](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/200). *Developer. pipeline-relevant.*
