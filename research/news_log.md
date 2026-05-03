# News Log

Shared log of news surfaced in morning briefings. Prevents repeating the same item across sessions and roles (PM, Scientist, Developer all log here).

**Format:**
- `## YYYY-MM-DD` — date-level section (one per day; new dates at the top)
- `### HH:MM UTC — Editor: <Role>` — time + editor attribution (one per session per day)
- `- **Item** (date) — keywords. → action. *Role. Signal.*` — one bullet per item, nested under the time section

Historical entries (before 2026-05-02) pre-date the time/editor sub-heading convention and are kept as-is.

---

## 2026-05-03

### 14:07 UTC — Editor: Scientist

- **Autogene cevumeran 6-year follow-up** (AACR 2026 abstract, Balachandran/MSK) — ~90% of immune-responders alive at ~6 years post-vaccination in Phase 1 PDAC trial; durable CD8+ T-cell repertoire. → No issue (DISCUSSION reference for clinical translation breadth). *Scientist. clinical-validation.*
- **Rojas et al. — original Phase 1 autogene cevumeran** (Nature 2023, foundational) — first PDAC personalized mRNA neoantigen vaccine trial; 8/16 immune responders. → No issue (foundational INTRODUCTION + DISCUSSION reference). *Scientist. clinical-validation.*
- **Nature 2025 follow-up — long-lived CD8+ T cells from autogene cevumeran** — 7.7-year mean estimated T-cell clone lifespan; RFS not reached for responders vs 13.4-mo non-responders. → No issue (DISCUSSION reference for vaccine durability mechanism). *Scientist. clinical-validation.*
- **Personalized Cancer Vaccines clinical trial pipeline review** (APJCO 2026) — field-wide landscape including TG4050 head/neck (no relapses at 16.2-mo median, Arm A). → No issue (field-landscape DISCUSSION reference). *Scientist. methodology-signal.*

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
