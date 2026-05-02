# News Log

Shared log of news surfaced in morning briefings. Prevents repeating the same item across sessions and roles (PM, Scientist, Developer all log here).

**Format:**
- `## YYYY-MM-DD` — date-level section (one per day; new dates at the top)
- `### HH:MM UTC — Editor: <Role>` — time + editor attribution (one per session per day)
- `- **Item** (date) — keywords. → action. *Role. Signal.*` — one bullet per item, nested under the time section

Historical entries (before 2026-05-02) pre-date the time/editor sub-heading convention and are kept as-is.

---

## 2026-05-02

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
