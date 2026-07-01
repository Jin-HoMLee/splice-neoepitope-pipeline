# Deep-research runs - immunogenicity-benchmark arc

Durable record of the multi-agent deep-research sweeps run for the splice-immunogenicity benchmark ([#680](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/680) registry, [#736](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/736) scoring harness).
Each run was a Workflow fan-out (many parallel sub-agents: discovery lenses / design stances / catalog families, then adversarial verification, then synthesis).
These briefs are the synthesized outputs; the actionable findings from each are also routed onto the board (see "Feeds" per run).

## Why these live here
The raw per-run Workflow output files (`tasks/<id>.output`) are **session-ephemeral** - they live in a scratch directory tied to one Claude Code session and are not durable.
These committed briefs are the preserved distillation.
Where a brief cites a `tasks/*.output` path or a `wf_*` / `w*` run ID, that is provenance only - the underlying file is not in the repo.

## Runs

| # | Date | Topic | Scale | Honest headline |
|---|------|-------|-------|-----------------|
| 1 | 2026-06-30 | Registry enrichment (hunt positives + hard negatives) | 36 agents, ~3.0M tok | **0 new hard true-negatives** (field still n=1); the win was a 5-peptide HLA-A\*24:02 cluster (Oka 2021, Tg-mouse surrogate) denting the A2 monoculture |
| 2 | 2026-06-30 | #736 scoring-harness design | 27 agents, ~1.6M tok | Ship #736 as a **no-scalar ranking/enrichment benchmark on A\*02:01** with a power-ledger red-light headline; specificity is structurally unmeasurable at n=1 (winner: uncertainty-first) |
| 3 | 2026-07-01 | Baseline-predictor landscape for #736 | 63 agents, ~0.56M tok (partial; usage-limit) | Core-5 CPU-only baseline set (MHCflurry, NetMHCpan-4.1, BigMHC_IM, PRIME 2.1, IEDB-Calis); leakage-circularity is the sharpest threat; a DTU publish-license clause touches 2 of the 5 |

## Files
- `run1_registry_enrichment_2026-06-30.md` - the 16 confirmed registry-ready rows, 4 needs-human-access PDFs, rejects, completeness critique.
- `run2_scoring_design_2026-06-30.md` - judged design ranking, recommended design, forbidden-claims manifest, 7 open questions, 11 implementation steps.
- `run3_predictor_landscape_2026-07-01.md` - ranked baseline set, splice-admissibility + availability table (19 predictors), methodology implications, decision forks.
- `run3_predictor_dossiers.json` - raw structured dossiers behind run #3 (per-predictor inputs, score semantics, availability, verify verdicts).
- `workflow_scripts/` - the Workflow orchestration scripts for runs #2 and #3 (method transparency; paths inside are session-specific).

## Where each run feeds on the board
- **Run #1** -> [#911](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/911) (hard-negative sourcing: the neoTST preprints), [#839](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/839) (A\*24:02 rebalance feeder), [#681](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/681) (presentation decoys + Courcelles), [#817](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/817) (Zhao 2025 sequence recovery).
- **Run #2** -> [#736](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/736) (the design it commits).
- **Run #3** -> [#736](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/736) (baseline set) + [#926](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/926) (leakage-annotation prerequisite).

## Strategic through-line
More literature searching will not move the hard-negative constraint (run #1 proved the field is dry); the remaining yield is PDF-fetch + wet-lab, not more mining.
So #736 must ship honest-about-being-underpowered (run #2), and the baseline comparison must be leakage-guarded because our positives partly come from the baselines' own training corpora (run #3).
