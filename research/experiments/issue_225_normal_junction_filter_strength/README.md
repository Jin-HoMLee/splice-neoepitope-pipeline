# Issue #225 — Normal-junction filter strength (chr22, patient_001)

**Issue:** [#225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225)
**Parent:** [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) (Experiment 3)
**Status:** in progress
**Design spec:** [docs/superpowers/specs/2026-05-21-issue-225-normal-junction-filter-strength-design.md](../../../docs/superpowers/specs/2026-05-21-issue-225-normal-junction-filter-strength-design.md)

## Goal

Quantify the marginal filtering value of AlphaGenome on patient_001's chr22 tumor junction set, alongside matched-normal and a Snaptron-derived chr22 GTEx pan-tissue proxy. Output the decision-rule numbers for the Exp 3 row of #203's decision table.

## Inputs

- **Tumor:** `results/patient_001_test/alignment/SRR9143066_test/junctions.tsv` (produced by test pipeline; not tracked in git)
- **Matched-normal:** `results/patient_001_test/alignment/SRR9143065_test/junctions.tsv` (produced for #224)
- **AlphaGenome predictions:** `research/notebooks/issue_224_alphagenome_exp1_outputs/chr22_stomach_predicted_junctions.parquet` *(cross-experiment dep — single consumer)*
- **GTEx panel:** built fresh from Snaptron hg38 GTEx v2 endpoint; cached locally

## Outputs

- `outputs/chr22_gtex_panel.parquet` — Snaptron pan-tissue chr22 panel (≥1 sample inclusion)
- `outputs/filter_overlap_table.tsv` — 3-way overlap + unique contributions
- `outputs/filter_venn_chr22.png` — 3-way Venn diagram

## Cross-experiment deps

- Reads `research/notebooks/issue_224_alphagenome_exp1_outputs/chr22_stomach_predicted_junctions.parquet` (single consumer; explicit path reference per the cross-experiment-sharing convention in CLAUDE.md).

## Caveats

- **Snaptron proxy ≠ #211 production panel.** When [Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211) lands, §2(c) of the notebook re-runs against the production GTEx panel.
- **chr22 only.** Test-config scope; full-genome scale-up is out of scope for #225.
- **Experiment 2 (germline-aware AG) deferred** to [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) pending WGS acquisition.

## Conda env

`splice-neoepitope-alphagenome` (defined in `workflow/envs/alphagenome.yaml`). Includes `matplotlib-venn` added by Task 1 of the implementation plan.
