# Issue #225 — Normal-junction filter strength on patient_001 (chr22) — design

**Issue:** [#225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) — *research(filter): comparative filter strength on patient_001 (matched-normal vs GTEx vs AlphaGenome)*
**Parent:** [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) — *rethink normal filtering with population panel + AlphaGenome fallback*
**Author:** Scientist
**Date:** 2026-05-21
**Status:** Design approved 2026-05-21; pending user review of this spec.

## Goal

Run Experiment 3 from #203: quantify the marginal filtering value of AlphaGenome (AG) on patient_001's chr22 tumor junction set, alongside matched-normal (MN) and a Snaptron-derived chr22 GTEx pan-tissue proxy panel. Output a 3-way overlap table, a Venn diagram, and the decision-rule numbers that determine whether AG is adopted as a 3rd always-on filter source (#203 decision table).

## Scope

- **In scope:** chr22 only (test-config harness); patient_001 only; three filters (MN, Snaptron-chr22-GTEx-proxy, AG-predicted-normal); 3-way set-overlap analysis; decision-rule outcome for the Exp 3 row of #203.
- **Out of scope:** Exp 2 (germline-aware AG, deferred to [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) pending WGS); production GTEx panel (deferred to [Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211)); whole-genome scale-up; tissue-matched (vs pan-tissue) GTEx.
- **Migration of existing notebooks:** a separate follow-up Issue, filed after #225 lands.

## Key design decisions

### D1. Per-Issue notebook, not appending to #224
The original #225 body said "Add a section to the Experiment 1+2 notebook (no separate notebook)" — written 2026-05-01 when Exp 1+2 were a single notebook. The 2026-05-16 audit carved Exp 2 out to [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381); #224 is now Exp 1 only. Folding Exp 3 in would re-couple work the audit deliberately split. #225 gets its own notebook. **The #225 body is updated to reflect this as part of the PR.**

### D2. New folder convention: `research/experiments/<issue_NNN>_<short-desc>/`
Today `research/notebooks/` mixes two species:
- Stable per-patient analyses (`patient_001_results.ipynb`, `patient_002_results.ipynb`) — canonical, manuscript-supporting, long-lived
- Per-issue experimental work (`issue_224_*`, `issue_299_*`) — scoped, time-bounded

Separating them mirrors `research/slides/issue_NNN_<short-desc>/slides.qmd`. Per-experiment artifacts (notebook + outputs + README) live together in a self-contained folder.

### D3. Cross-experiment data-sharing convention (new)
1. **Default:** each experiment owns its outputs in `<experiment>/outputs/`.
2. **Shared (≥2 consumers):** `research/experiments/_shared/`. Promote lazily — only when a 2nd consumer materializes. Filenames carry provenance (`gtex_panel_chr22_snaptron_v1.parquet`, not `gtex_panel.parquet`).
3. **Cross-experiment read, single consumer:** explicit path reference. Documented in the experiment's README under "Cross-experiment deps".
4. **Promoted to production:** moves to `resources/` when stable enough to be a pipeline input.

**Size guidance:**
- `< 10 MB` → check into git.
- `10–100 MB` → check into git with a regenerator script committed alongside.
- `> 100 MB` → keep out of git; store in `gs://splice-neoepitope-project/experiments/<issue>/`; commit a `data_manifest.yaml` in the experiment folder listing artifact paths + checksums + fetch commands.

**Bad practices to avoid (all sizes):** symlinks across experiments, copying artifacts, one experiment writing into another's `outputs/`.

This convention lands in CLAUDE.md as part of #225's PR.

### D4. GTEx axis = Snaptron chr22 proxy (not deferred, not blocking on #211)
[Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211) (real GTEx pan-tissue panel) is OPEN. Deferring the GTEx leg would leave #203's decision-rule incomplete. Building a chr22 proxy from Snaptron's hg38 GTEx endpoint reuses the same plumbing pattern as [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299) and stays consistent with #225's XS sizing. The proxy will be re-run with the production panel when #211 lands; only the §2(c) cell changes. **Caveat is flagged explicitly in the notebook + decision write-up.**

- Endpoint: `https://snaptron.cs.jhu.edu/gtexv2/snaptron` (hg38 GTEx v2)
- Coord scope: chr22 (matches the test-config harness)
- Tissue scope: pan-tissue (all 49 GTEx tissues — matches #203's vaccine-safety reasoning; "kept pan-tissue, not tissue-matched")
- Inclusion: sample-count ≥ 1 (most-conservative pan-tissue)

### D5. AlphaGenome threshold = F1-max from Exp 1
The AG parquet at `research/notebooks/issue_224_alphagenome_exp1_outputs/chr22_stomach_predicted_junctions.parquet` is a scored prediction set. The binary "in/out" call needs a threshold. Use the F1-max threshold computed in §4 of the #224 notebook (data-driven, defensible). If §4 didn't persist that value, recompute it inline at the top of §2(b) and print the chosen threshold + why.

### D6. Filter matching key = `(chrom, donor_0based, acceptor_0based_excl, strand)`, exact match
Same convention as `filter_junctions._parse_junction_id`. No strand collapsing. Matches the production pipeline's semantics so the comparative numbers translate cleanly to a production decision.

### D7. Tumor pipeline precondition is documented but not automated
Chr22 tumor junctions (`results/patient_001_test/alignment/SRR9143066_test/junctions.tsv`) are not yet produced. The notebook's §0 setup cell will fail-fast with the exact command to run if the file is missing. Not embedded in the notebook (running Snakemake from a notebook is fragile and slow to iterate on).

## Architecture

```
research/experiments/issue_225_normal_junction_filter_strength/
├── README.md                          # one-page: goal, parent issue link, status, outputs index, cross-experiment deps
├── notebook.ipynb                     # the analysis
└── outputs/
    ├── chr22_gtex_panel.parquet       # built by notebook §2(c), cached (Snaptron chr22 pan-tissue, ≥1 sample)
    ├── filter_overlap_table.tsv       # 3-way overlap, dropdown for #203 decision-rule writeup
    └── filter_venn_chr22.png          # 3-way Venn (also referenced from manuscript later)
```

**Inputs (referenced by path, not copied):**
- `results/patient_001_test/alignment/SRR9143066_test/junctions.tsv` (tumor; **precondition: run test pipeline for SRR9143066**)
- `results/patient_001_test/alignment/SRR9143065_test/junctions.tsv` (matched-normal; produced by #224)
- `research/notebooks/issue_224_alphagenome_exp1_outputs/chr22_stomach_predicted_junctions.parquet` (AG predictions; cross-experiment dep, single consumer)
- Snaptron hg38 GTEx endpoint (network)

## Components

The notebook has 7 sections after §0 setup:

| § | Purpose | Output |
|---|---|---|
| §0 | Setup — imports, paths, `REPO_ROOT`, path-existence sanity prints | (none) |
| §1 | Load chr22 tumor junctions (SRR9143066) via the same `load_pipeline_junctions` loader as #224 | `tumor` DataFrame |
| §2 | Load/build 3 filter sets: (a) matched-normal from `SRR9143065_test/junctions.tsv`, (b) AG predicted-normal from #224 parquet @ F1-max threshold, (c) Snaptron chr22 pan-tissue GTEx panel (hg38, ≥1 sample inclusion) | `mn_set`, `ag_set`, `gtex_set` + cached `outputs/chr22_gtex_panel.parquet` |
| §3 | Apply each filter separately: `tumor_caught_by_F = tumor ∩ F` (set-op on the matching key) | `caught_by_mn`, `caught_by_gtex`, `caught_by_ag` |
| §4 | Overlap analysis: pairwise + triple intersections, unique contributions per source | `outputs/filter_overlap_table.tsv` |
| §5 | 3-way Venn (`matplotlib-venn`) | `outputs/filter_venn_chr22.png` |
| §6 | Decision-rule outcome: render the #203 decision table with the computed numbers; flag which row triggers (adopt / fallback / no-go). Exp 2 delta row marked N/A (deferred to [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381)). | Markdown table in notebook |
| §7 | Caveats + forward-references: Snaptron proxy disclaimer, Exp 2 deferred note, "re-run §2(c) with production panel when [Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211) lands" | (none) |

## Data flow

```
results/patient_001_test/alignment/SRR9143066_test/junctions.tsv ──┐ (§1)
                                                                    │
results/patient_001_test/alignment/SRR9143065_test/junctions.tsv ──┤ (§2a)
                                                                    ├─► (§3) 3 filtered sets ─► (§4) overlap table
research/notebooks/issue_224_*/outputs/chr22_stomach_*.parquet ────┤ (§2b)                                  └─► (§5) Venn
                                                                    │                                       └─► (§6) decision rule
Snaptron hg38 GTEx endpoint  ─────────────────────────────────────┘ (§2c)
                                                                    │
                                                                    └─► outputs/chr22_gtex_panel.parquet (cached)
```

## Error handling

This is a notebook (interactive, not unattended). "Error handling" = sanity guards that fail fast with a clear message rather than producing a silently-wrong overlap table.

- **Path-existence guards (§0):** assert each input file exists; print which is missing with the fix command (e.g. "missing — run `snakemake --cores 4 --use-conda --configfile config/test_config.yaml`" for the tumor TSV).
- **Coord-system guards (§1 + §2a):** assert tumor + matched-normal junctions are chr22-only; assert non-zero junction counts (catches a failed alignment).
- **AG parquet schema guard (§2b):** assert expected columns exist (`chrom`, donor, acceptor, strand, score-or-similar); pick the F1-max threshold from §4 of the #224 notebook if persisted, else recompute inline with a clear print of the chosen threshold + why.
- **Snaptron query guard (§2c):** wrap the API call with a timeout + retry-once; assert non-empty result; assert chr22-only result; assert `(chrom, donor, acceptor, strand)` keys parse cleanly. **Cache only on success** — never write a partial parquet.
- **Set-arithmetic guards (§3):** assert each filter's `tumor ∩ F` is a subset of `tumor` (catches a coord-system mismatch where the intersection silently goes empty or explodes).
- **Decision-rule guard (§6):** if F1 (from §4 of #224) wasn't recomputable, raise rather than producing a half-rendered decision table. This is the headline output — partial = misleading.

## Testing

No new pytest. The notebook's correctness is verified inline by:

1. **Cell-level assertions** (see Error handling).
2. **Read-back checks** — after §4 writes the overlap TSV, reload it and verify counts match the in-memory values (catches silent dtype / serialization bugs).
3. **Caveats section (§7)** lists the assumptions a reviewer needs to verify: Snaptron sample-coverage threshold = 1, AG threshold = F1-max-from-Exp-1, GTEx pan-tissue not tissue-matched, chr22 only.

CI doesn't run notebooks today; the verification gate is "scientist re-executes the notebook end-to-end before PR" (standard per-experiment lab-notebook workflow).

## Deliverables

PR checklist (Test plan for the eventual PR):
- [ ] `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb` end-to-end clean run from a fresh kernel
- [ ] `outputs/chr22_gtex_panel.parquet` produced + cached
- [ ] `outputs/filter_overlap_table.tsv` matches notebook §4 in-memory values (read-back check)
- [ ] `outputs/filter_venn_chr22.png` rendered
- [ ] `README.md` in the experiment folder one-pager (goal, parent issue, status, outputs index, cross-experiment deps)
- [ ] CLAUDE.md updated with the `research/experiments/` convention (D2 + D3)
- [ ] #225 body updated to drop the stale "append to Exp 1+2 notebook" line
- [ ] Lab notebook entry per the shared always-in-effect rule
- [ ] #203 decision-rule outcome updated in the parent issue body (Exp 3 row)
- [ ] Decision-rule outcome triggers one of: adopt / fallback / no-go (#203 decision table)

## Open items / known cross-experiment deps

- AG parquet read-by-path coupling to `issue_224_*/outputs/`. Documented in the experiment README. Will be revisited when the migration Issue runs.
- The cross-experiment-sharing convention (D3) lands here in CLAUDE.md; future experiments inherit it without re-deciding.

## References

- Parent epic — [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) (decision rule + Exp 1/2/3 framing)
- Sibling — [Sub-Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) (Exp 1, closed; AG predictions cached)
- Sibling — [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) (Exp 2, deferred on WGS)
- Adjacent — [Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211) / [Issue #212](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/212) (production GTEx panel; #225 uses Snaptron proxy until #211 lands)
- Plumbing pattern — [Issue #299](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/299) (Snaptron query convention)
