# Issue #1161 - explicit junction depth floor

This artifact records the local, open-data characterization that motivated the
replacement of the legacy per-file mean gate. It uses only committed chr22
outputs from Issue #1162 and is reproducible without network access, conda, or
patient-specific restricted data:

```bash
python3 research/experiments/issue_1161_junction_depth/analyze_depth.py
```

## Delivered rule

`filtering.min_junction_reads` is now an explicit integer floor applied to both
tumor and every matched-normal junction file. A row is retained when
`mapped_reads >= min_junction_reads`. The production default is **10**. A
5-read floor is an explicit supported choice for shallow/local panels and is
selected in `config/test_config.yaml` for the 500K-read chr22 smoke fixture.
There is no per-file floating threshold. The raw `mean_reads` value remains a
descriptive distribution statistic, and `min_junction_reads_threshold` is
recorded in the filter stats output.

Using one floor on both arms preserves the existing origin hierarchy and
`normal_source` provenance for retained normal matches. Matched-normal evidence
below the floor is intentionally not a `normal_shared` match. The shared rule
also prevents an upstream read-removal filter from changing the depth gate
itself; any origin change caused by changed normal evidence remains explicit in
the retained `normal_source` semantics. The GTEx and annotated checks remain
downstream of depth filtering, in their existing priority order.

## Field-practice rationale

The default of 10 is not fitted to this fixture. It is the explicit
junction-coverage cutoff used by pVACsplice (`--junction-score`, default 10):
<https://pvactools.readthedocs.io/en/stable/pvacsplice/run.html>. Open SNAF
API documentation exposes a higher tumor cutoff of 20 and a normal cutoff of 5
(<https://snaf.readthedocs.io/en/latest/api.html>); IRIS exposes explicit
junction read-coverage cutoffs in its open screening implementation
(<https://github.com/Xinglab/IRIS/blob/master/IRIS/IRIS_screening_sjc.py>).
These precedents support a declared integer floor rather than a
right-tail-derived mean. Ten is the conservative middle choice for a general
pipeline default. Five remains available when the operator deliberately
accepts the shallow-panel tradeoff.

This is a depth-quality convention, not experimental or clinical validation.
The local fixture is a 500K-read chr22 smoke panel, not a production cohort.

## Measured depth sensitivity

The legacy gate is reproduced from the committed Issue #1162 classified TSVs.
The HISAT2 panel has 1,872 raw tumor junctions, mean 1.205, and an effective
`reads >= 2` gate, retaining 234 rows before origin classification. The STAR
panel has 3,379 raw junctions, mean 2.111, and an effective `reads >= 3` gate.
The STAR uniqueness-on control has mean 1.488 and would move the legacy gate
back to `reads >= 2`, demonstrating that removing upstream rows changes the
threshold itself. Its raw tumor/normal files are committed, but no matching
reference/GTEx-classified TSV is committed, so candidate and presenter counts
for that control are not estimable.

Candidate counts below are `tumor_exclusive` rows after the matched-normal and
GTEx origin filters. Explicit rows at or above the legacy effective threshold
are exact because all legacy-retained rows are committed. The analysis skips a
lower-than-legacy threshold when the discarded rows lack origin provenance.

| Local panel | Metric | Legacy mean | Explicit 5 | Explicit 10 | Explicit 20 |
|---|---|---:|---:|---:|---:|
| HISAT2 chr22 | retained junctions | 234 | 12 | 4 | 1 |
|  | tumor-exclusive candidates | 122 | 1 | 0 | 0 |
| STAR chr22 | retained junctions | 423 | 79 | 15 | 6 |
|  | tumor-exclusive candidates | 319 | 45 | 7 | 2 |

For the HISAT2 panel, the 5-read option places 1,860 of 1,872 raw rows below
the floor and leaves three unannotated classified rows, including one
tumor-exclusive junction. The default 10 places 1,868 rows below the floor and
removes all tumor-exclusive rows. For STAR, the 5-read option retains 79
junctions and 45 candidates, while the default retains 15 junctions and 7
candidates. The apparent severity is an honest property of the shallow 500K-read
fixture and is why the local test configuration explicitly selects 5 rather
than silently inheriting production 10.

The output TSV also reports raw rows, means, below-floor rows, retained rows,
annotated rows, retained unannotated origin counts, and threshold for each
panel and option. Annotated rows are absent from the committed legacy
classified TSV, so their count is inferred as retained rows minus classified
rows; no annotated row is reclassified as a candidate. The artifact can be
regenerated from tracked inputs with the command above.

## Predicted-presenter coverage and limitations

No committed MHCflurry prediction TSV is available for the Issue #1162 STAR or
HISAT2 classified panels, and the local reference FASTA/contig outputs are
ignored run products. Therefore before/after predicted-presenter counts are
**not estimable from the available local/open artifacts** without rerunning the
pipeline's alignment, contig, peptide, and CPU MHCflurry stages. The
`predicted_presenters` field is intentionally `not_available`, not zero.

A historical lab-notebook entry reports 4,045 unique peptides and 147 strong
plus 401 weak presenters for an earlier chr22 run, but it cannot be joined to
the committed Issue #1162 panel rows or reassigned to the 5/10 thresholds. It
is context only, not a before/after result. No claim of experimental or
clinical validation is made.
