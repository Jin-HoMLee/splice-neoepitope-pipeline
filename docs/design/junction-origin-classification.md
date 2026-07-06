# Junction origin classification

*Design explanation. Migrated verbatim (punctuation normalized to house style) from `CLAUDE.md` for [Issue #777](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/777).*

Normal samples are used to filter tumor junctions at the junction level, not at the prediction level. The hierarchy:

```
all junctions
  └─ annotated        (in GENCODE)            → discard
  └─ unannotated      (not in GENCODE)
       ├─ normal_shared  (also in normal)  → discard (kept in TSV for reference)
       └─ tumor_exclusive    (absent in normal) → neoepitope prediction
```

This is the clinically correct approach: a junction present in matched normal tissue is not tumor-specific and should not be a neoepitope target. The Fisher's exact test (end-of-pipeline statistical comparison) was removed in favour of this upstream filtering step.

When no normal sample is present, all unannotated junctions are labeled `tumor_exclusive` with a warning - the pipeline still runs.
