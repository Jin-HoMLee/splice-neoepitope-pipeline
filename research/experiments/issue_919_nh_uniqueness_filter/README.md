# Issue #919 - does the NH-uniqueness filter remove *spurious* junctions, or just *arbitrary* ones?

**Parent issue:** [#919 - opt-in NH-uniqueness prefilter for junction extraction](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/919)
**Shipped in:** [PR #1113](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1113)
**Status:** complete. **Result is a null**, and the null is the finding.
**Decision (AC 9):** keep the filter **opt-in, default off**.

## Goal

The filter (`alignment.uniqueness_filter.enabled`, a `samtools view -e '[NH]==1'` BAM prefilter feeding regtools) rests on a specific mechanism claim: HISAT2 emits one arbitrary copy of a multimapped read, so spliced reads over repeats produce junction calls at arbitrary repeat copies, and dropping `NH>1` alignments should remove those.

That claim is falsifiable. If the junctions the filter removes are **not** enriched for repeat overlap relative to the ones it keeps, then it is removing something else, and the mechanism story is wrong even when the counts look plausible.

## Result

**Within the unannotated pool - the pool the filter actually draws from - enrichment is 0.98x (tumor) and 1.00x (matched normal). There is none.**

| stratum (tumor) | junctions | splice site in a repeat |
|---|---|---|
| lost to the filter, annotated | 11 | 0.0% |
| lost to the filter, **unannotated** | 256 | **96.9%** |
| retained, annotated | 276 | 3.3% |
| retained, **unannotated** | 1,329 | **99.0%** |

An **unstratified** comparison instead reads 1.13x, which looks like a mild win and is a **composition artifact**: the lost set is annotation-poor and annotated junctions are repeat-poor. Read depth says the same thing (85.8% of lost are single-read vs 87.8% of retained). By both probes available, what the filter removes is indistinguishable from what it keeps.

## Why - the load-bearing finding

**`NH` is computed relative to the index, so chr22 test data cannot validate this filter in principle.**

With a chr22-only index, a read from a repeat copy on another chromosome has nowhere else to map. It aligns **uniquely** to a chr22 copy and is tagged `NH=1`. A single-chromosome index therefore converts genome-wide multimappers into apparent unique mappers, deflating `NH` and blinding the filter to exactly the population it exists to catch.

The data shows the damage. Only 8% of reads map at all, and the unannotated junction pool is repeat-**saturated**:

| reference rate | splice site in a repeat |
|---|---|
| random mappable chr22 position | 51.9% |
| GENCODE-annotated splice site | 17.5% |
| **unannotated junctions here** | **99.0%** |

There is no headroom for the filter to be "enriched" for repeats because essentially everything in the pool already is. This is not a weak test of the filter; it is **not a test of the filter**. No refinement of the A/B fixes it.

Consequence: [#1095](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1095) (whole-genome run) is promoted from a completeness check to **the only run that can decide the question**, and [#1116](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1116) carries the index-independent alternative (a Portcullis-style anchor-vs-intron Hamming test, which *would* work on this fixture).

## Outputs

| file | what it is |
|---|---|
| `outputs/raw_junctions.{tumor,normal}.filter_{off,on}.tsv` | the A/B inputs: junction sets from the chr22 run with the knob off and on |
| `outputs/junction_repeat_categorization.{tumor,normal}.tsv` | per-junction fate (lost/retained/gained) + repeat overlap at each splice site |
| `outputs/repeat_overlap_report.{tumor,normal}.md` | the rendered stratified report |

All committed: ~380 KB total, well inside the <10 MB band, and offline-regenerable from committed inputs plus the `rmsk` fetch below.

## Slide deck

`slides.qmd` (14 slides, experiment tier - see `docs/research_artifact_conventions.md`). Renders to reveal.js:

```bash
cd research/experiments/issue_919_nh_uniqueness_filter/
quarto render slides.qmd --to revealjs        # -> slides.html (gitignored)
```

Figures are **not** hand-drawn numbers: `figures/_regenerate_figures.py` recomputes every plotted value from `outputs/junction_repeat_categorization.*.tsv` plus the `rmsk` BED, so the analysis stays canonical and the deck can never drift from it. Re-run it after re-running the analysis:

```bash
research/.venv/bin/python research/experiments/issue_919_nh_uniqueness_filter/figures/_regenerate_figures.py
```

(The research venv, not the snakemake env - matplotlib lives there, per `research/requirements.txt`.)

## Reproducing

```bash
# 1. the RepeatMasker track (production rule; lands in gitignored references/)
snakemake --cores 1 -- references/rmsk/hg38/rmsk.chr22.bed

# 2. the A/B junction sets (already committed under outputs/; regenerate with)
snakemake --cores 4 --use-conda --configfile config/test_config.yaml -f -- \
  results/patient_001_test/alignment/SRR9143066_test/raw_junctions.tsv          # knob off
#   ... and again with alignment.uniqueness_filter.enabled: true for the ON set

# 3. the analysis
python research/experiments/issue_919_nh_uniqueness_filter/junction_repeat_overlap.py \
  --junctions-off outputs/raw_junctions.tumor.filter_off.tsv \
  --junctions-on  outputs/raw_junctions.tumor.filter_on.tsv \
  --rmsk references/rmsk/hg38/rmsk.chr22.bed \
  --annotated-bed resources/test/chr22_reference_junctions.bed \
  --label "tumor (chr22)"
```

**`--annotated-bed` is not optional in spirit.** Without it the tool reports only the unstratified enrichment - the 1.13x artifact this experiment exists to debunk - and it will say so loudly rather than let that number stand.

## Cross-experiment deps

- **Reads from production:** `resources/test/chr22_reference_junctions.bed` (the annotated-intron ground truth) and `workflow/scripts/filter_junctions.py` (`_load_reference_junctions`, reused so "annotated" has exactly one definition in the repo).
- **Reads from `references/`:** `references/rmsk/hg38/rmsk.chr22.bed`, fetched by the `download_rmsk_chrom` production rule (gitignored, 79,521 repeats).
- **Downstream consumers:** [#1095](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1095) (same analysis, whole-genome index) and [#1116](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1116). Per the promotion rule in `docs/research_artifact_conventions.md`, `junction_repeat_overlap.py` moves to `research/experiments/_shared/` only once a second consumer actually materializes - not before.

## Layout note

This experiment has no `notebook.ipynb`. The analysis is a **tested CLI tool** (`junction_repeat_overlap.py` + `tests/`, the same shape as `issue_547_immunogenicity_calibration/`), because it is re-run across issues rather than being a one-shot exploration, and because the numbers it produces are the evidence base for a production default. Per [Cookiecutter Data Science](https://cookiecutter-data-science.drivendata.org/opinions/), source code beats a notebook here: it "is more portable, can be tested more easily, and is easier to code review."
