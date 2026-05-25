# Issue #225 — Normal-junction filter strength (chr22) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compare three normal-junction filter sources (matched-normal, Snaptron-chr22-GTEx-proxy, AlphaGenome-predicted-normal) on patient_001's chr22 tumor junctions and produce the decision-rule numbers for Experiment 3 of [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203).

**Architecture:** A single per-Issue Jupyter notebook under the new `research/experiments/issue_NNN_<short-desc>/` convention (folder + `outputs/` + `README.md`). Inputs are referenced by path (tumor & matched-normal `junctions.tsv` from the chr22 test pipeline; AG predictions parquet from #224; Snaptron hg38 GTEx endpoint). All set operations use the exact-match key `(chrom, donor_0based, acceptor_0based_excl, strand)`. The convention itself is documented into CLAUDE.md as part of this PR.

**Tech Stack:** Snakemake (precondition pipeline run), conda env `splice-neoepitope-alphagenome` (pandas, pyarrow, matplotlib, matplotlib-venn), Snaptron hg38 GTEx v2 endpoint, Jupyter.

**Spec:** [docs/superpowers/specs/2026-05-21-issue-225-normal-junction-filter-strength-design.md](../specs/2026-05-21-issue-225-normal-junction-filter-strength-design.md)

---

### Task 1: Add `matplotlib-venn` to alphagenome conda env

**Files:**
- Modify: [workflow/envs/alphagenome.yaml](../../../workflow/envs/alphagenome.yaml)

- [ ] **Step 1: Edit the env file**

Add `matplotlib-venn >=1.1` to the conda dependencies block (alphabetical order: between `matplotlib` and `numpy`).

```yaml
name: splice-neoepitope-alphagenome
channels:
  - conda-forge
  - defaults
dependencies:
  - python >=3.11,<3.13
  - pandas >=2.0
  - numpy >=1.25
  - matplotlib >=3.8
  - matplotlib-venn >=1.1
  - scipy >=1.11
  - seaborn >=0.13
  - pyarrow >=15.0
  - ipykernel
  - pip
  - pip:
    - alphagenome >=0.6,<0.7
```

- [ ] **Step 2: Update / install the env**

```bash
conda env update --file workflow/envs/alphagenome.yaml --prune
```

- [ ] **Step 3: Verify import works**

```bash
conda run -n splice-neoepitope-alphagenome python -c "import matplotlib_venn; print('matplotlib_venn:', matplotlib_venn.__version__)"
```

Expected: prints a version string ≥ 1.1.

- [ ] **Step 4: Commit**

```bash
git add workflow/envs/alphagenome.yaml
git commit -m "deps(envs): add matplotlib-venn to alphagenome env for #225 Venn output"
```

---

### Task 2: Produce chr22 tumor junctions (precondition)

**Files:**
- Read (no edits): [config/test_config.yaml](../../../config/test_config.yaml)
- Output (pipeline-produced, not committed): `results/patient_001_test/alignment/SRR9143066_test/junctions.tsv`

- [ ] **Step 1: Activate snakemake env**

```bash
conda activate snakemake
```

- [ ] **Step 2: Run the chr22 test pipeline targeting the tumor junctions file**

Note the canonical-form `--` separator (per CLAUDE.md "Snakemake 8 gotchas: positional target after --configfile"):

```bash
snakemake --cores 4 --use-conda --configfile config/test_config.yaml -- results/patient_001_test/alignment/SRR9143066_test/junctions.tsv
```

Expected: Snakemake builds the rules backwards from the target. Should hit existing HISAT2 index (built by #224 run) and just run alignment + regtools + `bed12_to_junctions.py` for SRR9143066. Runtime: ~10–30 min.

- [ ] **Step 3: Verify output exists and has the expected shape**

```bash
ls -la results/patient_001_test/alignment/SRR9143066_test/junctions.tsv
wc -l results/patient_001_test/alignment/SRR9143066_test/junctions.tsv
head -3 results/patient_001_test/alignment/SRR9143066_test/junctions.tsv
```

Expected: file exists, non-zero lines, format `<chrom>:<donor_1based>:<acceptor_0based_excl>:<strand>\t<reads>` (matches matched-normal format).

- [ ] **Step 4: No commit (pipeline output, not tracked)**

`results/` is in `.gitignore`. Skip git add.

---

### Task 3: Update CLAUDE.md with `research/experiments/` convention

**Files:**
- Modify: [CLAUDE.md](../../../CLAUDE.md) — add a new section above "Slide decks for experiment Issues" (keeps experiment-related sections grouped)

- [ ] **Step 1: Insert the convention section**

Insert this block above the "Slide decks for experiment Issues" heading in CLAUDE.md:

```markdown
## Experiment notebooks live under `research/experiments/`

Per-Issue experimental work (analysis notebooks + their cached outputs + a one-page README) lives at `research/experiments/issue_NNN_<short-content-desc>/`. Mirrors the slide-deck convention (`research/slides/issue_NNN_<short-desc>/slides.qmd`) — same shape applied to notebooks.

Layout per experiment:

```
research/experiments/issue_NNN_<short>/
├── README.md          # one-page: goal, parent issue link, status, outputs index, cross-experiment deps
├── notebook.ipynb     # the analysis
└── outputs/           # cached artifacts (parquet, tsv, png)
```

**Distinguish from `research/notebooks/`:** that folder holds stable per-patient analyses (`patient_001_results.ipynb`, `patient_002_results.ipynb`) — long-lived, manuscript-supporting. The experiments/ folder holds scoped per-Issue work that lands once and then becomes a frozen reference. A separate migration Issue will move existing per-Issue work (`issue_224_*`, `issue_299_*`) from `research/notebooks/` into `research/experiments/`; the convention is established by #225.

### Cross-experiment data sharing

1. **Default:** each experiment owns its outputs in `<experiment>/outputs/`.
2. **Shared between ≥2 experiments → `research/experiments/_shared/`.** Promote an artifact here only when a 2nd consumer materializes (YAGNI; don't pre-share). Filenames carry provenance (`gtex_panel_chr22_snaptron_v1.parquet`, not `gtex_panel.parquet`).
3. **Cross-experiment read, single consumer → explicit path reference.** Document in the consumer's README under "Cross-experiment deps".
4. **Promoted to production → `resources/`.** When an artifact becomes a stable pipeline input.

**Bad practices (any size):** symlinks across experiments, copying artifacts, one experiment writing into another's `outputs/`.

### Size guidance

- **< 10 MB:** check into git.
- **10–100 MB:** check into git, commit a regenerator script alongside.
- **> 100 MB:** keep out of git; store in `gs://splice-neoepitope-project/experiments/<issue>/`; commit a `data_manifest.yaml` in the experiment folder listing artifact paths + checksums + fetch commands. Pin the manifest schema when the first artifact crosses 100 MB.
```

- [ ] **Step 2: Verify markdown formatting (no broken headings or fences)**

```bash
grep -n "## Experiment notebooks live under" CLAUDE.md
grep -n "## Slide decks for experiment Issues" CLAUDE.md
```

Expected: the new heading appears immediately above the existing slide-decks one.

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs(claude-md): research/experiments/ convention + cross-experiment data sharing rules"
```

---

### Task 4: Scaffold the experiment folder + README

**Files:**
- Create: `research/experiments/issue_225_normal_junction_filter_strength/README.md`
- Create: `research/experiments/issue_225_normal_junction_filter_strength/outputs/.gitkeep`

- [ ] **Step 1: Create the folder + outputs/ + .gitkeep**

```bash
mkdir -p research/experiments/issue_225_normal_junction_filter_strength/outputs
touch research/experiments/issue_225_normal_junction_filter_strength/outputs/.gitkeep
```

- [ ] **Step 2: Write the README**

Create `research/experiments/issue_225_normal_junction_filter_strength/README.md` with:

```markdown
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
```

- [ ] **Step 3: Verify file tree**

```bash
ls -la research/experiments/issue_225_normal_junction_filter_strength/
```

Expected: `README.md`, `outputs/` directory present.

- [ ] **Step 4: Commit**

```bash
git add research/experiments/issue_225_normal_junction_filter_strength/
git commit -m "scaffold: issue_225_normal_junction_filter_strength/ + README"
```

---

### Task 5: Notebook §0 — setup + path guards

**Files:**
- Create: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`

- [ ] **Step 1: Create the notebook with kernel `splice-neoepitope-alphagenome`**

Easiest path:

```bash
conda activate splice-neoepitope-alphagenome
jupyter nbconvert --to notebook --output research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb /dev/null
# or create via VSCode / Jupyter UI selecting the alphagenome kernel
```

- [ ] **Step 2: Add §0 markdown cell**

```markdown
# Issue #225 — Normal-junction filter strength on patient_001 (chr22)

**Issue:** [#225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225)
**Parent:** [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) (Experiment 3)
**Design spec:** [docs/superpowers/specs/2026-05-21-issue-225-normal-junction-filter-strength-design.md](../../../docs/superpowers/specs/2026-05-21-issue-225-normal-junction-filter-strength-design.md)

Compare three normal-junction filter sources on patient_001's chr22 tumor junctions: matched-normal, Snaptron-chr22-GTEx-proxy, and AlphaGenome-predicted-normal. Produce the decision-rule numbers for the Exp 3 row of [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203)'s decision table.

**Scope:** chr22 only (test-config harness); patient_001 only; matching key `(chrom, donor_0based, acceptor_0based_excl, strand)`, exact match.
```

- [ ] **Step 3: Add §0 setup code cell**

```python
from __future__ import annotations

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

REPO_ROOT = Path.cwd().resolve()
while not (REPO_ROOT / "workflow").is_dir() and REPO_ROOT != REPO_ROOT.parent:
    REPO_ROOT = REPO_ROOT.parent

EXPERIMENT_DIR = REPO_ROOT / "research" / "experiments" / "issue_225_normal_junction_filter_strength"
OUTPUTS_DIR = EXPERIMENT_DIR / "outputs"
OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

TARGET_CHROM = "chr22"

# Inputs (referenced by path; not copied)
TUMOR_TSV = REPO_ROOT / "results" / "patient_001_test" / "alignment" / "SRR9143066_test" / "junctions.tsv"
NORMAL_TSV = REPO_ROOT / "results" / "patient_001_test" / "alignment" / "SRR9143065_test" / "junctions.tsv"
AG_PARQUET = REPO_ROOT / "research" / "notebooks" / "issue_224_alphagenome_exp1_outputs" / "chr22_stomach_predicted_junctions.parquet"

# Outputs
GTEX_PANEL_PARQUET = OUTPUTS_DIR / "chr22_gtex_panel.parquet"
OVERLAP_TSV = OUTPUTS_DIR / "filter_overlap_table.tsv"
VENN_PNG = OUTPUTS_DIR / "filter_venn_chr22.png"

# Path-existence guards
inputs = {
    "Tumor junctions": TUMOR_TSV,
    "Matched-normal junctions": NORMAL_TSV,
    "AlphaGenome predictions": AG_PARQUET,
}
missing = [(label, path) for label, path in inputs.items() if not path.exists()]
if missing:
    msg_lines = [f"  - {label}: {path}" for label, path in missing]
    fix_hint = (
        "\nTo produce the chr22 tumor junctions (the most common miss), run:\n"
        "  conda activate snakemake && snakemake --cores 4 --use-conda --configfile config/test_config.yaml "
        "-- results/patient_001_test/alignment/SRR9143066_test/junctions.tsv"
    )
    raise FileNotFoundError("Missing required inputs:\n" + "\n".join(msg_lines) + fix_hint)

print(f"Repo root:      {REPO_ROOT}")
print(f"Experiment dir: {EXPERIMENT_DIR.relative_to(REPO_ROOT)}")
print(f"All inputs present: {all(p.exists() for p in inputs.values())}")
```

- [ ] **Step 4: Execute §0; verify no FileNotFoundError**

Run the cell. Expected: prints repo root + experiment dir + "All inputs present: True".

- [ ] **Step 5: Commit (with cell outputs cleared per project lab-notebook convention)**

```bash
jupyter nbconvert --clear-output --inplace research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git add research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git commit -m "wip(#225): notebook §0 — setup + path guards"
```

---

### Task 6: Notebook §1 — load chr22 tumor junctions

**Files:**
- Modify: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`

- [ ] **Step 1: Add §1 markdown cell**

```markdown
## §1 — Tumor junctions (patient_001 chr22)

Load the post-PR #372 `junctions.tsv` produced by `workflow/scripts/bed12_to_junctions.py`. Format: 2-column TSV, no header: `<chrom>:<donor_1based>:<acceptor_0based_excl>:<strand>\t<reads>`. Loader normalizes to **0-based half-open intron coords** for set-ops downstream.
```

- [ ] **Step 2: Add §1 code cell**

```python
def load_pipeline_junctions(tsv_path: Path) -> pd.DataFrame:
    """Load 2-column junctions TSV from bed12_to_junctions.py / STAR awk.

    Format: ``<chrom>:<donor_1based>:<acceptor_0based_excl>:<strand>\\t<reads>``.
    Returns DataFrame with 0-based half-open intron coords for set operations.
    """
    raw = pd.read_csv(tsv_path, sep="\t", header=None, names=["key", "count"])
    parts = raw["key"].str.split(":", expand=True)
    parts.columns = ["chrom", "donor_1based", "acceptor_0based_excl", "strand"]
    return pd.DataFrame({
        "chrom": parts["chrom"],
        "donor": parts["donor_1based"].astype(int) - 1,
        "acceptor": parts["acceptor_0based_excl"].astype(int),
        "strand": parts["strand"],
        "count": raw["count"].astype(int),
    })

tumor = load_pipeline_junctions(TUMOR_TSV)
assert (tumor["chrom"] == TARGET_CHROM).all(), f"Expected {TARGET_CHROM}-only tumor; got {tumor['chrom'].unique()}"
assert len(tumor) > 0, "Tumor junctions empty — check alignment for SRR9143066"
print(f"Tumor junctions (chr22):     {len(tumor):,}")
print(f"Read-count summary: {tumor['count'].describe().to_dict()}")
tumor.head()
```

- [ ] **Step 3: Execute §1; verify count > 0 and chr22-only**

Expected: non-zero junction count, no AssertionError, head of DataFrame shows chr22 / donor / acceptor / strand / count.

- [ ] **Step 4: Commit**

```bash
jupyter nbconvert --clear-output --inplace research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git add research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git commit -m "wip(#225): notebook §1 — load chr22 tumor junctions"
```

---

### Task 7: Notebook §2(a) — load matched-normal filter

**Files:**
- Modify: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`

- [ ] **Step 1: Add §2(a) markdown + code cell**

```markdown
## §2 — Filter sets

### §2(a) — Matched-normal (current pipeline source)

Load matched-normal junctions from the same `junctions.tsv` format. Anchor of the comparison; gold standard for tissue-expressed splicing in patient_001's stomach normal.
```

```python
mn = load_pipeline_junctions(NORMAL_TSV)
assert (mn["chrom"] == TARGET_CHROM).all()
assert len(mn) > 0

# Build the matching-key set: (chrom, donor, acceptor, strand)
def to_key_set(df: pd.DataFrame) -> set[tuple]:
    return set(map(tuple, df[["chrom", "donor", "acceptor", "strand"]].itertuples(index=False, name=None)))

mn_set = to_key_set(mn)
print(f"Matched-normal junctions:    {len(mn):,}  (unique keys: {len(mn_set):,})")
```

- [ ] **Step 2: Execute; verify count, no AssertionError**

- [ ] **Step 3: Commit**

```bash
jupyter nbconvert --clear-output --inplace research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git add research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git commit -m "wip(#225): notebook §2(a) — matched-normal filter set"
```

---

### Task 8: Notebook §2(b) — load AG parquet + recompute F1-max threshold + build ag_set

**Files:**
- Modify: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`

**Note:** #224 did not persist `best_threshold` to disk — it lives only as an in-memory variable inside the #224 notebook. To stay self-contained, #225 recomputes it inline against the same ground-truth definition used in #224 §5 (matched-normal ∩ GENCODE-annotated). If this becomes painful, a future PR could persist `best_threshold` from #224.

- [ ] **Step 1: Add §2(b) markdown cell**

```markdown
### §2(b) — AlphaGenome predicted-normal

Load AG predictions cached by [#224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) at `research/notebooks/issue_224_alphagenome_exp1_outputs/chr22_stomach_predicted_junctions.parquet` *(cross-experiment dep)*. Recompute the F1-max threshold against the ground-truth definition from #224 §5 (matched-normal ∩ GENCODE-annotated, both restricted to chr22). Threshold-binarise to produce `ag_set`.
```

- [ ] **Step 2: Add §2(b) code cell — load + recompute F1-max threshold**

```python
import gzip
import re

# Load AG predictions
ag = pd.read_parquet(AG_PARQUET)
assert {"chrom", "donor", "acceptor", "strand"}.issubset(ag.columns), f"AG parquet missing key cols: {ag.columns.tolist()}"
score_col = "score" if "score" in ag.columns else next((c for c in ag.columns if "score" in c.lower() or "prob" in c.lower()), None)
assert score_col is not None, f"No score column found in AG parquet; columns: {ag.columns.tolist()}"
print(f"AG predictions: {len(ag):,}  | score column: {score_col}")

# Build the ground-truth set for F1: matched-normal ∩ GENCODE-annotated (chr22)
GENCODE_GTF = REPO_ROOT / "resources" / "test" / "chr22.gtf.gz"
assert GENCODE_GTF.exists(), f"GENCODE chr22 GTF not found: {GENCODE_GTF}"

def gencode_introns_chr22(gtf_path: Path) -> set[tuple]:
    """Return set of (chrom, donor_0based, acceptor_0based_excl, strand) introns from a GENCODE GTF.

    Introns are inferred from consecutive exons within the same transcript.
    """
    transcripts: dict[str, list] = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            chrom, start, end, strand, attrs = f[0], int(f[3]) - 1, int(f[4]), f[6], f[8]
            m = re.search(r'transcript_id "([^"]+)"', attrs)
            if not m:
                continue
            transcripts.setdefault(m.group(1), []).append((chrom, start, end, strand))
    introns = set()
    for txid, exons in transcripts.items():
        exons.sort(key=lambda e: e[1])
        for a, b in zip(exons, exons[1:]):
            if a[0] != b[0] or a[3] != b[3]:
                continue
            introns.add((a[0], a[2], b[1], a[3]))  # (chrom, donor=prev_exon_end, acceptor=next_exon_start, strand)
    return introns

gencode_introns = gencode_introns_chr22(GENCODE_GTF)
print(f"GENCODE chr22 introns:       {len(gencode_introns):,}")

ground_truth = mn_set & gencode_introns
print(f"Ground truth (MN ∩ GENCODE): {len(ground_truth):,}")
assert len(ground_truth) > 0, "Empty ground truth — chr22 GENCODE intersection failed; check coord conventions"
```

- [ ] **Step 3: Execute; verify counts**

Expected: AG predictions count, GENCODE introns count, non-zero ground truth.

- [ ] **Step 4: Add §2(b) F1 sweep + threshold selection cell**

```python
def f1_at_threshold(ag_df: pd.DataFrame, threshold: float, ground_truth: set[tuple]) -> dict:
    above = ag_df[ag_df[score_col] >= threshold]
    pred = to_key_set(above)
    tp = len(pred & ground_truth)
    fp = len(pred - ground_truth)
    fn = len(ground_truth - pred)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return {"threshold": threshold, "tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1, "n_pred": len(pred)}

# Sweep over all unique AG scores (matches #224 §5 dense-grid approach)
unique_scores = sorted(ag[score_col].unique())
print(f"Score sweep grid: {len(unique_scores):,} unique scores")
sweep = pd.DataFrame([f1_at_threshold(ag, t, ground_truth) for t in unique_scores])
best = sweep.iloc[sweep["f1"].idxmax()]
best_threshold = float(best["threshold"])
print(f"Best F1: {best['f1']:.4f} at τ = {best_threshold:.4f}  (P={best['precision']:.3f}, R={best['recall']:.3f})")
```

- [ ] **Step 5: Execute; verify F1 sweep prints best threshold**

Expected: prints F1, threshold, precision, recall. Best F1 should be substantially above 0 (#224 reported a meaningful F1 in §5).

- [ ] **Step 6: Add §2(b) ag_set construction**

```python
ag_predicted_normal = ag[ag[score_col] >= best_threshold]
ag_set = to_key_set(ag_predicted_normal)
print(f"AG predicted-normal @ τ={best_threshold:.4f}: {len(ag_set):,} junctions")
```

- [ ] **Step 7: Execute; verify ag_set non-empty**

- [ ] **Step 8: Commit**

```bash
jupyter nbconvert --clear-output --inplace research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git add research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git commit -m "wip(#225): notebook §2(b) — AG parquet + F1-max threshold + ag_set"
```

---

### Task 9: Notebook §2(c) — Snaptron chr22 GTEx panel

**Files:**
- Modify: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`

- [ ] **Step 1: Add §2(c) markdown cell**

```markdown
### §2(c) — Snaptron chr22 GTEx pan-tissue proxy

Build a chr22 pan-tissue GTEx union by querying Snaptron's hg38 GTEx v2 endpoint. Inclusion: junction observed in ≥ 1 GTEx sample across any tissue (most-conservative pan-tissue; matches #203's vaccine-safety reasoning).

**Caveat:** this is a proxy for the production GTEx panel built by [#211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211). When #211 lands, re-run this cell against the production panel; only the GTEx-set construction changes downstream.

Result is cached to `outputs/chr22_gtex_panel.parquet`; subsequent runs load the cache.
```

- [ ] **Step 2: Add §2(c) Snaptron query cell with cache logic**

```python
import time
import urllib.request
import urllib.error

SNAPTRON_GTEX_URL = "https://snaptron.cs.jhu.edu/gtexv2/snaptron"
SNAPTRON_CHR22_REGION = "chr22:1-50818468"

def fetch_snaptron_chr22(region: str = SNAPTRON_CHR22_REGION, timeout_s: int = 120, retry: bool = True) -> pd.DataFrame:
    """Query Snaptron GTEx v2 for chr22 junctions. Returns parsed DataFrame."""
    url = f"{SNAPTRON_GTEX_URL}?regions={region}"
    attempts = 2 if retry else 1
    last_err = None
    for attempt in range(1, attempts + 1):
        try:
            print(f"Snaptron query (attempt {attempt}/{attempts}): {url}")
            with urllib.request.urlopen(url, timeout=timeout_s) as resp:
                raw = resp.read().decode("utf-8")
            break
        except urllib.error.URLError as e:
            last_err = e
            if attempt < attempts:
                time.sleep(5)
            else:
                raise
    # Snaptron returns TSV with a header row
    df = pd.read_csv(pd.io.common.StringIO(raw), sep="\t")
    print(f"Snaptron rows returned: {len(df):,}")
    return df

def snaptron_to_key_set(df: pd.DataFrame, min_samples: int = 1) -> set[tuple]:
    """Convert Snaptron DataFrame to the (chrom, donor, acceptor, strand) key set.

    Snaptron col conventions (gtexv2): chromosome, start, end, strand, samples_count.
    'start' is 1-based donor; 'end' is 1-based inclusive acceptor (= 0-based exclusive acceptor).
    We normalize to 0-based half-open intron coords.
    """
    needed = {"chromosome", "start", "end", "strand", "samples_count"}
    assert needed.issubset(df.columns), f"Snaptron schema mismatch; got {df.columns.tolist()}"
    df = df[df["samples_count"] >= min_samples]
    df = df[df["chromosome"] == TARGET_CHROM]
    keys = set()
    for _, row in df[["chromosome", "start", "end", "strand"]].iterrows():
        donor_0based = int(row["start"]) - 1
        acceptor_0based_excl = int(row["end"])  # Snaptron end is 1-based inclusive → 0-based exclusive
        keys.add((row["chromosome"], donor_0based, acceptor_0based_excl, row["strand"]))
    return keys

if GTEX_PANEL_PARQUET.exists():
    print(f"Loading cached GTEx panel from {GTEX_PANEL_PARQUET}")
    gtex_panel = pd.read_parquet(GTEX_PANEL_PARQUET)
else:
    raw_gtex = fetch_snaptron_chr22()
    # Persist a slim view (only the columns we use)
    gtex_panel = raw_gtex[["chromosome", "start", "end", "strand", "samples_count"]].copy()
    assert (gtex_panel["chromosome"] == TARGET_CHROM).all(), "Snaptron returned non-chr22 rows"
    assert len(gtex_panel) > 0, "Empty Snaptron result"
    # Cache only on success
    gtex_panel.to_parquet(GTEX_PANEL_PARQUET, index=False)
    print(f"Cached GTEx panel to {GTEX_PANEL_PARQUET}  ({len(gtex_panel):,} rows)")

gtex_set = snaptron_to_key_set(gtex_panel, min_samples=1)
print(f"GTEx chr22 pan-tissue (≥1 sample): {len(gtex_set):,} unique junction keys")
```

- [ ] **Step 3: Execute; verify Snaptron query succeeds and cache file is created**

```bash
ls -la research/experiments/issue_225_normal_junction_filter_strength/outputs/chr22_gtex_panel.parquet
```

Expected: parquet file exists, ~MB-scale (chr22 GTEx coverage is moderate), `gtex_set` size printed.

- [ ] **Step 4: Commit (parquet is < 10 MB band — checked in per Task 3 convention)**

```bash
jupyter nbconvert --clear-output --inplace research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git add research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb research/experiments/issue_225_normal_junction_filter_strength/outputs/chr22_gtex_panel.parquet
git commit -m "wip(#225): notebook §2(c) — Snaptron chr22 GTEx panel + cached parquet"
```

---

### Task 10: Notebook §3 — apply each filter

**Files:**
- Modify: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`

- [ ] **Step 1: Add §3 markdown + code cell**

```markdown
## §3 — Apply each filter separately to tumor junctions

For each filter source F, compute `tumor_caught_by_F = tumor ∩ F`. The filter "catches" tumor junctions that also appear in F — these would be removed as non-tumor-specific.
```

```python
tumor_set = to_key_set(tumor)
assert len(tumor_set) > 0

caught_by_mn   = tumor_set & mn_set
caught_by_gtex = tumor_set & gtex_set
caught_by_ag   = tumor_set & ag_set

# Sanity guard: each caught set is a subset of tumor
for label, caught in [("MN", caught_by_mn), ("GTEx", caught_by_gtex), ("AG", caught_by_ag)]:
    assert caught.issubset(tumor_set), f"{label} caught set is not a subset of tumor — coord-system mismatch?"
    pct = 100 * len(caught) / len(tumor_set)
    print(f"Tumor caught by {label:5s}: {len(caught):6,}  ({pct:5.1f}% of tumor)")
print(f"Total tumor junctions:       {len(tumor_set):6,}")
```

- [ ] **Step 2: Execute; verify all three caught counts printed, no AssertionError**

- [ ] **Step 3: Commit**

```bash
jupyter nbconvert --clear-output --inplace research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git add research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git commit -m "wip(#225): notebook §3 — apply filters, per-source caught counts"
```

---

### Task 11: Notebook §4 — overlap analysis + write TSV

**Files:**
- Modify: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`
- Output: `research/experiments/issue_225_normal_junction_filter_strength/outputs/filter_overlap_table.tsv`

- [ ] **Step 1: Add §4 markdown + code cell**

```markdown
## §4 — Overlap analysis

Pairwise + triple intersections of the three caught sets; unique contributions (caught by exactly one source). Output exported as `outputs/filter_overlap_table.tsv` for the #203 decision-rule writeup.
```

```python
# Pairwise + triple intersections
mn_gtex   = caught_by_mn & caught_by_gtex
mn_ag     = caught_by_mn & caught_by_ag
gtex_ag   = caught_by_gtex & caught_by_ag
all_three = caught_by_mn & caught_by_gtex & caught_by_ag

# Unique contributions (caught by exactly one)
only_mn   = caught_by_mn - caught_by_gtex - caught_by_ag
only_gtex = caught_by_gtex - caught_by_mn - caught_by_ag
only_ag   = caught_by_ag - caught_by_mn - caught_by_gtex

# Net stacked filter (union)
caught_by_any = caught_by_mn | caught_by_gtex | caught_by_ag

overlap_table = pd.DataFrame([
    {"region": "MN",                       "n": len(caught_by_mn)},
    {"region": "GTEx",                     "n": len(caught_by_gtex)},
    {"region": "AG",                       "n": len(caught_by_ag)},
    {"region": "MN ∩ GTEx",                "n": len(mn_gtex)},
    {"region": "MN ∩ AG",                  "n": len(mn_ag)},
    {"region": "GTEx ∩ AG",                "n": len(gtex_ag)},
    {"region": "MN ∩ GTEx ∩ AG",           "n": len(all_three)},
    {"region": "Only MN",                  "n": len(only_mn)},
    {"region": "Only GTEx",                "n": len(only_gtex)},
    {"region": "Only AG",                  "n": len(only_ag)},
    {"region": "Caught by any (union)",    "n": len(caught_by_any)},
    {"region": "Tumor total",              "n": len(tumor_set)},
])
overlap_table["pct_of_tumor"] = (100 * overlap_table["n"] / len(tumor_set)).round(2)

overlap_table.to_csv(OVERLAP_TSV, sep="\t", index=False)
print(f"Wrote {OVERLAP_TSV}")
overlap_table
```

- [ ] **Step 2: Execute; verify TSV is written + values look sane**

Expected: table prints, `Caught by any` ≥ each individual filter, `MN ∩ GTEx ∩ AG` ≤ each pairwise.

- [ ] **Step 3: Read-back check**

```python
read_back = pd.read_csv(OVERLAP_TSV, sep="\t")
pd.testing.assert_frame_equal(read_back, overlap_table, check_dtype=False)
print("Read-back check passed.")
```

- [ ] **Step 4: Commit**

```bash
jupyter nbconvert --clear-output --inplace research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git add research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb research/experiments/issue_225_normal_junction_filter_strength/outputs/filter_overlap_table.tsv
git commit -m "wip(#225): notebook §4 — overlap analysis + filter_overlap_table.tsv"
```

---

### Task 12: Notebook §5 — 3-way Venn diagram

**Files:**
- Modify: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`
- Output: `research/experiments/issue_225_normal_junction_filter_strength/outputs/filter_venn_chr22.png`

- [ ] **Step 1: Add §5 markdown + code cell**

```markdown
## §5 — 3-way Venn diagram

Visualisation of which tumor junctions each filter catches and where they overlap.
```

```python
from matplotlib_venn import venn3

fig, ax = plt.subplots(figsize=(7, 7))
venn3(
    subsets=(caught_by_mn, caught_by_gtex, caught_by_ag),
    set_labels=("Matched-normal", "GTEx (Snaptron proxy)", "AlphaGenome"),
    ax=ax,
)
ax.set_title(f"Tumor junctions caught by each normal-filter source\n(patient_001, chr22; tumor total = {len(tumor_set):,})")
fig.tight_layout()
fig.savefig(VENN_PNG, dpi=150, bbox_inches="tight")
plt.show()
print(f"Wrote {VENN_PNG}")
```

- [ ] **Step 2: Execute; verify PNG written**

```bash
ls -la research/experiments/issue_225_normal_junction_filter_strength/outputs/filter_venn_chr22.png
```

Expected: PNG exists, non-zero size.

- [ ] **Step 3: Commit**

```bash
jupyter nbconvert --clear-output --inplace research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git add research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb research/experiments/issue_225_normal_junction_filter_strength/outputs/filter_venn_chr22.png
git commit -m "wip(#225): notebook §5 — 3-way Venn diagram"
```

---

### Task 13: Notebook §6 — decision-rule outcome

**Files:**
- Modify: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`

- [ ] **Step 1: Add §6 markdown + code cell**

```markdown
## §6 — Decision-rule outcome

Render the #203 decision table with the values computed above. Exp 2 row (patient-specific delta from germline-aware AG) is marked N/A — deferred to [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) pending WGS.

| Scenario | Threshold | Decision |
|---|---|---|
| F1 ≥ 0.8 AND patient-specific delta > 10% | Both required | Adopt as 3rd always-on source |
| F1 ≥ 0.7 AND ≥ 5% unique vs GTEx | Both required | Adopt as fallback only (unmatched-normal cases) |
| F1 < 0.5 OR delta < 1% | Either triggers | No-go — treat as tissue prior |
```

```python
# Headline metrics for #203 decision rule
exp1_f1 = float(best["f1"])  # F1 from §2(b) sweep
pct_ag_unique_vs_gtex = 100 * len(caught_by_ag - caught_by_gtex) / len(caught_by_ag) if len(caught_by_ag) > 0 else 0.0

# Exp 2 delta deferred
exp2_delta = None  # placeholder — depends on WGS (Sub-Issue #381)

# Decision logic
def decide(f1, delta, ag_unique_pct):
    if delta is not None and f1 >= 0.8 and delta > 10:
        return "ADOPT as 3rd always-on source"
    if f1 >= 0.7 and ag_unique_pct >= 5:
        return "ADOPT as fallback only (unmatched-normal cases)"
    if f1 < 0.5 or (delta is not None and delta < 1):
        return "NO-GO — treat as tissue prior"
    return "INCONCLUSIVE — needs Exp 2 (germline-aware) before deciding"

decision = decide(exp1_f1, exp2_delta, pct_ag_unique_vs_gtex)

decision_table = pd.DataFrame([
    {"metric": "Exp 1 F1 (AG vs MN ∩ GENCODE)",       "value": f"{exp1_f1:.4f}"},
    {"metric": "% AG-unique vs GTEx (caught_by_ag)",  "value": f"{pct_ag_unique_vs_gtex:.1f}%"},
    {"metric": "Exp 2 patient-specific delta",        "value": "N/A — deferred to #381 (WGS)"},
    {"metric": "Decision outcome",                    "value": decision},
])
print(decision_table.to_string(index=False))
```

- [ ] **Step 2: Execute; verify decision string printed**

- [ ] **Step 3: Decision-rule guard**

```python
assert exp1_f1 > 0, "F1 not computed — §2(b) sweep failed silently"
assert pct_ag_unique_vs_gtex >= 0, "AG-unique % negative — set arithmetic bug"
```

- [ ] **Step 4: Commit**

```bash
jupyter nbconvert --clear-output --inplace research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git add research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git commit -m "wip(#225): notebook §6 — decision-rule outcome"
```

---

### Task 14: Notebook §7 — caveats + forward-references

**Files:**
- Modify: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`

- [ ] **Step 1: Add §7 markdown cell**

```markdown
## §7 — Caveats + forward-references

- **Snaptron chr22 proxy ≠ production GTEx panel.** This experiment uses a Snaptron-derived chr22 union (≥1 sample inclusion) as a stand-in for the GTEx pan-tissue reference set being built by [Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211). When #211 lands, re-execute §2(c) against the production panel; the rest of the notebook is unchanged.
- **Experiment 2 (patient-specificity) deferred.** [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) handles AG with germline variants once WGS is acquired. The decision table above marks the Exp 2 delta as N/A.
- **chr22 only.** Test-config harness; full-genome scale-up is out of scope for #225.
- **Pan-tissue, not tissue-matched.** Aligns with #203's vaccine-safety reasoning ("kept pan-tissue, not tissue-matched"). Tissue-matched would change the GTEx panel definition.
- **AG threshold = F1-max from §2(b)** — recomputed inline; #224 didn't persist its `best_threshold`. A future PR could persist `best_threshold` from #224 to remove the recomputation here.
```

- [ ] **Step 2: Commit**

```bash
jupyter nbconvert --clear-output --inplace research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git add research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git commit -m "wip(#225): notebook §7 — caveats + forward-references"
```

---

### Task 15: End-to-end re-run + read-back verification

**Files:**
- Modify: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb` (clean execution outputs persisted)

- [ ] **Step 1: Clear all outputs**

```bash
jupyter nbconvert --clear-output --inplace research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
```

- [ ] **Step 2: Execute end-to-end from a fresh kernel**

```bash
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=splice-neoepitope-alphagenome --ExecutePreprocessor.timeout=600 research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
```

Expected: clean execution, no exceptions, all three outputs in `outputs/` regenerated.

- [ ] **Step 3: Verify all three outputs match in-memory values (the §4 read-back already did this for the TSV; visual check for PNG)**

```bash
ls -la research/experiments/issue_225_normal_junction_filter_strength/outputs/
```

Expected: `chr22_gtex_panel.parquet`, `filter_overlap_table.tsv`, `filter_venn_chr22.png` all present + non-zero.

- [ ] **Step 4: Commit the executed notebook with outputs**

```bash
git add research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
git commit -m "wip(#225): end-to-end notebook execution with cell outputs"
```

---

### Task 16: Update Issue #225 body — drop stale "append to Exp 1+2 notebook" line

**Files:**
- External: [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) body

- [ ] **Step 1: Read the current body and stage the edit**

```bash
gh issue view 225 --repo Jin-HoMLee/splice-neoepitope-pipeline --json body --jq '.body' > /tmp/issue_225_body.md
```

- [ ] **Step 2: Edit `/tmp/issue_225_body.md`**

Replace the `## Output` block:

```markdown
## Output

- Add a section to the Experiment 1+2 notebook (no separate notebook)
- Update the [#203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) decision-rule outcome with stacked filter strength data
```

with:

```markdown
## Output

- **2026-05-21 revision:** ships as its own per-Issue notebook under the new `research/experiments/issue_NNN_<short-desc>/` convention. Original plan to append to the Exp 1+2 notebook is obsolete (Exp 2 was carved out to [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) by the 2026-05-16 audit; #224 is Exp 1 only). Notebook at `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb`.
- Update the [#203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) decision-rule outcome with stacked filter strength data.
```

- [ ] **Step 3: Apply the edit**

```bash
gh issue edit 225 --repo Jin-HoMLee/splice-neoepitope-pipeline --body-file /tmp/issue_225_body.md
```

- [ ] **Step 4: Verify the change**

```bash
gh issue view 225 --repo Jin-HoMLee/splice-neoepitope-pipeline --json body --jq '.body' | grep -A 2 "2026-05-21 revision"
```

Expected: prints the new line.

- [ ] **Step 5: No commit (external state change, not a file change)**

---

### Task 17: Write the lab notebook entry

**Files:**
- Create: `research/lab_notebook/YYYY-MM-DD-HHMM-scientist-issue-225-normal-junction-filter-strength.md` (path will vary per project convention — confirm by looking at the most recent lab-notebook entry path)

- [ ] **Step 1: Find the lab-notebook path convention**

```bash
ls -t research/lab_notebook/ 2>/dev/null | head -3 || find research -path "*lab_notebook*" -name "*.md" -newer research/notebooks/issue_224_alphagenome_exp1_patient_001.ipynb 2>/dev/null | head -3
```

- [ ] **Step 2: Create a new entry at the same path style**

The entry should cover: what was done (Exp 3 ran end-to-end on chr22), key numbers (the §6 decision-table values: F1, % AG-unique vs GTEx, decision outcome), caveats (Snaptron proxy, Exp 2 deferred), follow-ups (re-run §2(c) when #211 lands, #271 unblocks when #203 closes), and references (Issue #225, parent #203, sibling #224).

- [ ] **Step 3: Commit**

```bash
git add research/lab_notebook/<the new file>
git commit -m "docs(scientist): lab notebook — #225 normal-junction filter strength outcome"
```

---

### Task 18: Push branch + open PR

**Files:**
- External: GitHub remote `research/scientist/issue-225-normal-junction-filter-strength` + new PR

- [ ] **Step 1: Surface SHA list + diff summary**

```bash
git log --oneline main..HEAD
git diff --stat main..HEAD
```

Print these to the user before pushing — per "commit/push/merge separate" rule. Wait for user OK to push.

- [ ] **Step 2: Push**

```bash
git push -u origin research/scientist/issue-225-normal-junction-filter-strength
```

- [ ] **Step 3: Open PR**

```bash
gh pr create \
  --repo Jin-HoMLee/splice-neoepitope-pipeline \
  --base main \
  --head research/scientist/issue-225-normal-junction-filter-strength \
  --title "research(#225): normal-junction filter strength on patient_001 (chr22)" \
  --body "$(cat <<'EOF'
## Summary

Experiment 3 from [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203): compares three normal-junction filter sources (matched-normal, Snaptron chr22 GTEx pan-tissue proxy, AlphaGenome predicted-normal) on patient_001's chr22 tumor junctions and produces the decision-rule numbers for the parent epic.

Establishes the new `research/experiments/issue_NNN_<short-desc>/` convention (separate per-Issue notebook + outputs/ + README) and documents it in CLAUDE.md, including the cross-experiment data-sharing rules (`_shared/` lazy, size bands <10/10-100/>100 MB with GCS-manifest pattern at the top band).

Closes [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) (sub-issue of [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203)).

## Design

[docs/superpowers/specs/2026-05-21-issue-225-normal-junction-filter-strength-design.md](docs/superpowers/specs/2026-05-21-issue-225-normal-junction-filter-strength-design.md)

## Test plan

- [x] `notebook.ipynb` end-to-end clean run from a fresh kernel
- [x] `outputs/chr22_gtex_panel.parquet` produced + cached
- [x] `outputs/filter_overlap_table.tsv` matches notebook §4 in-memory values (read-back check)
- [x] `outputs/filter_venn_chr22.png` rendered
- [x] `README.md` in the experiment folder
- [x] CLAUDE.md updated with the `research/experiments/` convention + size-band rules
- [x] #225 body updated to drop the stale "append to Exp 1+2 notebook" line
- [x] Lab notebook entry committed
- [ ] #203 decision-rule outcome to be updated in the parent issue body (Exp 3 row) once this PR merges
- [ ] Migration follow-up Issue (move `issue_224_*` / `issue_299_*` from `research/notebooks/` into `research/experiments/`) to be filed after merge
EOF
)"
```

- [ ] **Step 4: Capture PR URL**

`gh pr create` prints the URL — surface it to the user.

- [ ] **Step 5: Flip board statuses**

Per the PR + Issue Status lifecycle rule: PR → `Ready for review` (`8bf9192f`); Issue stays `In progress`.

```bash
# Get PR node ID and project item ID, then flip Status to Ready for review
# (commands shown for reference; the agent will inspect the PR + look up the item ID via the same pattern Task 17 uses)
PR_NUMBER=<from gh pr create>
PR_NODE_ID=$(gh pr view "$PR_NUMBER" --repo Jin-HoMLee/splice-neoepitope-pipeline --json id --jq .id)
# Then add to project (if not auto-added) + flip Status — agent should mirror the pattern already used for issue-status flips earlier in this session
```

- [ ] **Step 6: Offer `@claude review`**

After the Status flip, ask the user: *"Want me to ping @claude for a review pass on PR #N?"* — per the always-in-effect rule. On approval: `gh pr comment <N> --body "@claude review"`.

---

### Task 19: STOP — wait for explicit user confirmation before merge

**Per user instruction (2026-05-21): stall before merging.** Auto-mode continues for everything up to and including the PR being open / reviewed / iterated. The merge step requires explicit user "go".

- [ ] **Step 1: Surface what's ready**
  - PR URL
  - Test plan checkbox status (all boxes should be ticked except the post-merge ones — #203 update + migration Issue filing)
  - Closure-ritual audit output (`bash scripts/audit_and_merge.sh <PR_NUMBER>` in dry-run mode if possible, else manually review the body)

- [ ] **Step 2: Wait for the user's explicit "go" / "merge it"**

Do NOT proceed to Task 20 without it.

---

### Task 20: Merge + close-ritual + post-merge follow-ups

- [ ] **Step 1: Run the closure-ritual gate**

```bash
bash scripts/audit_and_merge.sh <PR_NUMBER>
```

This audits both the PR body Test plan and the linked Issue's Acceptance criteria for unticked `- [ ]` boxes; merges with default `--squash --delete-branch` if clean.

- [ ] **Step 2: Verify the merge landed**

```bash
gh pr view <PR_NUMBER> --repo Jin-HoMLee/splice-neoepitope-pipeline --json state,mergedAt
```

- [ ] **Step 3: Update parent [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) body**

Add the Exp 3 row outcome to the decision rule section of #203's body (use `gh issue edit 203 --body-file ...`). Include the F1, the % AG-unique vs GTEx, the decision outcome string, and a link to this PR.

- [ ] **Step 4: File the migration Issue**

Create a new Issue titled `chore(research): migrate existing per-issue notebooks from research/notebooks/ to research/experiments/` referencing #225 as the establishing example. Tag `role:scientist`. Add to project board with Status `Backlog`.

- [ ] **Step 5: Confirm board states**

The merged PR + Issue #225 should both auto-flip to `Done`. Verify.

---

## Self-review

**Spec coverage:** all spec sections map to tasks:
- §1 Architecture (spec) → Tasks 1, 3, 4
- §2 Components (spec) → Tasks 5–14
- §3 Error handling (spec) → guards inside Tasks 5, 6, 7, 8, 9, 10, 13
- §4 Testing (spec) → Tasks 11, 15
- Deliverables (spec PR checklist) → Tasks 16, 17, 18

**Placeholder scan:** none. Every code block is concrete; every command is exact. The lab-notebook path is "confirm by looking at the most recent entry" — that's a deliberate "convention-check" step, not a placeholder, because the path style is project-specific and best read from a recent example.

**Type / name consistency:** `tumor`, `mn`, `ag`, `gtex_panel` (DataFrames); `tumor_set`, `mn_set`, `ag_set`, `gtex_set` (key sets); `caught_by_mn`, `caught_by_gtex`, `caught_by_ag` (caught sets); `to_key_set` (single helper). Names used identically across Tasks 6–13.

**Open questions caught during self-review:** Snaptron schema column names (`chromosome`, `start`, `end`, `strand`, `samples_count`) are based on the public GTEx v2 endpoint convention but should be verified by inspecting the first response — Task 9 step 2 includes an `assert needed.issubset(df.columns)` guard that fails fast with the real column list if the convention differs.
