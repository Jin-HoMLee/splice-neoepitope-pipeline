# Issue #225 Slide Deck — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship a lab-seminar-quality Quarto slide deck for [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) in [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452), co-located with the existing notebook, and document the slides-co-location convention in CLAUDE.md.

**Architecture:** Single deck source (`slides.qmd`) + matplotlib regenerator script (`_regenerate_figures.py`) + local `refs.bib` + rendered `slides.html`, all co-located at `research/experiments/issue_225_normal_junction_filter_strength/`. Shared scaffolding (`_template.qmd`, `nature.csl`) stays at `research/slides/`. Quarto renders to reveal.js HTML; PDF disabled per CLAUDE.md footnotehyper conflict.

**Tech Stack:** Quarto (HTML reveal.js), Python 3.13 (matplotlib, pandas, scikit-learn), `gh` CLI, Zotero `pyzotero` (via `zotero_add.py`).

**Spec:** [docs/superpowers/specs/2026-05-22-issue-225-slide-deck-design.md](../specs/2026-05-22-issue-225-slide-deck-design.md)

---

## Working environment

All paths in this plan are relative to repo root: `/Users/jin-holee/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-scientist/`. Branch: `research/scientist/issue-225-normal-junction-filter-strength`. The branch is already checked out from earlier in this session.

**Python interpreter:** `conda activate splice-neoepitope-alphagenome` (same env as the notebook + the prior-art `issue_393` deck regenerator). `which python` should return a path under `~/miniconda3/envs/splice-neoepitope-alphagenome/`. **Do not** use `workflow/tests/.venv/bin/python` (that env doesn't have AG-related deps) or `research/.venv/bin/python` (that env's the notebook-research one and may diverge).

**Quarto:** required for slides render; check with `quarto --version` (any 1.4+ is fine). Install via `brew install --cask quarto` if missing per [research/slides/README.md](research/slides/README.md).

---

## Task 1: Add Snaptron + GTEx citations to Zotero (with dedup check)

**Why first:** the local `refs.bib` (Task 2) hand-copies BibTeX from Zotero. Need the entries to exist there before we copy. Per `feedback_zotero_dedup_check`, always pre-check for existing DOI before calling `zotero_add.py`.

**Files:**
- Read: `research/scripts/zotero_add.py` (existing; one-off script)
- Possibly modify: Zotero collection Z38GTJNW (external state, not a file in this repo)

- [ ] **Step 1: Dedup check for Snaptron DOI**

Run:
```bash
conda activate snakemake
python -c "
from pyzotero import zotero
import os
from dotenv import load_dotenv
load_dotenv('.env')
zot = zotero.Zotero(os.environ['ZOTERO_USER_ID'], 'user', os.environ['ZOTERO_API_KEY'])
items = zot.collection_items_top('Z38GTJNW')
snaptron_doi = '10.1093/bioinformatics/bty025'
gtex_doi = '10.1126/science.aaz1776'
existing = {it['data'].get('DOI', '').lower(): it['data'].get('title', '') for it in items}
print(f'Snaptron ({snaptron_doi}): {\"EXISTS — \" + existing[snaptron_doi.lower()] if snaptron_doi.lower() in existing else \"NOT FOUND\"}')
print(f'GTEx v8 ({gtex_doi}): {\"EXISTS — \" + existing[gtex_doi.lower()] if gtex_doi.lower() in existing else \"NOT FOUND\"}')
"
```

Expected: prints one line per DOI showing either "EXISTS" or "NOT FOUND". Record the results.

- [ ] **Step 2: Add Snaptron if missing**

If Step 1 reported `NOT FOUND` for Snaptron:
```bash
python research/scripts/zotero_add.py 10.1093/bioinformatics/bty025 \
  --tags splicing-database gtex tcga snaptron \
  --note "Snaptron — interactive query interface and database for splicing junctions from RNA-seq projects (GTEx, TCGA, SRA). Used as chr22 GTEx pan-tissue proxy in Issue #225 Experiment 3."
```

If it exists, skip — note its Zotero key for Task 2.

- [ ] **Step 3: Add GTEx v8 if missing**

If Step 1 reported `NOT FOUND` for GTEx v8:
```bash
python research/scripts/zotero_add.py 10.1126/science.aaz1776 \
  --tags gtex transcriptome population-panel \
  --note "GTEx Consortium v8 release — provides the population-scale tissue panel that the Snaptron chr22 proxy stands in for. Production GTEx panel for the pipeline tracked in Issue #211."
```

If it exists, skip — note its Zotero key for Task 2.

- [ ] **Step 4: Capture the Zotero keys for both entries**

Re-run Step 1's check; record the Zotero `key` field for both entries. These will become the BibTeX keys' `note` field in Task 2 (matches prior-art `refs.bib` pattern of `note = {Zotero key XXX}`).

- [ ] **Step 5: Commit (skip if no Zotero state changed)**

No files changed in the repo on this task. Skip git commit; the Zotero side is external state.

---

## Task 2: Create local `refs.bib`

**Files:**
- Create: `research/experiments/issue_225_normal_junction_filter_strength/refs.bib`
- Read: `research/slides/issue_393_alphagenome_chr22_poc/refs.bib` (for reuse pattern)

- [ ] **Step 1: Create the `refs.bib` file**

Create `research/experiments/issue_225_normal_junction_filter_strength/refs.bib`:

```bibtex
% Bibliography for the Issue #225 normal-junction filter strength deck.
% Hand-rolled (matches research/slides/README.md pilot pattern). Keys mirror
% Zotero collection Z38GTJNW where the item exists.
%
% Reuse pattern from research/slides/issue_393_alphagenome_chr22_poc/refs.bib:
%   cotto2023regtools, kim2019hisat2, mudge2025gencode, pedregosa2011scikit,
%   avsec2026alphagenome

@article{cotto2023regtools,
  title     = {Integrated analysis of genomic and transcriptomic data for the discovery of splice-associated variants in cancer},
  author    = {Cotto, Kelsy C. and Feng, Yang-Yang and Ramu, Avinash and Richters, Megan and Freshour, Sharon L. and Skidmore, Zachary L. and Xia, Huiming and McMichael, Joshua F. and Kunisaki, Jason and Campbell, Katie M. and others},
  journal   = {Nature Communications},
  volume    = {14},
  number    = {1},
  pages     = {1589},
  year      = {2023},
  doi       = {10.1038/s41467-023-37266-6},
}

@article{kim2019hisat2,
  title     = {Graph-based genome alignment and genotyping with {HISAT2} and {HISAT-genotype}},
  author    = {Kim, Daehwan and Paggi, Joseph M. and Park, Chanhee and Bennett, Christopher and Salzberg, Steven L.},
  journal   = {Nature Biotechnology},
  volume    = {37},
  number    = {8},
  pages     = {907--915},
  year      = {2019},
  doi       = {10.1038/s41587-019-0201-4},
}

@article{mudge2025gencode,
  title     = {{GENCODE} 2025: reference gene annotation for human and mouse},
  author    = {Mudge, Jonathan M. and Carbonell-Sala, S{\'\i}lvia and Diekhans, Mark and Gonzalez Martinez, Jose and Hunt, Toby and Jungreis, Irwin and Loveland, Jane E. and others},
  journal   = {Nucleic Acids Research},
  volume    = {53},
  number    = {D1},
  pages     = {D966--D975},
  year      = {2025},
  doi       = {10.1093/nar/gkae1078},
  note      = {Describes GENCODE human release v47},
}

@article{pedregosa2011scikit,
  title     = {Scikit-learn: Machine Learning in {Python}},
  author    = {Pedregosa, Fabian and Varoquaux, Ga{\"e}l and Gramfort, Alexandre and Michel, Vincent and Thirion, Bertrand and Grisel, Olivier and Blondel, Mathieu and Prettenhofer, Peter and Weiss, Ron and Dubourg, Vincent and others},
  journal   = {Journal of Machine Learning Research},
  volume    = {12},
  pages     = {2825--2830},
  year      = {2011},
}

@article{avsec2026alphagenome,
  title     = {Advancing regulatory variant effect prediction with {AlphaGenome}},
  author    = {Avsec, {\v Z}iga and others},
  journal   = {Nature},
  year      = {2026},
  doi       = {10.1038/s41586-025-10014-0},
}

@article{wilks2018snaptron,
  title     = {{Snaptron}: querying splicing patterns across tens of thousands of {RNA-seq} samples},
  author    = {Wilks, Christopher and Gaddipati, Phani and Nellore, Abhinav and Langmead, Ben},
  journal   = {Bioinformatics},
  volume    = {34},
  number    = {1},
  pages     = {114--116},
  year      = {2018},
  doi       = {10.1093/bioinformatics/btx547},
  note      = {Zotero key 5MIKCMRQ},
}

@article{gtex2020v8,
  title     = {The {GTEx} Consortium atlas of genetic regulatory effects across human tissues},
  author    = {{GTEx Consortium}},
  journal   = {Science},
  volume    = {369},
  number    = {6509},
  pages     = {1318--1330},
  year      = {2020},
  doi       = {10.1126/science.aaz1776},
  note      = {Zotero key HGHF7Q8R},
}
```

- [ ] **Step 2: Verify the file is well-formed BibTeX**

Run:
```bash
python -c "
import re
with open('research/experiments/issue_225_normal_junction_filter_strength/refs.bib') as f:
    content = f.read()
entries = re.findall(r'@article\{(\w+),', content)
print(f'Entries found: {len(entries)}')
for e in entries:
    print(f'  - {e}')
assert len(entries) == 7, f'Expected 7 entries, got {len(entries)}'
"
```

Expected: `Entries found: 7` and the 7 keys listed (`cotto2023regtools`, `kim2019hisat2`, `mudge2025gencode`, `pedregosa2011scikit`, `avsec2026alphagenome`, `wilks2018snaptron`, `gtex2020v8`).

- [ ] **Step 3: Commit**

```bash
git add research/experiments/issue_225_normal_junction_filter_strength/refs.bib
git commit -m "$(cat <<'EOF'
docs(#225): slide-deck refs.bib — 7 hand-rolled entries

Reuses 5 keys from research/slides/issue_393_alphagenome_chr22_poc/refs.bib
(regtools, HISAT2, GENCODE, sklearn, AlphaGenome). Adds 2 new for #225:
Wilks-Snaptron (the GTEx panel construction source) and GTEx v8 (the
population-panel context the Snaptron chr22 proxy stands in for).

Refs #225

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Create figure regenerator `_regenerate_figures.py`

**Why no formal TDD:** matplotlib figure scripts in this codebase don't use unit tests (see [research/slides/issue_393_alphagenome_chr22_poc/figures/_regenerate_figures.py](../../research/slides/issue_393_alphagenome_chr22_poc/figures/_regenerate_figures.py)). The "test" is visual eyeball + non-zero file size. The figure IS the artifact under review.

**Files:**
- Create: `research/experiments/issue_225_normal_junction_filter_strength/figures/_regenerate_figures.py`
- Read: `research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb` (cells §2(b) and §3-§4 for source logic)
- Read: `research/slides/issue_393_alphagenome_chr22_poc/figures/_regenerate_figures.py` (style + structure reference)

- [ ] **Step 1: Create the regenerator script**

Create `research/experiments/issue_225_normal_junction_filter_strength/figures/_regenerate_figures.py`:

```python
"""Regenerate the 2 deck-only figures for Issue #225.

Run from repo root:
    conda activate splice-neoepitope-alphagenome
    python research/experiments/issue_225_normal_junction_filter_strength/figures/_regenerate_figures.py

Outputs (next to this script):
    pr_curve.png      — Precision-Recall curve over universe-restricted F1 sweep
    caught_bar.png    — Tumor caught (% of tumor total) per filter source + union

Logic mirrors research/experiments/issue_225_normal_junction_filter_strength/notebook.ipynb
sections 2(b) and 3-4 verbatim — same data sources, same universe construction.

Re-run after editing the notebook to keep slide figures in sync. Both the notebook
and this script remain reproducible from the same cached parquet.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve

REPO_ROOT = Path(__file__).resolve().parents[4]
EXPERIMENT_DIR = REPO_ROOT / "research" / "experiments" / "issue_225_normal_junction_filter_strength"
OUTPUTS_DIR = EXPERIMENT_DIR / "outputs"
FIGURES_DIR = Path(__file__).resolve().parent
TARGET_CHROM = "chr22"

MATCHED_NORMAL_TSV = REPO_ROOT / "results" / "patient_001_test" / "alignment" / "SRR9143065_test" / "junctions.tsv"
AG_PARQUET = REPO_ROOT / "research" / "notebooks" / "issue_224_alphagenome_exp1_outputs" / "chr22_stomach_predicted_junctions.parquet"
GENCODE_GTF = REPO_ROOT / "resources" / "test" / "chr22.gtf.gz"


def load_pipeline_junctions(tsv_path: Path) -> pd.DataFrame:
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


def gencode_introns_chr22(gtf_path: Path) -> set[tuple]:
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
    introns: set[tuple] = set()
    for exons in transcripts.values():
        exons.sort(key=lambda e: e[1])
        for a, b in zip(exons, exons[1:]):
            if a[0] != b[0] or a[3] != b[3]:
                continue
            donor, acceptor = a[2], b[1]
            if donor >= acceptor:
                continue
            introns.add((a[0], donor, acceptor, a[3]))
    return introns


def to_key_set(df: pd.DataFrame) -> set[tuple]:
    return set(map(tuple, df[["chrom", "donor", "acceptor", "strand"]].itertuples(index=False, name=None)))


def make_pr_curve() -> None:
    """Plot PR curve from universe-restricted F1 sweep (matches notebook §2(b))."""
    mn = load_pipeline_junctions(MATCHED_NORMAL_TSV)
    mn_set = to_key_set(mn)
    gencode_introns = gencode_introns_chr22(GENCODE_GTF)
    ag = pd.read_parquet(AG_PARQUET)
    score_col = "score" if "score" in ag.columns else next(
        (c for c in ag.columns if "score" in c.lower() or "prob" in c.lower()), None
    )
    assert score_col is not None, f"No score column in AG parquet; cols: {ag.columns.tolist()}"
    ag_dedup = ag.groupby(["chrom", "donor", "acceptor", "strand"], as_index=False)[score_col].max()
    ag_lookup = dict(zip(
        map(tuple, ag_dedup[["chrom", "donor", "acceptor", "strand"]].itertuples(index=False, name=None)),
        ag_dedup[score_col].to_numpy(),
    ))
    universe_keys = list(gencode_introns)
    scores = np.array([ag_lookup.get(k, 0.0) for k in universe_keys], dtype=float)
    labels = np.array([1 if k in mn_set else 0 for k in universe_keys], dtype=np.int8)

    p, r, t = precision_recall_curve(labels, scores)
    denom = p + r
    f1 = np.where(denom > 0, 2 * p * r / np.where(denom > 0, denom, 1.0), 0.0)
    best_idx = int(np.argmax(f1[:-1]))
    ap = float(average_precision_score(labels, scores))
    baseline = labels.sum() / len(labels)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(r, p, lw=2, color="#1e5ba8", label=f"AG vs MN ∩ GENCODE  (AP = {ap:.3f})")
    ax.axhline(baseline, ls="--", color="grey", alpha=0.6, label=f"Baseline (prevalence = {baseline:.3f})")
    ax.plot(r[best_idx], p[best_idx], "o", color="crimson", markersize=10,
            label=f"F1-max = {f1[best_idx]:.3f} at τ = {t[best_idx]:.3f}\n(P = {p[best_idx]:.3f}, R = {r[best_idx]:.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("AlphaGenome PR curve — chr22 universe-restricted F1 sweep\n(universe = GENCODE chr22 introns; positives = MN ∩ GENCODE)")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=10)
    fig.tight_layout()
    out = FIGURES_DIR / "pr_curve.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}  (F1={f1[best_idx]:.4f}, τ={t[best_idx]:.4f}, AP={ap:.4f})")


def make_caught_bar() -> None:
    """Bar chart of tumor caught per filter source + union (matches notebook §3-§4)."""
    df = pd.read_csv(OUTPUTS_DIR / "filter_overlap_table.tsv", sep="\t")
    labels = ["MN", "GTEx", "AG", "Caught by any (union)"]
    rows = df[df["region"].isin(labels)].set_index("region").loc[labels].reset_index()

    colors = ["#1e5ba8", "#2e8b57", "#b8860b", "#666666"]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    bars = ax.bar(rows["region"], rows["pct_of_tumor"], color=colors)
    for bar, n, pct in zip(bars, rows["n"], rows["pct_of_tumor"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.4,
                f"{int(n):,}\n({pct:.1f}%)", ha="center", va="bottom", fontsize=11)
    ax.set_ylabel("% of tumor junctions caught")
    ax.set_title("Tumor junctions caught by each normal-filter source\n(patient_001 chr22; tumor n = 1,872)")
    ax.set_ylim(0, max(rows["pct_of_tumor"]) * 1.18)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    out = FIGURES_DIR / "caught_bar.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


def main() -> None:
    print(f"Repo root:      {REPO_ROOT}")
    print(f"Experiment dir: {EXPERIMENT_DIR.relative_to(REPO_ROOT)}")
    print(f"Outputs dir:    {FIGURES_DIR.relative_to(REPO_ROOT)}")
    make_pr_curve()
    make_caught_bar()
    print("Done.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke-test the imports**

Run:
```bash
conda activate splice-neoepitope-alphagenome
python -c "
import sys
sys.path.insert(0, 'research/experiments/issue_225_normal_junction_filter_strength/figures')
import _regenerate_figures
print('Module loaded; functions:', [f for f in dir(_regenerate_figures) if not f.startswith('_')])
"
```

Expected: prints `Module loaded; functions: ['REPO_ROOT', 'EXPERIMENT_DIR', ..., 'make_pr_curve', 'make_caught_bar', 'main', ...]`. If `ImportError`, the env is wrong.

- [ ] **Step 3: Commit (script only — figures come in Task 4)**

```bash
git add research/experiments/issue_225_normal_junction_filter_strength/figures/_regenerate_figures.py
git commit -m "$(cat <<'EOF'
feat(#225): slide-deck figure regenerator (PR curve + caught bar)

Mirrors notebook sections 2(b) and 3-4 verbatim — same data sources,
same universe construction, same matplotlib styling pattern as
research/slides/issue_393_alphagenome_chr22_poc/figures/_regenerate_figures.py.

Produces two PNGs (pr_curve.png, caught_bar.png). Figures committed in
the next commit (after a successful regen run).

Refs #225

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Run regenerator + commit figures

**Files:**
- Create: `research/experiments/issue_225_normal_junction_filter_strength/figures/pr_curve.png`
- Create: `research/experiments/issue_225_normal_junction_filter_strength/figures/caught_bar.png`

- [ ] **Step 1: Run the regenerator**

```bash
conda activate splice-neoepitope-alphagenome
python research/experiments/issue_225_normal_junction_filter_strength/figures/_regenerate_figures.py
```

Expected output:
```
Repo root:      /Users/jin-holee/dev/GitHub/Jin-HoMLee/splice-neoepitope-pipeline-scientist
Experiment dir: research/experiments/issue_225_normal_junction_filter_strength
Outputs dir:    research/experiments/issue_225_normal_junction_filter_strength/figures
Wrote .../pr_curve.png  (F1=0.3000, τ=3.1633, AP=...)
Wrote .../caught_bar.png
Done.
```

The F1 value MUST be `0.3000` and τ MUST be `3.1633` — same as the notebook §2(b) output. If different, the regenerator's data sources don't match the notebook's (likely a path issue).

- [ ] **Step 2: Visual inspection of both PNGs**

Open the two PNGs in Preview:
```bash
open research/experiments/issue_225_normal_junction_filter_strength/figures/pr_curve.png
open research/experiments/issue_225_normal_junction_filter_strength/figures/caught_bar.png
```

Checks:
- `pr_curve.png`: blue curve descending from top-left, dashed grey baseline near y≈0.034, crimson F1-max dot. Title mentions "universe-restricted". Legend readable.
- `caught_bar.png`: 4 bars (MN, GTEx, AG, union) with counts + percentages on top. Heights match: MN 4.9 / GTEx 25.8 / AG 6.6 / union 26.9. Y-axis "% of tumor junctions caught".

If anything looks broken (wrong colors / unreadable text / overlapping labels), edit `_regenerate_figures.py`, re-run, re-eyeball, then proceed.

- [ ] **Step 3: Sanity-check file sizes**

Run:
```bash
ls -la research/experiments/issue_225_normal_junction_filter_strength/figures/*.png
```

Expected: both PNGs are 20–200 KB (matplotlib PNGs at dpi=150 land in this range). If <5 KB, the figure render failed silently.

- [ ] **Step 4: Commit**

```bash
git add research/experiments/issue_225_normal_junction_filter_strength/figures/pr_curve.png \
        research/experiments/issue_225_normal_junction_filter_strength/figures/caught_bar.png
git commit -m "$(cat <<'EOF'
feat(#225): slide-deck figures — PR curve + caught bar

Generated by figures/_regenerate_figures.py. Matches notebook headline
numbers (F1=0.3000 at τ=3.1633; MN 4.9 / GTEx 25.8 / AG 6.6 / union 26.9
% of 1,872 tumor junctions).

Refs #225

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Write `slides.qmd`

**Files:**
- Create: `research/experiments/issue_225_normal_junction_filter_strength/slides.qmd`
- Read: `research/slides/issue_393_alphagenome_chr22_poc/slides.qmd` (structure reference)
- Read: `research/slides/_template.qmd` (template reference)

- [ ] **Step 1: Create `slides.qmd`**

Create `research/experiments/issue_225_normal_junction_filter_strength/slides.qmd`:

```markdown
---
title: "Is AlphaGenome worth adding as a third normal-junction filter?"
subtitle: "Comparative filter strength on patient_001 (chr22) · Issue #225"
author: "Jin-Ho Lee"
institute: "Splice Neoepitope Pipeline · JH M Lee Lab"
date: 2026-05-22
date-format: "YYYY-MM-DD"
format:
  revealjs:
    theme: simple
    incremental: false
    slide-number: c/t
    progress: true
    fig-align: center
    fig-cap-location: bottom
    transition: fade
    background-transition: fade
    width: 1280
    height: 720
# PDF (beamer) render disabled — Quarto/pandoc-emitted preamble hits a
# `\makesavenoteenv{longtable}` interaction with footnotehyper that fails compile.
# For a PDF handout: open slides.html in Chrome and "Print → Save as PDF".
bibliography: refs.bib
csl: ../../slides/nature.csl
execute:
  echo: false
  warning: false
---

## The question

::: {.r-fit-text}
Is AlphaGenome worth
:::

::: {.r-fit-text}
adding as a third
:::

::: {.r-fit-text}
normal-junction filter?
:::

. . .

::: {.callout-note appearance="simple"}
**Why now:** parent [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) — rethink normal-junction filtering with population panel + AlphaGenome as candidate fallback. This deck is the Experiment 3 outcome ([Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225)), feeding back into the #203 decision rule.
:::

---

## Why this matters {.smaller}

**Current pipeline:** matched-normal RNA-seq filters tumor junctions → only tumor-exclusive junctions reach neoepitope prediction.

**Problem:** matched-normal samples are not always available (cost, logistics, retrospective cohorts).

**Three candidate replacements:**

::: {.incremental}
1. **Matched-normal (MN)** — current pipeline source. Patient-specific, but availability-limited.
2. **GTEx pan-tissue** — population-scale union via Snaptron [@wilks2018snaptron; @gtex2020v8] (production panel tracked in [Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211))
3. **AlphaGenome predicted-normal** [@avsec2026alphagenome] — reference-only, tissue-conditioned ML
:::

. . .

**Today's question:** when you already have GTEx, does AG add value on top? Or is it subsumed?

---

## Method — block diagram {.smaller}

```{mermaid}
%%| fig-width: 11
flowchart LR
    T["Tumor RNA-seq<br/>(SRR9143066)"] -->|"HISAT2 + regtools"| TJ["Tumor junctions<br/>n = 1,872"]
    MN["Matched-normal<br/>RNA-seq (SRR9143065)"] -->|"HISAT2 + regtools"| MS["MN junctions<br/>n = 1,714"]
    GT["Snaptron GTEx v2<br/>(hg38)"] -->|"≥1 sample,<br/>chr22"| GS["GTEx panel<br/>n = 880,769"]
    AG["AlphaGenome<br/>chr22 predictions"] -->|"τ = F1-max<br/>vs MN ∩ GENCODE"| AS["AG predicted-normal<br/>n = 448"]
    TJ --> R{{"Set intersections<br/>tumor ∩ each filter"}}
    MS --> R
    GS --> R
    AS --> R
    R --> V[["Venn + decision rule"]]

    style R fill:#ffe6a3,stroke:#b8860b,stroke-width:2px
    style V fill:#a3d8ff,stroke:#1e5ba8,stroke-width:2px
```

**Tooling:** `HISAT2` [@kim2019hisat2] + `regtools` [@cotto2023regtools] (alignment + junction extraction); GENCODE v47 [@mudge2025gencode] for the F1-sweep universe.

---

## Data scope {.smaller}

| Field | Value |
|---|---|
| Patient | patient_001 |
| Samples | Tumor `SRR9143066` + matched normal `SRR9143065` (gastric cancer, adjacent stomach) |
| Read depth | 500K reads / sample (test dataset) |
| Chromosome | **chr22 only** (test-config harness) |
| Reference | GRCh38 (UCSC hg38) |
| Tumor set | 1,872 junctions |
| MN set | 1,714 junctions |
| GTEx panel | 880,769 junctions (Snaptron chr22 pan-tissue, ≥1 sample) |
| AG predicted-normal | 448 junctions @ τ = F1-max |

::: {.callout-warning appearance="simple"}
**Snaptron proxy ≠ production GTEx.** The chr22 union from Snaptron stands in for the production GTEx panel (tracked in [Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211)). Pan-tissue not tissue-matched — matches the vaccine-safety reasoning from [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203).
:::

---

## AG threshold via F1 sweep {.smaller}

**Universe:** GENCODE chr22 introns (n = 7,731). **Positives:** matched-normal ∩ GENCODE (n = 259) — the confirmed tissue-expressed set. **Negatives:** 7,472. **AG-scored:** 5,728 / 74.1% of universe.

::: {.callout-note appearance="simple"}
Universe-restricted F1 (not full tumor set) — matches [Issue #224](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/224) §5 semantics so this F1 is directly comparable to #224's reported value.
:::

![Precision–Recall curve, universe-restricted; F1-max marked.](figures/pr_curve.png){width=72%}

**Result:** F1 = **0.300** at τ = **3.16** (P = 0.238, R = 0.405). Threshold-binarised AG predicted-normal set: **448 junctions**.

Compute via `sklearn.metrics.precision_recall_curve` [@pedregosa2011scikit] over the dense score grid (5,729 unique scores).

---

## Caught counts {.smaller}

For each filter F, `caught_by_F = tumor ∩ F` — i.e. tumor junctions the filter would remove as non-tumor-specific.

![Per-source caught counts. Tumor total = 1,872 (chr22).](figures/caught_bar.png){width=70%}

::: {.incremental}
- **MN catches 4.9%** — sparse because this is a 500K-read patient sample.
- **GTEx catches 25.8%** — population scale advantage shows immediately.
- **AG catches 6.6%** — looks marginal on its own.
- **Stacked union catches 26.9%** — only +1.1pp over GTEx alone. Hint at the killer slide.
:::

---

## The killer slide {.smaller}

![Three-way Venn. AG sits entirely inside GTEx.](outputs/filter_venn_chr22.png){width=68%}

::: {.fragment}
**Two over-determining findings:**

- **Only AG = 0.** Not a single AG-caught tumor junction is unique to AG.
- **GTEx ∩ AG = AG (124 = 124).** AG is a *strict subset* of GTEx at the F1-max threshold.
:::

. . .

**Implication:** stacking AG on top of GTEx adds zero new filtering signal. Even Exp 2 (germline-aware AG, deferred to [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) pending WGS) cannot rescue this — the AG set itself is contained in GTEx.

---

## Decision rule {.smaller}

| Scenario | Threshold | Decision |
|---|---|---|
| F1 ≥ 0.8 **AND** patient-specific delta > 10% | both required | Adopt as 3rd always-on source |
| F1 ≥ 0.7 **AND** ≥ 5% unique vs GTEx | both required | Adopt as fallback only |
| **F1 < 0.5 OR delta < 1%** | **either triggers** | **No-go — treat as tissue prior** |

::: {.incremental}
- F1 = **0.300** → trips the `F1 < 0.5` NO-GO clause directly.
- **0.0%** AG-unique-vs-GTEx → fails the fallback tier's `≥ 5%` test.
- Exp 2 delta deferred — **not load-bearing**. Exp 1's F1 alone is conclusive.
:::

. . .

::: {.callout-important}
The verdict is **over-determined** — F1 and unique-vs-GTEx independently land in NO-GO.
:::

---

## Decision

::: {.r-fit-text}
🔴 NO-GO
:::

::: {.r-fit-text}
treat as tissue prior
:::

::: {.fragment}
- AlphaGenome predicted-normal does **not** add filtering value over GTEx on chr22 patient_001.
- Filter stack stays at **MN + GTEx** (with AG dropped as a candidate third source).
:::

---

## Caveats — read before citing the headline {.smaller}

::: {.incremental}
1. **chr22 PoC scope.** Generalisation to full genome must be re-validated.
2. **Snaptron proxy ≠ production GTEx.** Production panel via [Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211); when it lands, re-run the GTEx-set construction cell. Conclusion is unlikely to flip — AG is subsumed by Snaptron's pan-tissue union; the production panel will be at least as inclusive.
3. **Exp 2 deferred.** [Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) handles AG with germline variants once WGS lands. Not load-bearing for this NO-GO.
4. **Pan-tissue not tissue-matched.** Matches [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) vaccine-safety framing.
5. **AG threshold = F1-max from #224 §2(b)** — recomputed inline because #224 didn't persist `best_threshold`.
:::

---

## Next steps {.smaller}

::: {.incremental}
- **Close [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225)** after PR #452 merges; carry the NO-GO row into parent [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) Exp 3 (immediate post-merge).
- **[Issue #211](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/211) (production GTEx)** — when it lands, re-run §2(c) cell against the production panel; sanity-check the NO-GO holds.
- **[Sub-Issue #381](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/381) (Exp 2: germline-aware AG)** — patient_001 with WGS; nice-to-have follow-up, but the Exp 1 NO-GO already binds.
- **[Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) (notebooks/slides migration)** — moves `issue_224_*` + `issue_393_*` into the `research/experiments/` convention established by this work.
:::

---

## References

::: {#refs}
:::

<!-- bibliography renders here, auto-populated by Pandoc from refs.bib via [@cite-keys] in the text above -->
```

- [ ] **Step 2: Verify file structure**

Run:
```bash
grep -c '^---$' research/experiments/issue_225_normal_junction_filter_strength/slides.qmd
grep -c '^## ' research/experiments/issue_225_normal_junction_filter_strength/slides.qmd
```

Expected: at least one `---` (YAML closer) + ~12 slide separators (`---` between slides — count varies but should be 10+), and 12 H2 headings (one per slide title). If the count is wrong, the slide separators or H2 headings got dropped.

- [ ] **Step 3: Commit**

```bash
git add research/experiments/issue_225_normal_junction_filter_strength/slides.qmd
git commit -m "$(cat <<'EOF'
docs(#225): slide-deck source (slides.qmd) — 12 slides incl. references

Co-located form: lives at research/experiments/issue_225_*/slides.qmd
next to the notebook + outputs (vs the parallel research/slides/issue_NNN/
form used by the issue_393 pilot). Shared scaffolding (_template.qmd,
nature.csl) still referenced from research/slides/ via ../../slides/.

Citations: 7 keys (5 reused from issue_393 + 2 new for #225). Figures
sourced from figures/ (regenerator) + outputs/ (notebook canonical).

Refs #225

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Render `slides.qmd` to `slides.html`

**Files:**
- Create: `research/experiments/issue_225_normal_junction_filter_strength/slides.html`
- Create: `research/experiments/issue_225_normal_junction_filter_strength/slides_files/` (directory of reveal.js assets)

- [ ] **Step 1: Render via Quarto**

```bash
cd research/experiments/issue_225_normal_junction_filter_strength
quarto render slides.qmd --to revealjs
cd -
```

Expected:
- New file `slides.html` (~50–200 KB typical for a 12-slide reveal.js deck)
- New directory `slides_files/` with reveal.js assets

If `quarto render` errors:
- **YAML parse error:** check the `slides.qmd` front matter — common issue is unbalanced quotes
- **Bibliography error:** check `refs.bib` keys match the `[@cite-keys]` in the body
- **CSL not found:** confirm `research/slides/nature.csl` exists (it should — committed earlier)
- **Mermaid render error:** Quarto 1.4+ has built-in mermaid; if a flag is needed, the prior-art deck didn't need it. If error persists, simplify the diagram or move mermaid to a code block.

- [ ] **Step 2: Open `slides.html` and click through all 12 slides**

```bash
open research/experiments/issue_225_normal_junction_filter_strength/slides.html
```

Click-through checks:
- Title slide renders with framing question, r-fit-text scales correctly
- Why-this-matters: incremental bullets reveal correctly
- Method block diagram: mermaid renders (left-to-right flow visible)
- Data scope: table renders with the 9 rows
- F1 sweep: `pr_curve.png` displays
- Caught counts: `caught_bar.png` displays
- Venn slide: `outputs/filter_venn_chr22.png` displays (relative path from `slides.qmd` is `outputs/filter_venn_chr22.png`, which `outputs/` is one level down from `slides.qmd`)
- Decision rule table renders
- Decision slide: 🔴 r-fit-text displays
- Caveats: 5 incremental bullets
- Next steps: 4 incremental bullets
- References: bibliography auto-populated

If any figure shows as a broken-image icon, the relative path is wrong. Fix and re-render.

- [ ] **Step 3: Commit rendered output**

```bash
git add research/experiments/issue_225_normal_junction_filter_strength/slides.html \
        research/experiments/issue_225_normal_junction_filter_strength/slides_files/
git commit -m "$(cat <<'EOF'
build(#225): render slides.html — committed for PR review ergonomics

Divergence from research/slides/README.md gitignore-HTML rule: this deck
is added to PR #452 explicitly to support human review, so reviewers
should be able to open the deck without a Quarto install. Tracked as a
follow-up decision under #455 (either harmonise by committing HTML for
migrated decks too, or revert here once the convention is settled).

Refs #225

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Update CLAUDE.md with slides-co-location convention

**Files:**
- Modify: [CLAUDE.md](../../CLAUDE.md) — "Slide decks for experiment Issues" section

- [ ] **Step 1: Read the current section**

```bash
grep -n "Slide decks for experiment Issues" CLAUDE.md
```

Note the line number. The section runs ~10 lines below. The current text reads (verbatim, as of pre-this-task):

> Every experiment-tier Issue ships a Quarto slide deck alongside the per-patient notebook + manuscript work. Decks live at `research/slides/issue_NNN_<short-content-desc>/slides.qmd`. See [`research/slides/README.md`](research/slides/README.md) for the full convention (Quarto rationale, render commands, figure-source pattern, Zotero linkage, install). Scope: **one deck per experiment Issue**, not per sub-issue; no decks for closure tasks, doc updates, or single-fix PRs. Audience: lab seminar / external talk. Figures regenerate from `research/experiments/issue_NNN_<short>/outputs/*.parquet` so the **notebook stays canonical**.

- [ ] **Step 2: Apply the Edit**

Use the Edit tool with:
- **old_string** (must match verbatim, including the `research/slides/issue_NNN_...` path):
  ```
  Every experiment-tier Issue ships a Quarto slide deck alongside the per-patient notebook + manuscript work. Decks live at `research/slides/issue_NNN_<short-content-desc>/slides.qmd`. See [`research/slides/README.md`](research/slides/README.md) for the full convention (Quarto rationale, render commands, figure-source pattern, Zotero linkage, install). Scope: **one deck per experiment Issue**, not per sub-issue; no decks for closure tasks, doc updates, or single-fix PRs. Audience: lab seminar / external talk. Figures regenerate from `research/experiments/issue_NNN_<short>/outputs/*.parquet` so the **notebook stays canonical**.
  ```
- **new_string**:
  ```
  Every experiment-tier Issue ships a Quarto slide deck alongside the per-patient notebook + manuscript work. Decks live **co-located with the notebook** at `research/experiments/issue_NNN_<short-content-desc>/slides.qmd`. Shared scaffolding (`_template.qmd`, `nature.csl`) stays centralized at `research/slides/` and is referenced via `csl: ../../slides/nature.csl` from each deck. See [`research/slides/README.md`](research/slides/README.md) for the full convention (Quarto rationale, render commands, figure-source pattern, Zotero linkage, install). Scope: **one deck per experiment Issue**, not per sub-issue; no decks for closure tasks, doc updates, or single-fix PRs. Audience: lab seminar / external talk. Figures regenerate from `research/experiments/issue_NNN_<short>/outputs/*.parquet` (notebook outputs) and from a local `figures/_regenerate_figures.py` (deck-only figures) so the **notebook stays canonical**. Rationale for co-location: notebook + outputs + deck rename / archive / migrate as a unit; figure paths shorten from `../../experiments/issue_NNN/outputs/...` to `outputs/...`. (`research/slides/issue_393_alphagenome_chr22_poc/` predates this rule and is migrated under [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455).)
  ```

- [ ] **Step 3: Verify the change applied**

```bash
grep -A 1 "Slide decks for experiment Issues" CLAUDE.md | head -20
```

Expected: the section now starts with "Every experiment-tier Issue ships a Quarto slide deck **co-located with the notebook**..." (with the bold marker).

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md
git commit -m "$(cat <<'EOF'
docs(claude-md): slides-co-location convention — decks live with notebooks

Decks now live at research/experiments/issue_NNN_<short>/slides.qmd
co-located with the notebook + outputs they document. Shared scaffolding
(_template.qmd, nature.csl) stays at research/slides/.

Rationale (per the 2026-05-21 brainstorm captured in lab notebook):
notebook + outputs + deck rename/archive/migrate as a unit; deck figure
paths shorten from ../../experiments/issue_NNN/outputs/... to outputs/...

Pilot deck research/slides/issue_393_* predates this rule; migrated
under #455.

Refs #225

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Update [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) body — descope handled items

**Files:**
- Modify: [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) body (external — via `gh issue edit`)

- [ ] **Step 1: Fetch current Issue #455 body**

```bash
gh issue view 455 --json body --jq .body > /tmp/issue_455_body.md
wc -l /tmp/issue_455_body.md
```

Expected: ~30 lines.

- [ ] **Step 2: Edit `/tmp/issue_455_body.md` to descope handled items**

Open `/tmp/issue_455_body.md` in your editor and make these changes:

1. **Under "## Convention docs":** strike (or remove) the line:
   ```
   - [ ] Update CLAUDE.md *"Experiment notebooks live under `research/experiments/`"* section to document the slides-co-location rule…
   ```
   Replace with:
   ```
   - [x] ~~Update CLAUDE.md to document slides-co-location rule~~ — handled in [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) ([Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225))
   ```

2. **Under "## Out of scope":** strike the line:
   ```
   - A new deck for [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) itself — that is a separate follow-up to file once the slides-co-location convention is documented.
   ```
   Replace with:
   ```
   - ~~A new deck for Issue #225 itself~~ — shipped in [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) alongside the notebook.
   ```

3. **Under "## Acceptance criteria":** the line `- [ ] CLAUDE.md updated with the slides-co-location convention` is now handled. Change to:
   ```
   - [x] CLAUDE.md updated with the slides-co-location convention — handled in [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452)
   ```

- [ ] **Step 3: Push the edit back to GitHub**

```bash
gh issue edit 455 --body-file /tmp/issue_455_body.md
```

Expected: prints `https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455` on success.

- [ ] **Step 4: Verify the changes**

```bash
gh issue view 455 --json body --jq .body | grep -E "PR #452|slides-co-location|new deck for"
```

Expected: shows the 3 lines mentioning PR #452 (one per descoped item).

- [ ] **Step 5: No commit needed**

Issue body lives on GitHub, not in the repo. No file change.

---

## Task 9: Update [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) body — add deck Test plan items

**Files:**
- Modify: PR #452 body (external — via `gh pr edit`)

- [ ] **Step 1: Fetch current PR #452 body**

```bash
gh pr view 452 --json body --jq .body > /tmp/pr_452_body.md
```

- [ ] **Step 2: Append deck-related Test plan items**

Open `/tmp/pr_452_body.md` and find the `## Test plan` section. Append the following items at the end of the existing list (after the last `- [x] Migration follow-up Issue — filed as Issue #455…` line):

```
- [ ] `slides.qmd` renders cleanly to `slides.html` via `quarto render --to revealjs`
- [ ] `figures/_regenerate_figures.py` produces both `pr_curve.png` + `caught_bar.png` from notebook outputs; F1=0.300 at τ=3.16 (matches notebook)
- [ ] All figure refs in `slides.qmd` resolve (Venn from `outputs/`, PR curve + bar from `figures/`)
- [ ] CLAUDE.md slides-co-location rule documented
- [ ] Issue #455 body updated to descope CLAUDE.md convention-docs item + #225 deck item
- [ ] Lab notebook entry for the 2026-05-22 deck-creation session
```

- [ ] **Step 3: Update the PR body**

```bash
gh pr edit 452 --body-file /tmp/pr_452_body.md
```

Expected: prints the PR URL on success.

- [ ] **Step 4: Verify**

```bash
gh pr view 452 --json body --jq .body | grep -c "slides"
```

Expected: 4 or more matches (the new Test plan items mention `slides` multiple times).

- [ ] **Step 5: No commit needed**

PR body lives on GitHub.

---

## Task 10: Lab notebook entry for the deck-creation session

**Files:**
- Modify: [research/lab_notebook/scientist.md](../../research/lab_notebook/scientist.md)

- [ ] **Step 1: Determine entry timestamp**

Get current UTC time:
```bash
date -u "+%Y-%m-%d %H:%M UTC"
```

Use this in the entry header.

- [ ] **Step 2: Read the existing lab notebook to find the right insertion point**

```bash
head -30 research/lab_notebook/scientist.md
```

The file is reverse-chronological (newest at top). Insert after the first `---` separator line (so the new entry becomes the newest entry).

- [ ] **Step 3: Add the new entry**

Use the Edit tool to insert below the first `---` separator. The entry text:

```markdown

### <TIMESTAMP_FROM_STEP_1> — Editor: Scientist

#### [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) (filter strength deck) — slide deck added to PR #452 + slides-co-location convention documented

User asked for the [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) slide deck to ship in [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) (vs deferring to [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455)) on the grounds that notebook + deck "belong together and are necessary for human review". Agreed; established the slides-co-location convention informally on 2026-05-21 19:38 UTC (lab notebook entry above) is now formalised in CLAUDE.md as part of this PR.

**Scope shipped here:**

- `research/experiments/issue_225_*/slides.qmd` (12 slides, lab-seminar-quality, reveal.js HTML render)
- `figures/_regenerate_figures.py` + `figures/pr_curve.png` + `figures/caught_bar.png`
- `outputs/filter_venn_chr22.png` re-used from the notebook (notebook stays canonical per CLAUDE.md)
- `refs.bib` — 7 entries (5 reused from issue_393 deck, 2 new: Wilks-Snaptron + GTEx v8)
- `slides.html` committed for PR-review ergonomics (divergence from `research/slides/README.md`'s gitignore rule — tracked as a follow-up under [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455))
- CLAUDE.md "Slide decks for experiment Issues" section updated with the co-location rule + rationale
- [Issue #455](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/455) body descoped: removed the CLAUDE.md convention-docs item and the "[Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) deck follow-up" Out-of-scope line (both shipped here)
- PR #452 Test plan amended with 6 deck-related items

**Process note.** Followed the brainstorming → writing-plans → executing-plans superpowers flow (vs free-form). Spec at `docs/superpowers/specs/2026-05-22-issue-225-slide-deck-design.md`; plan at `docs/superpowers/plans/2026-05-22-issue-225-slide-deck.md`. The plan structure paid off — small bite-sized tasks made the per-task commits clean, and the spec's "Divergence from prior art" section flagged the committed-`slides.html` decision up front instead of discovering it mid-render.

**Headline numbers unchanged from notebook re-run on 2026-05-21:** F1=0.300, % AG-unique vs GTEx = 0.0%, decision = NO-GO. The deck just re-renders the same story in lab-seminar form.

**Ready to merge after this commit lands.** Closure-ritual gate (test plan + Issue #225 ACs + #203 carrier step) will re-check before invoking `bash scripts/audit_and_merge.sh 452`.

---
```

Replace `<TIMESTAMP_FROM_STEP_1>` with the actual timestamp from Step 1.

- [ ] **Step 4: Commit**

```bash
git add research/lab_notebook/scientist.md
git commit -m "$(cat <<'EOF'
docs(scientist): lab notebook 2026-05-22 — #225 deck shipped in PR #452

Captures scope shipped here (vs original #455 deferral), the
slides-co-location convention formalised into CLAUDE.md, and the
brainstorming/writing-plans/executing-plans flow that produced it.

Headline numbers unchanged (F1=0.300, NO-GO); deck re-renders the
notebook story in lab-seminar form for human review.

Refs #225

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: Push + tick PR #452 Test plan boxes

**Files:**
- Modify: PR #452 body (external)

- [ ] **Step 1: Confirm with user before pushing**

The branch has accumulated 6 commits (spec + refs + regenerator + figures + slides source + slides.html + CLAUDE.md + lab notebook). Surface a summary and **wait for the user's explicit OK before `git push`** (per the never-chain commit/push memory rule).

Run:
```bash
git log --oneline origin/research/scientist/issue-225-normal-junction-filter-strength..HEAD
git diff --stat origin/research/scientist/issue-225-normal-junction-filter-strength..HEAD
```

Surface both to the user with the message: *"Ready to push N commits to the PR branch. OK to `git push`?"* Wait for user's `yes` / `go` / equivalent.

- [ ] **Step 2: Push**

After user OK:
```bash
git push
```

- [ ] **Step 3: Re-fetch PR #452 body and tick the deck Test plan items**

```bash
gh pr view 452 --json body --jq .body > /tmp/pr_452_body.md
# Manually edit /tmp/pr_452_body.md: change each `- [ ]` on the 6 deck-related items to `- [x]`
# Items to tick (all 6 added in Task 9):
#   - slides.qmd renders cleanly
#   - _regenerate_figures.py produces both PNGs
#   - figure refs resolve
#   - CLAUDE.md updated
#   - Issue #455 descoped
#   - Lab notebook entry
gh pr edit 452 --body-file /tmp/pr_452_body.md
```

- [ ] **Step 4: Verify all Test plan items ticked**

```bash
gh pr view 452 --json body --jq .body | grep -c "^\- \[ \]"
```

Expected: `0` — no unticked boxes remain.

- [ ] **Step 5: Verify [Issue #225](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/225) ACs still all ticked**

```bash
gh issue view 225 --json body --jq .body | grep -c "^\- \[ \]"
```

Expected: `0`. (They were already all `[x]` before this deck work; the deck doesn't introduce new ACs.)

---

## Task 12: Closure-ritual gate + (optional) merge

**This step is gated on explicit user direction.** Per the never-chain memory rule, do not run `audit_and_merge.sh` without the user saying go.

- [ ] **Step 1: Confirm merge intent with user**

Ask: *"All 6 deck Test plan items now ticked, Issue #225 ACs clean, push complete. Ready to run `bash scripts/audit_and_merge.sh 452` to merge?"*

- [ ] **Step 2: On user OK, run the gate**

```bash
bash scripts/audit_and_merge.sh 452
```

Defaults: `--squash --delete-branch`. The script audits PR + linked Issue #225 ACs; if any unticked box remains it prints to stderr and exits 1 without merging.

Expected on a clean audit: PR merges, branch deleted, Issue #225 auto-closes via "Closes" reference in the PR body.

- [ ] **Step 3: Immediate post-merge — update Issue #203 Exp 3 row**

Per [PR #452](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/452) Test plan: the "Issue #203 Exp 3 row update" carrier step fires now.

```bash
gh issue view 203 --json body --jq .body > /tmp/issue_203_body.md
# Edit /tmp/issue_203_body.md — find the Exp 3 row, write in the NO-GO verdict + headline numbers:
#   F1=0.300, % AG-unique vs GTEx = 0.0%, decision = NO-GO — treat as tissue prior
# Reference the merged PR: "See PR #452 + research/experiments/issue_225_*"
gh issue edit 203 --body-file /tmp/issue_203_body.md
```

- [ ] **Step 4: Flip Issue #203 board Status if appropriate**

If all three of [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203)'s decision-rule experiments now have outcomes (Exp 1 done in #224, Exp 3 just done here, Exp 2 deferred to #381), [Issue #203](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/203) is ready to enter `In review` — check the parent before flipping.

---

## Self-review (run before declaring plan complete)

**Spec coverage check:**

| Spec scope item | Implementation task |
|---|---|
| Deck source at `slides.qmd` | Task 5 |
| Figure regenerator | Task 3 |
| Two new figures `pr_curve.png` + `caught_bar.png` | Task 4 |
| Re-used Venn from `outputs/` | Task 5 Step 1 (referenced; no file copy) |
| Local `refs.bib` | Task 2 |
| Rendered `slides.html` committed | Task 6 |
| CLAUDE.md slides-co-location update | Task 7 |
| Issue #455 body update | Task 8 |
| PR #452 Test plan amendment | Task 9 |
| Lab notebook entry | Task 10 |

All 10 spec items covered.

**Placeholder scan:** no TBD, TODO, "implement later", "fill in" placeholders in any task. Code blocks are complete. Commands have expected output. Test plan items are explicit. ✓

**Type / identifier consistency:**
- `_regenerate_figures.py` (Task 3) — name matches prior art; called by name in Task 4
- `make_pr_curve()` / `make_caught_bar()` (Task 3) — function names; `main()` calls both
- BibTeX keys (Task 2): `cotto2023regtools`, `kim2019hisat2`, `mudge2025gencode`, `pedregosa2011scikit`, `avsec2026alphagenome`, `wilks2018snaptron`, `gtex2020v8` — all 7 referenced from `slides.qmd` (Task 5)
- `pr_curve.png`, `caught_bar.png` filenames consistent across Tasks 3 / 4 / 5
- `filter_venn_chr22.png` path: lives in `outputs/`, referenced from `slides.qmd` as `outputs/filter_venn_chr22.png` (1 level down from deck) — consistent

All consistent.
