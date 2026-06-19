# #708 KDE + Centered-Isotonic Immunogenicity Calibrator — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone, importable calibrator that maps MHCflurry `genotype_presentation_score` → `calibrated_immunogenicity_log_odds` via adaptive KDE → density-ratio + true-base-rate prior → isotonic → centered-isotonic regression, validated with leave-one-cohort-out on the #707 harmonized cohorts.

**Architecture:** Approach A — all code under `research/experiments/issue_547_immunogenicity_calibration/`: a pure importable module (`calibrator.py`), a one-time cached MHCflurry precompute (`score_cohort.py`), an orchestration notebook (`notebook.ipynb`), and local pytest. #709 (Developer) consumes the module + serialized artifact and decides production promotion.

**Tech Stack:** Python 3.14 (`research/.venv`), numpy, scipy, scikit-learn, pandas, joblib, matplotlib, mhcflurry (`Class1PresentationPredictor`), pytest; R + `cir` package used **once** to generate frozen test fixtures.

## Global Constraints

- Reimplement NeoGuider's calibration **from the paper** — do NOT vendor code from `github.com/XuegongLab/neoguider` (AGPL-3.0 + non-commercial; paid-commercial deps). Verbatim from spec.
- Output column name is exactly `calibrated_immunogenicity_log_odds` (a presentation-anchored prior, **not** a validated probability).
- Do NOT modify `workflow/scripts/run_mhcflurry.py` — replicate + test-pin its genotype formula instead.
- Genotype score formula (verbatim): `genotype_presentation_score = 1 − ∏ᵢ(1 − wᵢ·pᵢ)`, `pᵢ` = per-allele `presentation_score`, `wᵢ` = HLA-C weight for HLA-C alleles else 1.0.
- Label mapping: NeoRanking `response_type` `CD8`→1 / `negative`→0 / drop `not_tested`; IMPROVE `response` already 1/0.
- HLA canonical form for MHCflurry: `HLA-A*01:01` (normalize NeoRanking `A*01:01` and IMPROVE `HLA-A01:01`).
- Prior odds use **true pre-subsample** assayed counts, never subsample counts.
- Keep ALL positives in any subsample; only negatives may be thinned.
- Tests run locally with `research/.venv/bin/python -m pytest` (these are research-folder tests, distinct from the `workflow/tests/` pipeline suite; not wired into CI in this issue).
- Commit frequently; never chain commit+push (surface SHA, await OK before pushing).

## File Structure

```
research/experiments/issue_547_immunogenicity_calibration/
├── calibrator.py            # centered_isotonic() + PresentationCalibrator (fit/transform/save/load)
├── score_cohort.py          # MHCflurry precompute → outputs/scored_cohort_subsample.parquet (cached)
├── notebook.ipynb           # LOCO + within-cohort k-fold + diagnostics + final artifact
├── tests/
│   ├── test_centered_isotonic.py   # port vs frozen R cir fixtures + edge cases
│   ├── test_calibrator.py          # monotonicity, prior, clipping, save/load
│   ├── test_score_formula.py       # genotype combine formula vs run_mhcflurry.py
│   ├── generate_cir_fixtures.R     # one-time R cir oracle generator (documents R/cir version)
│   └── fixtures/cir_fixtures.json  # frozen R cir outputs (committed)
└── outputs/
    ├── scored_cohort_subsample.parquet
    ├── calibrator_v1.joblib
    └── *.png                       # reliability / monotonicity / shift-gap / KDE-compare
research/requirements.txt    # ensure: scipy, scikit-learn, joblib, pytest, mhcflurry
```

The branch `research/scientist/issue-708-immunogenicity-calibrator` already exists with the design spec committed.

---

### Task 1: Centered-isotonic port + R-oracle fixtures

**Files:**
- Create: `research/experiments/issue_547_immunogenicity_calibration/calibrator.py`
- Create: `research/experiments/issue_547_immunogenicity_calibration/tests/generate_cir_fixtures.R`
- Create: `research/experiments/issue_547_immunogenicity_calibration/tests/fixtures/cir_fixtures.json`
- Create: `research/experiments/issue_547_immunogenicity_calibration/tests/test_centered_isotonic.py`

**Interfaces:**
- Produces: `centered_isotonic(x, y_iso, w) -> (cx: np.ndarray, cy: np.ndarray)` — collapses each flat level-set of an isotonic fit to its weight-centroid x, keeps its y; returns sorted centred knots.

- [ ] **Step 1: Ensure research env + test deps**

Run (one-time, per-clone):
```bash
cd research && pyenv local 3.14.4 && python -m venv .venv && cd ..
research/.venv/bin/pip install -r research/requirements.txt
research/.venv/bin/pip install scipy scikit-learn joblib pytest   # if absent from requirements.txt
research/.venv/bin/python -c "import numpy, scipy, sklearn, joblib, pytest; print('ok')"
```
Expected: `ok`. Add any newly-installed packages to `research/requirements.txt`.

- [ ] **Step 2: Generate the R `cir` oracle fixtures**

Write `tests/generate_cir_fixtures.R`:
```r
# Generates frozen centered-isotonic-regression oracle outputs from the R `cir`
# package for the Python port test. Run once: Rscript generate_cir_fixtures.R
# Requires: install.packages("cir"); jsonlite. Records package versions.
library(cir); library(jsonlite)
cases <- list(
  # Oron & Flournoy 2017 style: monotone-with-plateau dose-response
  onf_example = list(x=c(1,2,3,4,5,6), y=c(0.1,0.2,0.2,0.2,0.5,0.9), w=c(10,10,10,10,10,10)),
  multi_levelset = list(x=c(1,2,3,4,5,6,7,8), y=c(0.1,0.1,0.1,0.4,0.4,0.7,0.7,0.7), w=c(5,8,3,10,6,4,9,2)),
  boundary_plateau = list(x=c(1,2,3,4,5), y=c(0.3,0.3,0.3,0.6,0.9), w=c(7,7,7,7,7)),
  no_collapse = list(x=c(1,2,3,4), y=c(0.1,0.3,0.6,0.9), w=c(4,4,4,4))
)
out <- list(r_version=R.version.string, cir_version=as.character(packageVersion("cir")), cases=list())
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- cirPAVA(y=cs$y, x=cs$x, wt=cs$w, full=TRUE)  # centered isotonic (CIR)
  out$cases[[nm]] <- list(x=cs$x, y=cs$y, w=cs$w, cir_x=as.numeric(fit$x), cir_y=as.numeric(fit$y))
}
write_json(out, "fixtures/cir_fixtures.json", auto_unbox=TRUE, digits=10, pretty=TRUE)
```
Run:
```bash
cd research/experiments/issue_547_immunogenicity_calibration/tests
Rscript generate_cir_fixtures.R   # requires R + install.packages("cir")
```
Expected: `fixtures/cir_fixtures.json` written with 4 cases + version strings. **If R/`cir` is unavailable**, hand-enter the Oron & Flournoy 2017 published worked-example knots into `fixtures/cir_fixtures.json` for the `onf_example` case and note the provenance in the file header comment; keep the R script committed for future regeneration.

- [ ] **Step 3: Write the failing test**

`tests/test_centered_isotonic.py`:
```python
import json
from pathlib import Path
import numpy as np
import pytest
from sklearn.isotonic import IsotonicRegression
from calibrator import centered_isotonic

FIX = json.loads((Path(__file__).parent / "fixtures" / "cir_fixtures.json").read_text())

@pytest.mark.parametrize("name", list(FIX["cases"].keys()))
def test_matches_r_cir(name):
    c = FIX["cases"][name]
    x, y, w = np.array(c["x"], float), np.array(c["y"], float), np.array(c["w"], float)
    y_iso = IsotonicRegression(increasing=True, out_of_bounds="clip").fit_transform(x, y, sample_weight=w)
    cx, cy = centered_isotonic(x, y_iso, w)
    # the port must reproduce R cir's centred knots within tolerance
    np.testing.assert_allclose(cx, np.array(c["cir_x"], float), atol=1e-6)
    np.testing.assert_allclose(cy, np.array(c["cir_y"], float), atol=1e-6)

def test_centred_knots_are_monotone():
    x = np.array([1, 2, 3, 4, 5, 6.])
    y_iso = np.array([0.1, 0.2, 0.2, 0.2, 0.5, 0.9])
    w = np.ones_like(x)
    cx, cy = centered_isotonic(x, y_iso, w)
    assert np.all(np.diff(cy) >= 0)
    assert np.all(np.diff(cx) > 0)
```

- [ ] **Step 4: Run test to verify it fails**

Run: `cd research/experiments/issue_547_immunogenicity_calibration && ../../../research/.venv/bin/python -m pytest tests/test_centered_isotonic.py -v`
Expected: FAIL — `ImportError: cannot import name 'centered_isotonic'`.

- [ ] **Step 5: Implement `centered_isotonic` in `calibrator.py`**

```python
"""Immunogenicity calibrator: KDE -> density-ratio + prior -> isotonic ->
centered isotonic (Oron & Flournoy 2017). Reimplemented from the NeoGuider
paper (Wei et al. 2026); no NeoGuider code vendored."""
import numpy as np


def centered_isotonic(x, y_iso, w):
    """Centered isotonic regression (Oron & Flournoy 2017).

    Collapse each flat level-set of the isotonic fit `y_iso` to a single knot
    whose x is the weight-centroid of the pooled level-set points; keep the
    isotonic y. Returns sorted (centre_x, centre_y).
    """
    x = np.asarray(x, float); y_iso = np.asarray(y_iso, float); w = np.asarray(w, float)
    order = np.argsort(x)
    x, y_iso, w = x[order], y_iso[order], w[order]
    cx, cy = [], []
    i, n = 0, len(x)
    while i < n:
        j = i
        while j + 1 < n and np.isclose(y_iso[j + 1], y_iso[i]):
            j += 1
        sl = slice(i, j + 1)
        wsum = w[sl].sum()
        cx.append((x[sl] * w[sl]).sum() / wsum if wsum > 0 else x[sl].mean())
        cy.append(y_iso[i])
        i = j + 1
    return np.array(cx), np.array(cy)
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `cd research/experiments/issue_547_immunogenicity_calibration && ../../../research/.venv/bin/python -m pytest tests/test_centered_isotonic.py -v`
Expected: PASS (all parametrized cases + monotonicity).

- [ ] **Step 7: Commit**

```bash
git add research/experiments/issue_547_immunogenicity_calibration/calibrator.py \
        research/experiments/issue_547_immunogenicity_calibration/tests/ research/requirements.txt
git commit -m "feat(scoring): centered-isotonic port + R cir oracle fixtures (#708)"
```

---

### Task 2: `PresentationCalibrator.fit` / `.transform`

**Files:**
- Modify: `research/experiments/issue_547_immunogenicity_calibration/calibrator.py`
- Create: `research/experiments/issue_547_immunogenicity_calibration/tests/test_calibrator.py`

**Interfaces:**
- Consumes: `centered_isotonic` (Task 1).
- Produces:
  - `adaptive_kde(samples) -> callable(grid)->density` (Abramson variable-bandwidth) and `fixed_kde(samples) -> callable`.
  - `class PresentationCalibrator(kde_mode="adaptive", n_grid=512)` with `fit(scores, labels, n_pos_true=None, n_neg_true=None) -> self` and `transform(scores) -> np.ndarray` (the `calibrated_immunogenicity_log_odds`). After fit: attributes `cx_, cy_, prior_, score_range_`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_calibrator.py`:
```python
import numpy as np
import pytest
from calibrator import PresentationCalibrator

def _synthetic(seed=0, n_pos=200, n_neg=4000):
    rng = np.random.default_rng(seed)
    pos = np.clip(rng.beta(6.0, 2.2, n_pos), 0, 1)   # immunogenic skew high
    neg = np.clip(rng.beta(2.0, 4.0, n_neg), 0, 1)   # non-immunogenic skew low
    scores = np.concatenate([pos, neg])
    labels = np.concatenate([np.ones(n_pos), np.zeros(n_neg)]).astype(int)
    return scores, labels

@pytest.mark.parametrize("mode", ["adaptive", "fixed"])
def test_transform_is_monotone_nondecreasing(mode):
    scores, labels = _synthetic()
    cal = PresentationCalibrator(kde_mode=mode).fit(scores, labels)
    grid = np.linspace(0, 1, 200)
    out = cal.transform(grid)
    assert np.all(np.diff(out) >= -1e-9)

def test_higher_score_higher_logodds():
    scores, labels = _synthetic()
    cal = PresentationCalibrator().fit(scores, labels)
    assert cal.transform([0.9])[0] > cal.transform([0.1])[0]

def test_out_of_range_clips_to_boundary():
    scores, labels = _synthetic()
    cal = PresentationCalibrator().fit(scores, labels)
    lo, hi = cal.score_range_
    assert cal.transform([hi + 5.0])[0] == pytest.approx(cal.transform([hi])[0])
    assert cal.transform([lo - 5.0])[0] == pytest.approx(cal.transform([lo])[0])

def test_prior_uses_true_counts_not_subsample():
    scores, labels = _synthetic(n_pos=200, n_neg=4000)
    # same data, but declare the TRUE base rate as much rarer
    cal_sub = PresentationCalibrator().fit(scores, labels)
    cal_true = PresentationCalibrator().fit(scores, labels, n_pos_true=200, n_neg_true=400000)
    # rarer true prior => uniformly lower log-odds intercept
    g = np.linspace(0.2, 0.8, 50)
    assert np.all(cal_true.transform(g) < cal_sub.transform(g))
```

- [ ] **Step 2: Run to verify failure**

Run: `cd research/experiments/issue_547_immunogenicity_calibration && ../../../research/.venv/bin/python -m pytest tests/test_calibrator.py -v`
Expected: FAIL — `ImportError: cannot import name 'PresentationCalibrator'`.

- [ ] **Step 3: Implement KDEs + fit/transform**

Append to `calibrator.py`:
```python
from scipy.stats import gaussian_kde
from sklearn.isotonic import IsotonicRegression

_EPS = 1e-12


def fixed_kde(samples):
    """Fixed-bandwidth Gaussian KDE (Scott's rule) — the baseline."""
    kde = gaussian_kde(np.asarray(samples, float))
    return lambda grid: kde(np.atleast_1d(np.asarray(grid, float)))


def adaptive_kde(samples):
    """Abramson sample-point variable-bandwidth Gaussian KDE.

    Local bandwidth ∝ 1/sqrt(pilot density), normalised to the geometric mean.
    """
    s = np.asarray(samples, float)
    pilot = gaussian_kde(s)
    f_pilot = np.clip(pilot(s), _EPS, None)
    g = np.exp(np.mean(np.log(f_pilot)))
    lam = np.sqrt(g / f_pilot)                 # per-sample local factor
    std = np.std(s, ddof=1) if len(s) > 1 else 1.0
    bw = pilot.factor * std * lam              # per-sample bandwidth
    bw = np.clip(bw, _EPS, None)

    def density(grid):
        x = np.atleast_1d(np.asarray(grid, float))
        u = (x[:, None] - s[None, :]) / bw[None, :]
        k = np.exp(-0.5 * u ** 2) / (np.sqrt(2 * np.pi) * bw[None, :])
        return k.mean(axis=1)
    return density


class PresentationCalibrator:
    """Calibrate genotype_presentation_score -> calibrated_immunogenicity_log_odds."""

    def __init__(self, kde_mode="adaptive", n_grid=512):
        if kde_mode not in ("adaptive", "fixed"):
            raise ValueError("kde_mode must be 'adaptive' or 'fixed'")
        self.kde_mode = kde_mode
        self.n_grid = n_grid

    def fit(self, scores, labels, n_pos_true=None, n_neg_true=None, fit_cohorts=None):
        scores = np.asarray(scores, float)
        labels = np.asarray(labels, int)
        pos, neg = scores[labels == 1], scores[labels == 0]
        if len(pos) < 2 or len(neg) < 2:
            raise ValueError("need >=2 samples per class")
        n_pos_true = len(pos) if n_pos_true is None else n_pos_true
        n_neg_true = len(neg) if n_neg_true is None else n_neg_true
        self.prior_ = float(np.log(n_pos_true / n_neg_true))

        make = adaptive_kde if self.kde_mode == "adaptive" else fixed_kde
        kde_pos, kde_neg = make(pos), make(neg)

        lo, hi = float(scores.min()), float(scores.max())
        self.score_range_ = (lo, hi)
        grid = np.linspace(lo, hi, self.n_grid)
        d_pos = np.clip(kde_pos(grid), _EPS, None)
        d_neg = np.clip(kde_neg(grid), _EPS, None)
        raw = np.log(d_pos) - np.log(d_neg) + self.prior_         # log-odds(grid)
        # local data support as isotonic weights (true-count-scaled densities)
        w = n_pos_true * d_pos + n_neg_true * d_neg
        iso = IsotonicRegression(increasing=True, out_of_bounds="clip")
        y_iso = iso.fit_transform(grid, raw, sample_weight=w)
        self.cx_, self.cy_ = centered_isotonic(grid, y_iso, w)
        self.fit_cohorts_ = list(fit_cohorts) if fit_cohorts else None
        return self

    def transform(self, scores):
        scores = np.atleast_1d(np.asarray(scores, float))
        return np.interp(scores, self.cx_, self.cy_)  # np.interp clips to endpoints
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd research/experiments/issue_547_immunogenicity_calibration && ../../../research/.venv/bin/python -m pytest tests/test_calibrator.py -v`
Expected: PASS (4 tests / 5 cases incl. both KDE modes).

- [ ] **Step 5: Commit**

```bash
git add research/experiments/issue_547_immunogenicity_calibration/calibrator.py \
        research/experiments/issue_547_immunogenicity_calibration/tests/test_calibrator.py
git commit -m "feat(scoring): PresentationCalibrator fit/transform with adaptive+fixed KDE (#708)"
```

---

### Task 3: Calibrator serialization (`save`/`load`)

**Files:**
- Modify: `research/experiments/issue_547_immunogenicity_calibration/calibrator.py`
- Modify: `research/experiments/issue_547_immunogenicity_calibration/tests/test_calibrator.py`

**Interfaces:**
- Consumes: a fitted `PresentationCalibrator` (Task 2).
- Produces: `cal.save(path)` (joblib dict: `cx, cy, prior, score_range, kde_mode, fit_cohorts, version`) and classmethod `PresentationCalibrator.load(path) -> PresentationCalibrator` whose `transform` reproduces the original exactly.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_calibrator.py`:
```python
def test_save_load_roundtrip(tmp_path):
    scores, labels = _synthetic()
    cal = PresentationCalibrator().fit(scores, labels, fit_cohorts=["NCI", "TESLA"])
    p = tmp_path / "cal.joblib"
    cal.save(p)
    cal2 = PresentationCalibrator.load(p)
    g = np.linspace(0, 1, 100)
    np.testing.assert_allclose(cal.transform(g), cal2.transform(g))
    assert cal2.prior_ == cal.prior_
    assert cal2.fit_cohorts_ == ["NCI", "TESLA"]
```

- [ ] **Step 2: Run to verify failure**

Run: `cd research/experiments/issue_547_immunogenicity_calibration && ../../../research/.venv/bin/python -m pytest tests/test_calibrator.py::test_save_load_roundtrip -v`
Expected: FAIL — `AttributeError: 'PresentationCalibrator' object has no attribute 'save'`.

- [ ] **Step 3: Implement save/load**

Add `import joblib` near the top of `calibrator.py`, then append methods to the class:
```python
    CALIBRATOR_VERSION = "v1"

    def save(self, path):
        joblib.dump({
            "cx": self.cx_, "cy": self.cy_, "prior": self.prior_,
            "score_range": self.score_range_, "kde_mode": self.kde_mode,
            "fit_cohorts": self.fit_cohorts_, "version": self.CALIBRATOR_VERSION,
        }, path)

    @classmethod
    def load(cls, path):
        d = joblib.load(path)
        obj = cls(kde_mode=d["kde_mode"])
        obj.cx_, obj.cy_ = d["cx"], d["cy"]
        obj.prior_ = d["prior"]
        obj.score_range_ = tuple(d["score_range"])
        obj.fit_cohorts_ = d["fit_cohorts"]
        return obj
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd research/experiments/issue_547_immunogenicity_calibration && ../../../research/.venv/bin/python -m pytest tests/test_calibrator.py -v`
Expected: PASS (all calibrator tests).

- [ ] **Step 5: Commit**

```bash
git add research/experiments/issue_547_immunogenicity_calibration/calibrator.py \
        research/experiments/issue_547_immunogenicity_calibration/tests/test_calibrator.py
git commit -m "feat(scoring): calibrator joblib save/load (knots + metadata) (#708)"
```

---

### Task 4: Cohort scoring precompute (`score_cohort.py`)

**Files:**
- Create: `research/experiments/issue_547_immunogenicity_calibration/score_cohort.py`
- Create: `research/experiments/issue_547_immunogenicity_calibration/tests/test_score_formula.py`
- Modify: `research/requirements.txt` (add `mhcflurry`)

**Interfaces:**
- Produces:
  - `genotype_score(per_allele_scores, alleles, hla_c_weight=0.5) -> float` — the `1 − ∏(1 − wᵢ·pᵢ)` combiner (test-pinned against `run_mhcflurry.py`'s documented formula).
  - `normalize_hla(allele) -> str` → canonical `HLA-A*01:01`.
  - `build_scored_cohort(out_parquet, n_neg=50000, seed=0) -> pd.DataFrame` — loads cohorts, filters/labels, subsamples (all positives + stratified negatives), scores via `Class1PresentationPredictor`, writes parquet with true-count metadata; skip-if-cached.

- [ ] **Step 1: Install mhcflurry + models**

Run:
```bash
research/.venv/bin/pip install mhcflurry
research/.venv/bin/python -m mhcflurry.downloads fetch models_class1_presentation
echo "mhcflurry" >> research/requirements.txt   # if not already present
```
Expected: models downloaded to MHCflurry's platformdirs cache.

- [ ] **Step 2: Write the failing formula + HLA tests**

`tests/test_score_formula.py`:
```python
import numpy as np
import pytest
from score_cohort import genotype_score, normalize_hla

def test_genotype_score_matches_breadth_formula():
    # 1 - prod(1 - w_i * p_i); HLA-C weighted at 0.5
    p = {"HLA-A*02:01": 0.8, "HLA-B*07:02": 0.3, "HLA-C*07:02": 0.6}
    expected = 1.0 - ((1 - 1.0 * 0.8) * (1 - 1.0 * 0.3) * (1 - 0.5 * 0.6))
    got = genotype_score(p, hla_c_weight=0.5)
    assert got == pytest.approx(expected, abs=1e-9)

def test_genotype_score_single_allele():
    assert genotype_score({"HLA-A*02:01": 0.5}) == pytest.approx(0.5, abs=1e-9)

@pytest.mark.parametrize("raw,canon", [
    ("A*01:01", "HLA-A*01:01"),
    ("HLA-A01:01", "HLA-A*01:01"),
    ("HLA-A*01:01", "HLA-A*01:01"),
    ("C*07:02", "HLA-C*07:02"),
])
def test_normalize_hla(raw, canon):
    assert normalize_hla(raw) == canon
```

- [ ] **Step 3: Run to verify failure**

Run: `cd research/experiments/issue_547_immunogenicity_calibration && ../../../research/.venv/bin/python -m pytest tests/test_score_formula.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'score_cohort'`.

- [ ] **Step 4: Implement the helpers + builder**

`score_cohort.py`:
```python
"""One-time MHCflurry precompute: score the assayed calibration cohorts'
peptides at the genotype level (genotype_presentation_score) and cache to
parquet. Subsample-first (all positives + stratified negatives). Cached.

Run: ../../../research/.venv/bin/python score_cohort.py
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "outputs" / "scored_cohort_subsample.parquet"
HLA_C_WEIGHT = 0.5   # mirror config default used by run_mhcflurry.py


def genotype_score(per_allele_scores, hla_c_weight=HLA_C_WEIGHT):
    """1 - prod(1 - w_i * p_i); w_i = hla_c_weight for HLA-C else 1.0.
    Mirrors workflow/scripts/run_mhcflurry.py (do not import it)."""
    prod = 1.0
    for allele, p in per_allele_scores.items():
        w = hla_c_weight if allele.startswith("HLA-C") else 1.0
        prod *= (1.0 - w * float(p))
    return round(1.0 - prod, 6)


def normalize_hla(allele):
    """-> canonical 'HLA-A*01:01'."""
    a = allele.strip()
    if not a.startswith("HLA-"):
        a = "HLA-" + a
    # ensure the gene/number separator '*' after the gene letter(s)
    m = re.match(r"^HLA-([A-Z]+\d?)\*?(\d+:\d+.*)$", a)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}"
    return a


def _load_neoranking():
    df = pd.read_csv(DATA / "neoranking" / "Neopep_data_org.txt", sep="\t",
                     usecols=["mutant_seq", "response_type", "mutant_best_alleles", "dataset", "patient"])
    df = df[df["response_type"].isin(["CD8", "negative"])].copy()
    df["label"] = (df["response_type"] == "CD8").astype(int)
    df["cohort"] = df["dataset"]
    df = df.rename(columns={"mutant_seq": "peptide"})
    return df[["peptide", "label", "cohort", "patient", "mutant_best_alleles"]]


def _load_improve():
    df = pd.read_csv(DATA / "improve_borch" / "In_house_neoepitope_for_CV.tsv", sep="\t")
    df = df.rename(columns={"Mut_peptide": "peptide", "HLA_allele": "mutant_best_alleles",
                            "response": "label"})
    df["cohort"] = "IMPROVE"
    df["patient"] = "IMPROVE"   # no per-patient genotype; score single best allele
    return df[["peptide", "label", "cohort", "patient", "mutant_best_alleles"]]


def _hla_genotypes():
    raw = pd.read_csv(DATA / "neoranking" / "HLA_allotypes.txt", sep="\t",
                      header=0, names=["patient", "alleles"])
    return {r.patient: [normalize_hla(a) for a in str(r.alleles).split(",")]
            for r in raw.itertuples(index=False)}


def build_scored_cohort(out_parquet=OUT, n_neg=50000, seed=0):
    if Path(out_parquet).exists():
        return pd.read_parquet(out_parquet)
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()

    cohort = pd.concat([_load_neoranking(), _load_improve()], ignore_index=True)
    genotypes = _hla_genotypes()

    # true pre-subsample counts (for the calibrator prior)
    true_counts = (cohort.groupby(["cohort", "label"]).size()
                   .unstack(fill_value=0).rename(columns={0: "n_neg", 1: "n_pos"}))

    pos = cohort[cohort["label"] == 1]
    neg = cohort[cohort["label"] == 0]
    # cohort-stratified negative subsample, fixed seed; keep ALL positives
    frac = min(1.0, n_neg / max(len(neg), 1))
    neg_sub = neg.groupby("cohort", group_keys=False).apply(
        lambda g: g.sample(frac=frac, random_state=seed))
    work = pd.concat([pos, neg_sub], ignore_index=True)

    rows = []
    for r in work.itertuples(index=False):
        alleles = genotypes.get(r.patient, [normalize_hla(r.mutant_best_alleles)])
        # per-allele presentation_score via MHCflurry (genotype call)
        pred = predictor.predict(peptides=[r.peptide], alleles={a: [a] for a in alleles},
                                 verbose=0)
        # one row per allele in pred -> map allele -> presentation_score
        per_allele = {normalize_hla(a): s for a, s in
                      zip(pred["best_allele"], pred["presentation_score"])} \
            if "best_allele" in pred else {}
        # robust path: score each allele individually
        per_allele = {}
        for a in alleles:
            pp = predictor.predict(peptides=[r.peptide], alleles={a: [a]}, verbose=0)
            per_allele[a] = float(pp["presentation_score"].iloc[0])
        rows.append({
            "peptide": r.peptide, "cohort": r.cohort, "label": int(r.label),
            "genotype_presentation_score": genotype_score(per_allele),
            "best_presentation_score": max(per_allele.values()),
        })
    out = pd.DataFrame(rows)
    Path(out_parquet).parent.mkdir(parents=True, exist_ok=True)
    out.attrs["true_counts"] = true_counts.to_dict()
    out.to_parquet(out_parquet)
    true_counts.to_csv(Path(out_parquet).with_suffix(".true_counts.csv"))
    return out


if __name__ == "__main__":
    df = build_scored_cohort()
    print(df.groupby(["cohort", "label"]).size())
    print(df["genotype_presentation_score"].describe())
```

- [ ] **Step 5: Run the formula/HLA tests to verify they pass**

Run: `cd research/experiments/issue_547_immunogenicity_calibration && ../../../research/.venv/bin/python -m pytest tests/test_score_formula.py -v`
Expected: PASS (formula + single-allele + 4 HLA cases). *(These tests don't import mhcflurry — only the pure helpers.)*

- [ ] **Step 6: Run the actual precompute (one-time, ~30–90 min CPU)**

Run: `cd research/experiments/issue_547_immunogenicity_calibration && ../../../research/.venv/bin/python score_cohort.py`
Expected: prints per-cohort label counts (all positives kept: NCI/TESLA/HiTIDE CD8 + IMPROVE pos; ~50K negatives total) and a `genotype_presentation_score` summary; writes `outputs/scored_cohort_subsample.parquet` + `.true_counts.csv`.

- [ ] **Step 7: Commit**

```bash
git add research/experiments/issue_547_immunogenicity_calibration/score_cohort.py \
        research/experiments/issue_547_immunogenicity_calibration/tests/test_score_formula.py \
        research/experiments/issue_547_immunogenicity_calibration/outputs/scored_cohort_subsample.true_counts.csv \
        research/requirements.txt
# NOTE: the .parquet is gitignored (data/) — confirm it is NOT staged; commit only the small .true_counts.csv
git commit -m "feat(scoring): cohort MHCflurry precompute + genotype-score formula test (#708)"
```

---

### Task 5: Validation notebook + final artifact

**Files:**
- Create: `research/experiments/issue_547_immunogenicity_calibration/notebook.ipynb`
- Create (notebook outputs): `outputs/calibrator_v1.joblib`, `outputs/*.png`

**Interfaces:**
- Consumes: `outputs/scored_cohort_subsample.parquet` (Task 4), `PresentationCalibrator` (Tasks 2–3).
- Produces: `outputs/calibrator_v1.joblib` (final fit on all 4 cohorts) + reliability / monotonicity / shift-gap / KDE-compare PNGs.

Build the notebook cell-by-cell (narrate each cell before writing per the notebook-narration rule). Cells, in order:

- [ ] **Step 1: Setup cell** — imports (`numpy, pandas, matplotlib, joblib`, `from calibrator import PresentationCalibrator`), load `outputs/scored_cohort_subsample.parquet` + `.true_counts.csv`, print per-cohort label counts. Markdown lead: goal + the proxy-not-direct-test caveat (copy from spec).

- [ ] **Step 2: LOCO cell (primary).** For each held-out cohort in `["NCI","TESLA","HiTIDE","IMPROVE"]`: fit `PresentationCalibrator` on the other three (pass true `n_pos_true/n_neg_true` summed over the training cohorts), `transform` the held-out scores, compute reliability (bin predicted log-odds → observed positive rate, with Wilson CIs), monotonicity check, and AUPRC + positives-in-top-20. Collect into a results table.

- [ ] **Step 3: Within-cohort k-fold contrast cell.** Stratified 5-fold within each cohort; same metrics. Compute and display the **LOCO-vs-within-cohort gap** per cohort (the cohort-shift penalty).

- [ ] **Step 4: KDE comparison cell.** Re-run LOCO with `kde_mode="fixed"` vs `"adaptive"`; overlay reliability curves; one sentence on whether adaptive helps.

- [ ] **Step 5: Diagnostics figures cell.** Save `outputs/pr_reliability.png`, `outputs/monotonicity.png`, `outputs/shift_gap.png`, `outputs/kde_compare.png`.

- [ ] **Step 6: Final-fit cell.** Fit `PresentationCalibrator(kde_mode=<winner>)` on ALL four cohorts with the full true counts, `cal.save("outputs/calibrator_v1.joblib")`. Print the knots + prior + score range.

- [ ] **Step 7: Execute the notebook end-to-end.**

Run: `cd research/experiments/issue_547_immunogenicity_calibration && ../../../research/.venv/bin/jupyter nbconvert --to notebook --execute --inplace notebook.ipynb`
Expected: executes with no errors; `outputs/calibrator_v1.joblib` + 4 PNGs written.

- [ ] **Step 8: Sanity-check the artifact.**

Run:
```bash
cd research/experiments/issue_547_immunogenicity_calibration
../../../research/.venv/bin/python -c "from calibrator import PresentationCalibrator as C; import numpy as np; c=C.load('outputs/calibrator_v1.joblib'); g=np.linspace(0,1,5); print('cohorts',c.fit_cohorts_,'prior',round(c.prior_,3)); print(list(zip(g.round(2), c.transform(g).round(3))))"
```
Expected: prints the 4 cohorts, the (negative) prior, and a monotone non-decreasing log-odds sequence across the grid.

- [ ] **Step 9: Commit**

```bash
git add research/experiments/issue_547_immunogenicity_calibration/notebook.ipynb \
        research/experiments/issue_547_immunogenicity_calibration/outputs/calibrator_v1.joblib \
        research/experiments/issue_547_immunogenicity_calibration/outputs/*.png
# (calibrator_v1.joblib is small + reproducible; confirm not caught by a data/ ignore rule)
git commit -m "feat(scoring): LOCO + k-fold validation notebook + calibrator_v1 artifact (#708)"
```

---

### Task 6: README outputs index + epic doc

**Files:**
- Modify: `research/experiments/issue_547_immunogenicity_calibration/README.md`

- [ ] **Step 1: Update the README** — under "Layout" replace the `(later) notebook.ipynb + outputs/` line with the real files (`calibrator.py`, `score_cohort.py`, `notebook.ipynb`, `tests/`, `outputs/calibrator_v1.joblib`). Add an "Outputs index" entry for the artifact + figures, and a "How to reproduce" block (env setup → `score_cohort.py` → execute notebook). Note the LOCO design + proxy caveat in one line, linking the spec.

- [ ] **Step 2: Commit**

```bash
git add research/experiments/issue_547_immunogenicity_calibration/README.md
git commit -m "docs(scoring): #708 calibrator outputs index + reproduce steps (#708)"
```

---

## Post-implementation (not code tasks)

- **Open PR** → flip PR Status to `Ready for review`; offer `@claude review`.
- **Lab-notebook entry** (`research/lab_notebook/scientist.md`) — written AFTER review, BEFORE merge: per-cohort calibration diagnostics, the LOCO-vs-within-cohort shift gap, adaptive-vs-fixed verdict, referencing the PR/#708.
- **Tick** #708 ACs + PR Test plan; **merge** via `scripts/audit_and_merge.sh <PR>`.
- **Board lifecycle** already at In progress; PR + Issue → In review on review request; Done auto on merge.

## Self-Review

- **Spec coverage:** centered-isotonic + R oracle (Task 1), adaptive+fixed KDE module on `genotype_presentation_score` (Task 2), serialized artifact (Tasks 3, 5), MHCflurry precompute/subsample/true-count prior (Task 4), LOCO + within-cohort contrast + per-cohort diagnostics (Task 5), output column name (Task 2 docstring + #709 consumes), README/reproduce (Task 6), lab notebook (post-impl). All ACs mapped.
- **Placeholders:** none — every code step has full code; the only conditional is the documented R-unavailable fixture fallback in Task 1 Step 2.
- **Type consistency:** `centered_isotonic(x, y_iso, w)->(cx,cy)` used identically in Tasks 1/2; `PresentationCalibrator(kde_mode, n_grid)` + `fit(scores, labels, n_pos_true, n_neg_true, fit_cohorts)` + `transform` + `save`/`load` consistent across Tasks 2/3/5; `genotype_score` / `normalize_hla` / `build_scored_cohort` consistent in Task 4 + tests.
- **Known soft spot:** `score_cohort.py`'s MHCflurry call pattern (per-allele scoring loop) should be verified against the installed MHCflurry API on first run (Task 4 Step 6); the genotype combine + HLA-normalize helpers are API-independent and unit-tested separately.
