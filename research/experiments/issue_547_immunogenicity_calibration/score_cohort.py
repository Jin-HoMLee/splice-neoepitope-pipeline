"""Cohort scoring module: genotype_score, normalize_hla, and build_scored_cohort.

`genotype_score` and `normalize_hla` are pure-Python (no mhcflurry dependency)
and can be imported freely from any Python environment.

`build_scored_cohort` imports mhcflurry lazily — it is NOT called at module
import time. This lets the formula/HLA unit tests run under research/.venv
(Python 3.14, no mhcflurry).
"""
import logging
import math
import re
from pathlib import Path

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# HLA locus weights (A/B = full, C = half)
# ---------------------------------------------------------------------------
_LOCUS_WEIGHT = {"A": 1.0, "B": 1.0, "C": 0.5}
_DEFAULT_WEIGHT = 1.0   # fallback for any unlisted locus


def normalize_hla(allele: str) -> str:
    """Normalise an HLA allele string to the canonical MHCflurry form ``HLA-X*GG:FF``.

    Handles four incoming variants:
    - ``A*01:01``       → ``HLA-A*01:01``   (bare gene*group:protein)
    - ``HLA-A01:01``    → ``HLA-A*01:01``   (HLA-prefix but missing star)
    - ``HLA-A*01:01``   → ``HLA-A*01:01``   (already canonical, no-op)
    - ``C*07:02``       → ``HLA-C*07:02``   (C locus, bare)
    """
    s = allele.strip()
    # Already canonical
    if re.match(r"^HLA-[A-Z]\*\d{2}:\d{2}", s):
        return s
    # HLA-A01:01 form — has HLA- prefix but no star
    m = re.match(r"^HLA-([A-Z])(\d{2}:\d{2}.*)$", s)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}"
    # A*01:01 or C*07:02 — bare gene*group:protein
    m = re.match(r"^([A-Z])\*(\d{2}:\d{2}.*)$", s)
    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}"
    # Unrecognised — return as-is (let MHCflurry raise on its own)
    return s


def _locus_weight(canonical_allele: str) -> float:
    """Return the per-locus HLA weight for a canonical allele string."""
    m = re.match(r"^HLA-([A-Z])", canonical_allele)
    if m:
        return _LOCUS_WEIGHT.get(m.group(1), _DEFAULT_WEIGHT)
    return _DEFAULT_WEIGHT


def genotype_score(per_allele: dict) -> float:
    """Compute the genotype presentation score from per-allele MHCflurry scores.

    Formula: log(1 + sum(w_i * s_i))

    where w_i is the locus weight (A=1.0, B=1.0, C=0.5) and s_i is the
    presentation_score for allele i.

    Parameters
    ----------
    per_allele : dict[str, float]
        Mapping of canonical allele string → MHCflurry presentation_score.
        Keys should be in ``HLA-X*GG:FF`` form (as returned by normalize_hla).

    Returns
    -------
    float
        genotype_presentation_score (non-negative; 0.0 for empty dict).
    """
    if not per_allele:
        return 0.0
    weighted_sum = sum(_locus_weight(a) * s for a, s in per_allele.items())
    return math.log1p(weighted_sum)


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def _load_neoranking(data_dir: Path):
    """Load NeoRanking cohort TSV; return DataFrame with standardised columns.

    Columns returned: peptide, allele (canonical), label (1=immunogenic), cohort.
    Only rows with response_type in {'CD8', 'negative'} are included;
    'not_tested' rows are excluded.
    """
    import pandas as pd

    path = data_dir / "neoranking" / "Neopep_data_org.txt"
    df = pd.read_csv(path, sep="\t", low_memory=False,
                     usecols=["patient", "dataset", "response_type",
                               "mutant_seq", "mutant_best_alleles"])
    # Keep only labelled rows
    df = df[df["response_type"].isin(["CD8", "negative"])].copy()
    # Use the best allele for each peptide (first allele in comma-sep list)
    df["allele"] = df["mutant_best_alleles"].str.split(",").str[0].str.strip()
    df["allele"] = df["allele"].apply(normalize_hla)
    df["label"] = (df["response_type"] == "CD8").astype(int)
    df = df.rename(columns={"mutant_seq": "peptide", "dataset": "cohort"})
    df["patient"] = df["patient"].astype(str)
    return df[["peptide", "allele", "patient", "label", "cohort"]].reset_index(drop=True)


def _load_improve(data_dir: Path):
    """Load IMPROVE cohort TSV; return DataFrame with standardised columns.

    Columns returned: peptide, allele (canonical), label, cohort.
    IMPROVE has no patient column — assign cohort='IMPROVE'.
    """
    import pandas as pd

    path = data_dir / "improve_borch" / "In_house_neoepitope_for_CV.tsv"
    df = pd.read_csv(path, sep="\t")
    df["allele"] = df["HLA_allele"].apply(normalize_hla)
    df["label"] = df["response"].astype(int)
    df["cohort"] = "IMPROVE"
    df["patient"] = None
    df = df.rename(columns={"Mut_peptide": "peptide"})
    return df[["peptide", "allele", "patient", "label", "cohort"]].reset_index(drop=True)


def _hla_genotypes(data_dir: Path) -> dict:
    """Load NeoRanking HLA_allotypes.txt; return {patient_str: [canonical_allele, ...]}.

    Only NeoRanking has per-patient genotypes; IMPROVE records are single-allele.
    """
    import pandas as pd

    path = data_dir / "neoranking" / "HLA_allotypes.txt"
    df = pd.read_csv(path, sep="\t")
    result = {}
    for _, row in df.iterrows():
        patient = str(row["Patient"]).strip()
        alleles = [normalize_hla(a.strip()) for a in str(row["Alleles"]).split(",")]
        result[patient] = alleles
    return result


# ---------------------------------------------------------------------------
# Main scoring entry point
# ---------------------------------------------------------------------------

def build_scored_cohort(
    data_dir: Path,
    output_parquet: Path,
    n_neg: int = 50_000,
    seed: int = 42,
) -> None:
    """Score the combined NeoRanking + IMPROVE cohort with MHCflurry and write parquet.

    For each (peptide, allele) pair:
      1. Run Class1PresentationPredictor.predict(peptides=[peptide], alleles={a:[a]})
         for each allele in the patient's genotype.
      2. Extract presentation_score per allele → per_allele dict.
      3. Compute genotype_score(per_allele).

    Subsampling strategy:
      - ALL positives (label=1) are kept.
      - n_neg negatives (label=0) are sampled cohort-proportionally with fixed seed.
      - True per-cohort counts (before subsampling) are written to a companion
        ``<output_parquet>.true_counts.csv``.

    The function is idempotent: if output_parquet already exists it is skipped
    (use skip-if-cached semantics for re-runs).

    Parameters
    ----------
    data_dir : Path
        Directory containing ``neoranking/`` and ``improve_borch/`` subdirs.
    output_parquet : Path
        Destination for the scored cohort parquet file.
    n_neg : int
        Number of negatives to retain after subsampling.
    seed : int
        Random seed for reproducible subsampling.
    """
    # Lazy mhcflurry import — only pulled in when actually scoring
    from mhcflurry import Class1PresentationPredictor  # noqa: PLC0415

    import pandas as pd
    import numpy as np

    output_parquet = Path(output_parquet)
    if output_parquet.exists():
        logger.info("skip-if-cached: %s already exists", output_parquet)
        return

    output_parquet.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load cohorts
    # ------------------------------------------------------------------
    logger.info("Loading NeoRanking cohort ...")
    neo = _load_neoranking(data_dir)
    logger.info("  NeoRanking rows: %d", len(neo))

    logger.info("Loading IMPROVE cohort ...")
    imp = _load_improve(data_dir)
    logger.info("  IMPROVE rows: %d", len(imp))

    df = pd.concat([neo, imp], ignore_index=True)
    logger.info("Combined rows: %d", len(df))

    # ------------------------------------------------------------------
    # 2. Record TRUE per-cohort label counts (BEFORE subsampling)
    # ------------------------------------------------------------------
    true_counts = (
        df.groupby(["cohort", "label"])
        .size()
        .reset_index(name="count")
    )
    counts_path = str(output_parquet) + ".true_counts.csv"
    true_counts.to_csv(counts_path, index=False)
    logger.info("True counts written to %s", counts_path)

    # ------------------------------------------------------------------
    # 3. Subsample: keep ALL positives + n_neg cohort-stratified negatives
    # ------------------------------------------------------------------
    positives = df[df["label"] == 1].copy()
    negatives = df[df["label"] == 0].copy()

    neg_cohort_counts = negatives["cohort"].value_counts()
    total_neg = len(negatives)
    n_neg_actual = min(n_neg, total_neg)

    rng = np.random.default_rng(seed)
    sampled_neg_parts = []
    for cohort, cohort_neg in negatives.groupby("cohort"):
        frac = len(cohort_neg) / total_neg
        n_cohort = max(1, round(n_neg_actual * frac))
        n_cohort = min(n_cohort, len(cohort_neg))
        idx = rng.choice(len(cohort_neg), size=n_cohort, replace=False)
        sampled_neg_parts.append(cohort_neg.iloc[idx])

    sampled_neg = pd.concat(sampled_neg_parts, ignore_index=True)
    df_sub = pd.concat([positives, sampled_neg], ignore_index=True)
    logger.info(
        "After subsampling: %d rows (%d pos, %d neg)",
        len(df_sub), len(positives), len(sampled_neg),
    )

    # ------------------------------------------------------------------
    # 4. Load per-patient HLA genotypes (NeoRanking only)
    # ------------------------------------------------------------------
    genotypes = _hla_genotypes(data_dir)

    # ------------------------------------------------------------------
    # 5. Score each row
    # ------------------------------------------------------------------
    predictor = Class1PresentationPredictor.load()
    logger.info("MHCflurry predictor loaded.")

    presentation_scores = []
    genotype_scores = []

    for i, row in df_sub.iterrows():
        peptide = row["peptide"]
        # Determine alleles to score against
        if row["patient"] is not None and row["patient"] in genotypes:
            alleles = genotypes[row["patient"]]
        else:
            # IMPROVE / patients without genotype: use single allele from row
            alleles = [row["allele"]]

        # Per-allele scoring: one predict() call per allele
        per_allele = {}
        for a in alleles:
            try:
                result = predictor.predict(
                    peptides=[peptide],
                    alleles={a: [a]},
                    verbose=0,
                )
                score = float(result["presentation_score"].iloc[0])
            except Exception as exc:  # noqa: BLE001
                logger.warning("predict failed for peptide=%s allele=%s: %s", peptide, a, exc)
                score = 0.0
            per_allele[a] = score

        gs = genotype_score(per_allele)
        best_ps = max(per_allele.values()) if per_allele else 0.0
        presentation_scores.append(best_ps)
        genotype_scores.append(gs)

        if (i + 1) % 5000 == 0:
            logger.info("  scored %d / %d rows", i + 1, len(df_sub))

    df_sub = df_sub.copy()
    df_sub["presentation_score"] = presentation_scores
    df_sub["genotype_presentation_score"] = genotype_scores

    # ------------------------------------------------------------------
    # 6. Write parquet
    # ------------------------------------------------------------------
    df_sub.to_parquet(output_parquet, index=False)
    logger.info("Scored cohort written to %s", output_parquet)
