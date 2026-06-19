"""Cohort scoring module: genotype_score, normalize_hla, and build_scored_cohort.

`genotype_score` and `normalize_hla` are pure-Python (no mhcflurry dependency)
and can be imported freely from any Python environment.

`build_scored_cohort` imports mhcflurry lazily — it is NOT called at module
import time. This lets the formula/HLA unit tests run under research/.venv
(Python 3.14, no mhcflurry).
"""
import logging
import re
from pathlib import Path

logger = logging.getLogger(__name__)

def normalize_hla(allele: str) -> str:
    """Normalise an HLA allele string to the canonical MHCflurry form ``HLA-X*GG:FF``.

    Handles four incoming variants:
    - ``A*01:01``       → ``HLA-A*01:01``   (bare gene*group:protein)
    - ``HLA-A01:01``    → ``HLA-A*01:01``   (HLA-prefix but missing star)
    - ``HLA-A*01:01``   → ``HLA-A*01:01``   (already canonical, no-op)
    - ``C*07:02``       → ``HLA-C*07:02``   (C locus, bare)

    Only single-character Class I gene names (A, B, C) are handled. Multi-character loci
    (e.g. Class II ``DRB1``) match none of the patterns and are returned **unchanged** —
    safe for the current all-Class-I cohorts, but would need extending for Class II reuse.
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


def genotype_score(per_allele: dict, hla_c_weight: float = 0.5) -> float:
    """genotype_presentation_score = 1 - prod(1 - w_i * p_i); w_i = hla_c_weight for HLA-C else 1.0.
    Mirrors workflow/scripts/run_mhcflurry.py.

    Parameters
    ----------
    per_allele : dict[str, float]
        Mapping of canonical allele string → MHCflurry presentation_score.
        Keys should be in ``HLA-X*GG:FF`` form (as returned by normalize_hla).
    hla_c_weight : float
        Weight applied to HLA-C alleles (default 0.5). All other loci use 1.0.

    Returns
    -------
    float
        genotype_presentation_score in [0, 1]; 0.0 for empty dict.
    """
    if not per_allele:
        return 0.0
    product = 1.0
    for allele, p in per_allele.items():
        w = hla_c_weight if allele.startswith("HLA-C") else 1.0
        product *= (1.0 - w * float(p))
    # round to 6 dp for stable parquet diffs in this precompute; the production
    # run_mhcflurry.py this mirrors does NOT round (offset <1e-6, scores in [0, 1]).
    return round(1.0 - product, 6)


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
    # Best allele per peptide (first in the comma-sep list). NOTE: for NeoRanking this
    # `allele` column is a FALLBACK only — build_scored_cohort scores each row via the full
    # per-patient HLA genotype. The single-allele proxy is the genuine path only for IMPROVE
    # rows (no patient genotype), so this asymmetry is intentional.
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
    counts_path = output_parquet.parent / (output_parquet.name + ".true_counts.csv")
    true_counts.to_csv(counts_path, index=False)
    logger.info("True counts written to %s", counts_path)

    # ------------------------------------------------------------------
    # 3. Subsample: keep ALL positives + n_neg cohort-stratified negatives
    # ------------------------------------------------------------------
    positives = df[df["label"] == 1].copy()
    negatives = df[df["label"] == 0].copy()

    total_neg = len(negatives)
    n_neg_actual = min(n_neg, total_neg)

    rng = np.random.default_rng(seed)
    sampled_neg_parts = []
    # Per-cohort round() can make the parts sum to n_neg_actual ± a couple (classic
    # integer-rounding drift); harmless at this scale — don't add a strict sum==n_neg assert.
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
    # 5. Score each row — batched: one predict() call per unique allele
    # ------------------------------------------------------------------
    predictor = Class1PresentationPredictor.load()
    logger.info("MHCflurry predictor loaded.")

    # Single pass: build per-row allele lists AND the allele→peptide map together.
    # pd.notna guards against NaN patient IDs (IMPROVE rows set patient=None,
    # which becomes NaN after pd.concat — "is not None" misses that case).
    from collections import defaultdict
    row_alleles = []
    allele_peptides: dict = defaultdict(set)
    for _, row in df_sub.iterrows():
        if pd.notna(row["patient"]) and row["patient"] in genotypes:
            alleles = genotypes[row["patient"]]
        else:
            alleles = [row["allele"]]
        row_alleles.append(alleles)
        for a in alleles:
            allele_peptides[a].add(row["peptide"])

    # One predict() call per unique allele
    allele_peptide_score: dict = {}  # (allele, peptide) -> score
    unique_alleles = sorted(allele_peptides)
    logger.info("Batching MHCflurry: %d unique alleles", len(unique_alleles))
    for idx, allele in enumerate(unique_alleles, 1):
        peptides_for_allele = sorted(allele_peptides[allele])
        try:
            result = predictor.predict(
                peptides=peptides_for_allele,
                alleles={allele: [allele]},
                verbose=0,
            )
            scores = dict(zip(result["peptide"], result["presentation_score"]))
        except Exception as exc:  # noqa: BLE001
            logger.warning("predict failed for allele=%s: %s", allele, exc)
            scores = {p: 0.0 for p in peptides_for_allele}
        for peptide, score in scores.items():
            allele_peptide_score[(allele, peptide)] = float(score)
        if idx % 50 == 0 or idx == len(unique_alleles):
            logger.info("  scored allele %d / %d", idx, len(unique_alleles))

    # Assemble per-row results from the lookup. NOTE: row_alleles[row_idx] was built in
    # this same df_sub.iterrows() order above — the two loops are coupled by position, so
    # df_sub must not be reordered/mutated between them.
    presentation_scores = []
    genotype_scores = []
    for row_idx, (_, row) in enumerate(df_sub.iterrows()):
        alleles = row_alleles[row_idx]
        peptide = row["peptide"]
        per_allele = {
            a: allele_peptide_score.get((a, peptide), 0.0)
            for a in alleles
        }
        gs = genotype_score(per_allele)
        best_ps = max(per_allele.values()) if per_allele else 0.0
        presentation_scores.append(best_ps)
        genotype_scores.append(gs)

        if (row_idx + 1) % 5000 == 0:
            logger.info("  assembled %d / %d rows", row_idx + 1, len(df_sub))

    df_sub = df_sub.copy()
    df_sub["presentation_score"] = presentation_scores
    df_sub["genotype_presentation_score"] = genotype_scores

    # ------------------------------------------------------------------
    # 6. Write parquet
    # ------------------------------------------------------------------
    df_sub.to_parquet(output_parquet, index=False)
    logger.info("Scored cohort written to %s", output_parquet)


if __name__ == "__main__":
    # Reproducible entry point: scores the combined cohort and writes the parquet
    # (+ companion .true_counts.csv) under outputs/. Run with the mhcflurry env:
    #   conda activate mhcflurry-scoring && python score_cohort.py
    # (Do NOT use `conda run` — it buffers stdout and hides the per-allele progress
    # logs during a 30+ min run; see CLAUDE.md.) Paths are resolved relative to this
    # file so the command works from any cwd.
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    _here = Path(__file__).resolve().parent
    build_scored_cohort(
        data_dir=_here / "data",
        output_parquet=_here / "outputs" / "scored_cohort_subsample.parquet",
    )
