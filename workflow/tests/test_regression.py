#!/usr/bin/env python3
"""test_regression.py — Pipeline regression tests against a golden snapshot.

What is a regression test?
---------------------------
A regression test checks that a code change has not accidentally broken
something that was working before.  Here, "broken" means producing different
pipeline outputs (different junction counts, different peptide sets, different
strong-binder predictions) from the same input data.

How it works
------------
After every successful test-dataset run, key output metrics are stored in a
committed JSON file called the *golden snapshot* (golden/test_metrics.json).
Whenever the pipeline code changes, this test re-runs the test dataset and
compares the new outputs against the golden snapshot.  Any difference causes
the test to fail, which forces the developer to make a deliberate decision:

  * Unexpected difference  → bug introduced; fix the code.
  * Expected difference    → update the golden snapshot with --update-golden
                             and commit it so the PR diff shows the intended
                             change in outputs clearly.

What is compared
----------------
  junctions   total unannotated junctions and breakdown by origin
              (tumor_exclusive vs normal_shared)
  contigs     number of assembled 50 nt contigs
  peptides    total and unique junction-spanning 9-mers
  predictions binder-class counts (strong / weak / non) and the exact set of
              strong-binder peptide sequences

Usage
-----
Run as a pytest (CI / pre-push check):
    pytest workflow/tests/test_regression.py

Manually compare current results to the golden snapshot:
    python workflow/tests/test_regression.py

Update the golden snapshot after an intentional pipeline change:
    python workflow/tests/test_regression.py --update-golden

The golden file must be committed to git so that output changes are visible
as a diff in PR reviews alongside the code changes that caused them.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import yaml

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[2]
GOLDEN_DIR = Path(__file__).parent / "golden"

# Read the first patient_id from the test samplesheet — the same source of
# truth the pipeline uses, so this stays in sync automatically.
def _first_patient_id() -> str:
    tsv = REPO_ROOT / "config" / "samples" / "patient_001_test.tsv"
    with tsv.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pid = (row.get("patient_id") or "").strip()
            if not pid or pid.startswith("#"):
                continue
            return pid
    raise RuntimeError(f"No patient_id found in {tsv}")

PATIENT_ID = _first_patient_id()
RESULTS = REPO_ROOT / "results" / PATIENT_ID
GOLDEN_FILE = GOLDEN_DIR / "test_metrics.json"


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override into base (same semantics as Snakemake configfile merge)."""
    result = dict(base)
    for key, val in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(val, dict):
            result[key] = _deep_merge(result[key], val)
        else:
            result[key] = val
    return result


def _load_merged_config() -> dict:
    """Load config.yaml merged with test_config.yaml, mimicking Snakemake's merge order."""
    base = yaml.safe_load((REPO_ROOT / "config" / "config.yaml").read_text())
    test = yaml.safe_load((REPO_ROOT / "config" / "test_config.yaml").read_text())
    return _deep_merge(base, test)


def _read_run_config() -> dict:
    """Collect the key input parameters that produced the current test outputs."""
    cfg = _load_merged_config()
    samples_tsv = REPO_ROOT / cfg.get("samples_tsv", "config/test_samples.tsv")
    samples = [
        {"sample_id": r["sample_id"], "sample_type": r["sample_type"]}
        for r in csv.DictReader(samples_tsv.open(), delimiter="\t")
    ]
    ref = cfg.get("reference", {})
    asm = cfg.get("assembly", {})
    mhc = cfg.get("mhcflurry", {})
    hla = cfg.get("hla", {})

    return {
        "patient_id": PATIENT_ID,
        "samples": samples,
        "reference_genome": ref.get("genome", "unknown"),
        "gencode_gtf": Path(ref.get("gencode_gtf", "")).name,
        "hla_typing_enabled": hla.get("enabled", False),
        "upstream_nt": asm.get("upstream_nt", 26),
        "downstream_nt": asm.get("downstream_nt", 24),
        "ic50_strong_nM": mhc.get("ic50_strong", 50),
        "ic50_weak_nM": mhc.get("ic50_weak", 500),
    }


# ---------------------------------------------------------------------------
# Metric extraction
# ---------------------------------------------------------------------------

def _read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def collect_metrics() -> dict:
    """Read current test-run outputs and return a metrics dict.

    Each key maps to a sub-dict of counts derived from one pipeline output
    file.  Only counts and sorted peptide lists are stored — not raw rows —
    so the snapshot stays small and human-readable.
    """
    # Junctions ---------------------------------------------------------------
    junc_rows = _read_tsv(RESULTS / "junctions" / "novel_junctions.tsv")
    junc_origins: dict[str, int] = {}
    for r in junc_rows:
        origin = r["junction_origin"]
        junc_origins[origin] = junc_origins.get(origin, 0) + 1

    # Contigs -----------------------------------------------------------------
    contigs_fa = RESULTS / "contigs" / "contigs.fa"
    n_contigs = sum(1 for line in contigs_fa.open() if line.startswith(">"))

    # Peptides ----------------------------------------------------------------
    pep_rows = _read_tsv(RESULTS / "peptides" / "peptides.tsv")
    n_unique_peptides = len(set(r["peptide"] for r in pep_rows))

    # Predictions -------------------------------------------------------------
    pred_rows = _read_tsv(RESULTS / "predictions" / "mhc_affinity.tsv")
    binder_counts: dict[str, int] = {}
    for r in pred_rows:
        bc = r["binder_class"]
        binder_counts[bc] = binder_counts.get(bc, 0) + 1
    # Store the exact strong-binder peptide set so we catch substitutions
    # (same count, different peptides) as well as count changes.
    strong_peptides = sorted(set(r["peptide"] for r in pred_rows if r["binder_class"] == "strong"))

    return {
        "run_config": _read_run_config(),
        "junctions": {
            "total_unannotated": len(junc_rows),
            "by_origin": junc_origins,
        },
        "contigs": {
            "n_contigs": n_contigs,
        },
        "peptides": {
            "n_total": len(pep_rows),
            "n_unique": n_unique_peptides,
        },
        "predictions": {
            "by_class": binder_counts,
            "strong_peptides": strong_peptides,
        },
    }


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def compare(current: dict, golden: dict) -> list[str]:
    """Return a list of human-readable difference strings (empty = identical).

    Each string describes one metric that changed, showing the golden value
    and the current value side-by-side so it is clear what shifted.
    """
    diffs: list[str] = []

    def _check(label: str, cur, gold):
        if cur != gold:
            diffs.append(f"  {label}:\n    golden:  {gold!r}\n    current: {cur!r}")

    # Junctions
    _check("junctions.total_unannotated",
           current["junctions"]["total_unannotated"],
           golden["junctions"]["total_unannotated"])
    for origin in sorted(set(current["junctions"]["by_origin"]) | set(golden["junctions"]["by_origin"])):
        _check(f"junctions.by_origin[{origin}]",
               current["junctions"]["by_origin"].get(origin, 0),
               golden["junctions"]["by_origin"].get(origin, 0))

    # Contigs
    _check("contigs.n_contigs",
           current["contigs"]["n_contigs"],
           golden["contigs"]["n_contigs"])

    # Peptides
    _check("peptides.n_total",  current["peptides"]["n_total"],  golden["peptides"]["n_total"])
    _check("peptides.n_unique", current["peptides"]["n_unique"], golden["peptides"]["n_unique"])

    # Predictions — binder-class counts
    for cls in sorted(set(current["predictions"]["by_class"]) | set(golden["predictions"]["by_class"])):
        _check(f"predictions.by_class[{cls}]",
               current["predictions"]["by_class"].get(cls, 0),
               golden["predictions"]["by_class"].get(cls, 0))

    # Predictions — exact strong-binder peptide set
    # A count match alone is not sufficient: the same number of strong binders
    # could contain completely different peptide sequences.
    cur_strong = set(current["predictions"]["strong_peptides"])
    gold_strong = set(golden["predictions"]["strong_peptides"])
    lost   = sorted(gold_strong - cur_strong)
    gained = sorted(cur_strong - gold_strong)
    if lost:
        diffs.append(f"  predictions.strong_peptides — lost {len(lost)}: {lost}")
    if gained:
        diffs.append(f"  predictions.strong_peptides — gained {len(gained)}: {gained}")

    return diffs


# ---------------------------------------------------------------------------
# Pytest entry point
# ---------------------------------------------------------------------------

def test_results_match_golden():
    """Pytest: current pipeline outputs must match the committed golden snapshot.

    This test requires local pipeline outputs (results/test/) and is skipped
    automatically in CI where those files are not present.  Run it locally
    after executing the test pipeline:

        snakemake --cores 4 --use-conda --configfile config/test_config.yaml
        pytest workflow/tests/test_regression.py

    If this test fails after a deliberate pipeline change, update the snapshot:
        python workflow/tests/test_regression.py --update-golden
    then commit golden/test_metrics.json together with the code change.
    """
    import pytest

    sentinel = RESULTS / "junctions" / "novel_junctions.tsv"
    if not sentinel.exists():
        pytest.skip("Pipeline results not found — run the test pipeline first")

    assert GOLDEN_FILE.exists(), (
        f"Golden snapshot not found: {GOLDEN_FILE}\n"
        "Run the test pipeline first, then:\n"
        "  python workflow/tests/test_regression.py --update-golden"
    )
    golden = json.loads(GOLDEN_FILE.read_text())
    current = collect_metrics()
    diffs = compare(current, golden)
    assert not diffs, (
        "Pipeline outputs differ from the golden snapshot:\n"
        + "\n".join(diffs)
        + "\n\nIf this change is intentional, update the snapshot:\n"
        + "  python workflow/tests/test_regression.py --update-golden\n"
        + "then commit golden/test_metrics.json alongside your code change."
    )


# ---------------------------------------------------------------------------
# CLI: manual comparison and golden update
# ---------------------------------------------------------------------------

def update_golden() -> None:
    """Overwrite the golden snapshot with the current test-run outputs."""
    metrics = collect_metrics()
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    GOLDEN_FILE.write_text(json.dumps(metrics, indent=2) + "\n")
    print(f"Golden snapshot written to {GOLDEN_FILE}")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare current pipeline outputs against the golden snapshot.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python workflow/tests/test_regression.py             # compare\n"
            "  python workflow/tests/test_regression.py --update-golden  # refresh\n"
        ),
    )
    parser.add_argument(
        "--update-golden", action="store_true",
        help="Overwrite the golden snapshot with the current test-run outputs.",
    )
    args = parser.parse_args()

    if args.update_golden:
        update_golden()
    else:
        if not GOLDEN_FILE.exists():
            print(f"No golden snapshot at {GOLDEN_FILE}")
            print("Run the test pipeline, then: python workflow/tests/test_regression.py --update-golden")
            sys.exit(1)
        golden = json.loads(GOLDEN_FILE.read_text())
        current = collect_metrics()
        diffs = compare(current, golden)
        if diffs:
            print("DIFFERENCES from golden snapshot:")
            print("\n".join(diffs))
            sys.exit(1)
        else:
            print("OK — outputs match the golden snapshot.")
