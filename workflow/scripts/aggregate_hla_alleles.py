#!/usr/bin/env python3
"""aggregate_hla_alleles.py — Aggregate per-sample OptiType HLA calls into a
per-patient alleles TSV.

Normal-first policy
-------------------
HLA alleles are germline, so normal tissue gives a cleaner signal than tumor
(which may carry CNV or LOH events at the HLA locus). For each locus (A, B, C)
the aggregator:

  1. Collects calls from Solid Tissue Normal / Blood Derived Normal samples
     that pass the min_reads threshold.
  2. Falls back to the Primary Tumor sample if no normal call is available.
  3. Falls back to config-defined fallback alleles if no sample call is
     available or confident.

QC output
---------
A side TSV records per-locus:
  - source      (normal / tumor / fallback)
  - reads       (total HLA reads from OptiType for the chosen sample)
  - discrepancy (if normal and tumor calls differ)

Input (via Snakemake)
---------------------
  snakemake.input.result_tsvs  — list of {sample}_result.tsv paths from OptiType
  snakemake.output.alleles     — output alleles TSV path
  snakemake.output.qc          — output QC TSV path
  snakemake.params.samples_tsv — path to the samples TSV
  snakemake.params.min_reads   — minimum total HLA reads for a confident call
  snakemake.params.fallback    — dict of fallback alleles per locus {A: ..., B: ..., C: ...}
"""

import csv
import logging
from pathlib import Path


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

_LOCI = ["A", "B", "C"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_samples(tsv_path: str) -> dict[str, str]:
    """Return {sample_id: sample_type} from the samples TSV."""
    result: dict[str, str] = {}
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sid = (row.get("sample_id") or "").strip()
            if not sid or sid.startswith("#"):
                continue
            result[sid] = (row.get("sample_type") or "Unknown").strip()
    return result


def load_optitype_result(tsv_path: str) -> dict[str, tuple[list[str], int]]:
    """Load an OptiType *_result.tsv file.

    OptiType writes one result row per solution. We take the first row.
    Alleles are in short format (e.g. ``A*02:01``); we prepend ``HLA-``.

    Returns a dict mapping locus letter to (allele_list, reads):
        {'A': (['HLA-A*02:01', 'HLA-A*24:02'], 1234), 'B': ..., 'C': ...}

    Returns an empty dict when the file is missing, empty, or has no data row.
    """
    path = Path(tsv_path)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    try:
        with path.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                # OptiType TSV has an unnamed index column; skip blank rows
                if not any(row.get(k, "").strip() for k in ("A1", "B1", "C1")):
                    continue
                reads = int(float(row.get("Reads", 0) or 0))
                result: dict[str, tuple[list[str], int]] = {}
                for locus in _LOCI:
                    a1 = normalise_allele(row.get(f"{locus}1", ""))
                    a2 = normalise_allele(row.get(f"{locus}2", ""))
                    if a1:
                        result[locus] = ([a1, a2] if a2 else [a1], reads)
                return result
    except Exception as exc:
        log.warning("Could not parse OptiType TSV %s: %s", tsv_path, exc)
    return {}


def normalise_allele(raw: str) -> str:
    """Return allele in HLA-X*YY:ZZ format (4-digit resolution).

    OptiType returns alleles like ``A*02:01``; MHCflurry expects ``HLA-A*02:01``.
    Truncates to the first two colon-separated fields.
    """
    raw = raw.strip()
    if not raw:
        return ""
    if not raw.startswith("HLA-"):
        raw = f"HLA-{raw}"
    parts = raw.split(":")
    return ":".join(parts[:2])


# ---------------------------------------------------------------------------
# Main aggregation
# ---------------------------------------------------------------------------

def aggregate(
    result_tsvs: list[str],
    samples_tsv: str,
    min_reads: int,
    fallback: dict[str, str],
    out_alleles: str,
    out_qc: str,
) -> None:
    """Aggregate OptiType results into per-patient alleles.tsv and hla_qc.tsv."""
    sample_types = read_samples(samples_tsv)
    log.info("Loaded %d samples from %s", len(sample_types), samples_tsv)

    # Map sample_id → result TSV path
    tsv_by_sample: dict[str, str] = {}
    for tpath in result_tsvs:
        # Strip the _result.tsv suffix to get sample_id
        sample_id = Path(tpath).name.replace("_result.tsv", "")
        tsv_by_sample[sample_id] = tpath

    normal_samples = [
        s for s in tsv_by_sample
        if "normal" in sample_types.get(s, "").lower()
    ]
    tumor_samples = [
        s for s in tsv_by_sample
        if "tumor" in sample_types.get(s, "").lower()
    ]
    log.info("Normal samples: %s", normal_samples)
    log.info("Tumor samples:  %s", tumor_samples)

    allele_rows: list[dict] = []
    qc_rows: list[dict] = []

    for locus in _LOCI:
        def _get_calls(sample_list: list[str]) -> dict[str, tuple[list[str], int]]:
            calls: dict[str, tuple[list[str], int]] = {}
            for sid in sample_list:
                result = load_optitype_result(tsv_by_sample[sid])
                if locus not in result:
                    continue
                alleles, reads = result[locus]
                if reads >= min_reads:
                    calls[sid] = (alleles, reads)
                else:
                    log.warning(
                        "Sample %s locus %s: only %d total HLA reads (min %d), skipping",
                        sid, locus, reads, min_reads,
                    )
            return calls

        normal_calls = _get_calls(normal_samples)
        tumor_calls  = _get_calls(tumor_samples)

        # Normal-first policy
        if normal_calls:
            best = max(normal_calls, key=lambda s: normal_calls[s][1])
            chosen_alleles, reads = normal_calls[best]
            source = "normal"
            log.info("Locus %s: using normal sample %s (%d reads)", locus, best, reads)
        elif tumor_calls:
            best = max(tumor_calls, key=lambda s: tumor_calls[s][1])
            chosen_alleles, reads = tumor_calls[best]
            source = "tumor"
            log.warning("Locus %s: no normal call, using tumor sample %s", locus, best)
        else:
            fb = fallback.get(locus, "")
            chosen_alleles = [fb] if fb else []
            reads = 0
            source = "fallback"
            log.warning("Locus %s: no confident call, using fallback '%s'", locus, fb)

        allele1 = chosen_alleles[0] if len(chosen_alleles) > 0 else fallback.get(locus, "")
        allele2 = chosen_alleles[1] if len(chosen_alleles) > 1 else allele1

        allele_rows.append({"locus": locus, "allele1": allele1, "allele2": allele2})

        # QC: normal/tumor concordance check
        discrepancy = ""
        if normal_calls and tumor_calls:
            n_best = max(normal_calls, key=lambda s: normal_calls[s][1])
            t_best = max(tumor_calls,  key=lambda s: tumor_calls[s][1])
            n_set = set(normal_calls[n_best][0])
            t_set = set(tumor_calls[t_best][0])
            if n_set != t_set:
                discrepancy = f"normal={sorted(n_set)} tumor={sorted(t_set)}"
                log.warning("Locus %s: normal/tumor discrepancy — %s", locus, discrepancy)

        qc_rows.append({
            "locus": locus, "allele1": allele1, "allele2": allele2,
            "source": source, "reads": reads, "discrepancy": discrepancy,
        })

    Path(out_alleles).parent.mkdir(parents=True, exist_ok=True)
    with open(out_alleles, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["locus", "allele1", "allele2"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(allele_rows)
    log.info("Wrote alleles to %s", out_alleles)

    with open(out_qc, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["locus", "allele1", "allele2", "source", "reads", "discrepancy"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(qc_rows)
    log.info("Wrote QC to %s", out_qc)


# ---------------------------------------------------------------------------
# Snakemake entry point
# ---------------------------------------------------------------------------

def _snakemake_main() -> None:
    logging.getLogger().addHandler(logging.FileHandler(snakemake.log[0]))  # noqa: F821  # type: ignore[name-defined]
    aggregate(
        result_tsvs=list(snakemake.input.result_tsvs),  # noqa: F821
        samples_tsv=snakemake.params.samples_tsv,
        min_reads=snakemake.params.min_reads,
        fallback=snakemake.params.fallback,
        out_alleles=snakemake.output.alleles,
        out_qc=snakemake.output.qc,
    )


try:
    snakemake  # noqa: F821  # type: ignore[name-defined]
    _snakemake_main()
except NameError:
    pass  # imported directly (e.g. for testing)
