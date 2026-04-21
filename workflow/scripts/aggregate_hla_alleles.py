#!/usr/bin/env python3
"""aggregate_hla_alleles.py — Aggregate per-sample OptiType HLA calls into a
per-patient alleles TSV.

Priority order
--------------
HLA alleles are germline, so reliable sources are preferred over noisier ones.
For each locus (A, B, C) the aggregator applies this cascade:

  1. OptiType from Primary Tumor — reflects what the tumor cell actually
     presents, including any HLA loss-of-heterozygosity (LOH) events.
  2. Clinical serology (Red Cross / WGS) — germline ground truth, more
     accurate than RNA-seq-based typing. Read from serology_X1/X2 columns
     in the samples TSV (germline/normal row only; tumor row should be empty).
  3. OptiType from Solid Tissue Normal / Blood Derived Normal — used when no
     serology is available.
  4. Config-defined fallback alleles when no sample call is available or
     confident (reads < min_reads_per_locus).

Null alleles (alleles ending in 'N', e.g. A*01:11N) are not expressed at the
cell surface and are excluded from the prediction allele list, but are stored
in the QC output for transparency.

QC output
---------
A side TSV records per-locus:
  - source              (tumor / serology / normal / fallback)
  - reads               (total HLA reads from OptiType for the chosen sample)
  - discrepancy         (if normal and tumor OptiType calls differ)
  - serology_allele1/2  (known alleles from serology, if provided)
  - serology_validation (match / match (null allele excluded) / mismatch: ...)

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


def load_serology(tsv_path: str) -> dict[str, tuple[str, str]]:
    """Read clinical serology alleles from the samples TSV.

    Returns {locus: (allele1, allele2)} using the first non-empty serology
    values found across all rows. Serology is per-patient; only the germline
    (normal/blood) row is expected to carry values.

    Null alleles (e.g. A*01:11N) are preserved as-is for QC display; callers
    filter them before using alleles for prediction.
    """
    result: dict[str, tuple[str, str]] = {}
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sid = (row.get("sample_id") or "").strip()
            if not sid or sid.startswith("#"):
                continue
            for locus in _LOCI:
                if locus in result:
                    continue
                allele1 = normalise_allele(row.get(f"serology_{locus}1", "") or "")
                allele2 = normalise_allele(row.get(f"serology_{locus}2", "") or "")
                if allele1:
                    result[locus] = (allele1, allele2 if allele2 else allele1)
    if result:
        log.info("Loaded serology for loci: %s", sorted(result))
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
                    allele1 = normalise_allele(row.get(f"{locus}1", ""))
                    allele2 = normalise_allele(row.get(f"{locus}2", ""))
                    if allele1:
                        result[locus] = ([allele1, allele2] if allele2 else [allele1], reads)
                return result
    except Exception as exc:
        log.warning("Could not parse OptiType TSV %s: %s", tsv_path, exc)
    return {}


def normalise_allele(raw: str) -> str:
    """Return allele in HLA-X*YY:ZZ format (4-digit resolution).

    OptiType returns alleles like ``A*02:01``; MHCflurry expects ``HLA-A*02:01``.
    Truncates to the first two colon-separated fields (preserves trailing N for
    null alleles, e.g. ``HLA-A*01:11N``).
    """
    raw = raw.strip()
    if not raw:
        return ""
    if not raw.startswith("HLA-"):
        raw = f"HLA-{raw}"
    parts = raw.split(":")
    return ":".join(parts[:2])


def _validate_vs_serology(
    optitype_alleles: list[str] | None,
    serology: tuple[str, str],
) -> str:
    """Compare an OptiType call against known serology alleles.

    Null alleles (ending in 'N') are not expressed and will not appear in
    OptiType calls; they are excluded from the serology set before comparison.

    Returns one of: 'match', 'match (null allele excluded)', 'mismatch: ...', ''.
    """
    if not optitype_alleles:
        return ""
    ser_allele1, ser_allele2 = serology
    ser_set = {ser_allele1, ser_allele2}
    ser_non_null = {allele for allele in ser_set if not allele.endswith("N")}
    optitype_set = set(optitype_alleles)

    if optitype_set == ser_set:
        return "match"
    if ser_non_null and optitype_set == ser_non_null:
        return "match (null allele excluded)"
    return f"mismatch: optitype={sorted(optitype_set)} serology={sorted(ser_set)}"


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
    serology = load_serology(samples_tsv)
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

        tumor_calls  = _get_calls(tumor_samples)
        normal_calls = _get_calls(normal_samples)

        # Priority: tumor OptiType > serology > normal OptiType > fallback
        if tumor_calls:
            best = max(tumor_calls, key=lambda s: tumor_calls[s][1])
            chosen_alleles, reads = tumor_calls[best]
            source = "tumor"
            log.info("Locus %s: using tumor sample %s (%d reads)", locus, best, reads)
        elif serology and locus in serology:
            ser_allele1, ser_allele2 = serology[locus]
            # Exclude null alleles from prediction (not expressed at cell surface)
            expressed = [a for a in [ser_allele1, ser_allele2] if a and not a.endswith("N")]
            chosen_alleles = expressed if expressed else [ser_allele1]
            reads = 0
            source = "serology"
            log.info("Locus %s: using serology (%s / %s)", locus, ser_allele1, ser_allele2)
        elif normal_calls:
            best = max(normal_calls, key=lambda s: normal_calls[s][1])
            chosen_alleles, reads = normal_calls[best]
            source = "normal"
            log.info("Locus %s: using normal sample %s (%d reads)", locus, best, reads)
        else:
            fb = fallback.get(locus, "")
            chosen_alleles = [fb] if fb else []
            reads = 0
            source = "fallback"
            log.warning("Locus %s: no confident call, using fallback '%s'", locus, fb)

        allele1 = chosen_alleles[0] if len(chosen_alleles) > 0 else fallback.get(locus, "")
        allele2 = chosen_alleles[1] if len(chosen_alleles) > 1 else allele1

        allele_rows.append({"locus": locus, "allele1": allele1, "allele2": allele2})

        # QC: normal/tumor OptiType concordance check
        discrepancy = ""
        if normal_calls and tumor_calls:
            n_best = max(normal_calls, key=lambda s: normal_calls[s][1])
            t_best = max(tumor_calls,  key=lambda s: tumor_calls[s][1])
            n_set = set(normal_calls[n_best][0])
            t_set = set(tumor_calls[t_best][0])
            if n_set != t_set:
                discrepancy = f"normal={sorted(n_set)} tumor={sorted(t_set)}"
                log.warning("Locus %s: normal/tumor discrepancy — %s", locus, discrepancy)

        # QC: OptiType vs serology validation
        ser_allele1 = ser_allele2 = ""
        serology_validation = ""
        if serology and locus in serology:
            ser_allele1, ser_allele2 = serology[locus]
            # Compare the best available OptiType call against serology
            best_optitype: list[str] | None = None
            if tumor_calls:
                best_optitype = tumor_calls[max(tumor_calls, key=lambda s: tumor_calls[s][1])][0]
            elif normal_calls:
                best_optitype = normal_calls[max(normal_calls, key=lambda s: normal_calls[s][1])][0]
            serology_validation = _validate_vs_serology(best_optitype, (ser_allele1, ser_allele2))
            if serology_validation:
                log.info("Locus %s serology validation: %s", locus, serology_validation)

        qc_rows.append({
            "locus": locus, "allele1": allele1, "allele2": allele2,
            "source": source, "reads": reads, "discrepancy": discrepancy,
            "serology_allele1": ser_allele1, "serology_allele2": ser_allele2,
            "serology_validation": serology_validation,
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
            fieldnames=[
                "locus", "allele1", "allele2", "source", "reads", "discrepancy",
                "serology_allele1", "serology_allele2", "serology_validation",
            ],
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
