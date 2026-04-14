#!/usr/bin/env python3
"""aggregate_hla_alleles.py — Aggregate per-sample arcasHLA calls into a
per-patient alleles TSV.

Normal-first policy
-------------------
HLA alleles are germline, so normal tissue gives a cleaner signal than tumor
(which may carry CNV or LOH events at the HLA locus). For each locus (A, B, C)
the aggregator:

  1. Collects typed alleles from the Solid Tissue Normal sample if it passes
     the min_reads threshold.
  2. Falls back to the Primary Tumor sample if no normal call is available.
  3. Falls back to config-defined fallback alleles if no sample call is
     available or confident.

QC output
---------
A side TSV is written with one row per locus recording:
  - source     (normal / tumor / fallback)
  - reads      (read depth at that locus from arcasHLA genotype.log)
  - discrepancy (if normal and tumor calls differ)

Input (via Snakemake)
---------------------
  snakemake.input.jsons    — list of {sample}.genotype.json paths
  snakemake.input.genologs — list of {sample}.genotype.log paths
  snakemake.output.alleles — output alleles TSV path
  snakemake.output.qc      — output QC TSV path
  snakemake.params.samples_tsv  — path to the samples TSV
  snakemake.params.loci         — list of loci, e.g. ["A", "B", "C"]
  snakemake.params.min_reads    — minimum reads per locus for confidence
  snakemake.params.fallback     — dict of fallback alleles per locus
"""

import csv
import json
import logging
import re
from pathlib import Path


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_samples(tsv_path: str) -> dict[str, str]:
    """Return {sample_id: sample_type} from the samples TSV."""
    result = {}
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sid = (row.get("sample_id") or "").strip()
            if not sid or sid.startswith("#"):
                continue
            result[sid] = (row.get("sample_type") or "Unknown").strip()
    return result


def load_genotype(json_path: str) -> dict[str, list[str]]:
    """Load arcasHLA 0.5.x genotype JSON.

    arcasHLA 0.5.x writes keys like "HLA-A" with values like
    ["A*02:01:01", "A*24:02:01"]. Returns {} on missing / empty / invalid file.
    """
    path = Path(json_path)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    try:
        data = json.loads(path.read_text())
        return data if isinstance(data, dict) else {}
    except json.JSONDecodeError as exc:
        log.warning("Could not parse %s: %s", json_path, exc)
        return {}


def read_locus_reads(log_path: str, locus: str) -> int:
    """Parse arcasHLA genotype.log for read depth at a given locus.

    arcasHLA 0.5.x writes lines like:
      [genotype] HLA-A: 1234 reads
    Returns 0 if the locus line is not found.
    """
    path = Path(log_path)
    if not path.exists():
        return 0
    pattern = re.compile(rf"HLA-{re.escape(locus)}.*?(\d+)\s+reads", re.IGNORECASE)
    for line in path.read_text().splitlines():
        m = pattern.search(line)
        if m:
            return int(m.group(1))
    return 0


def normalise_allele(raw: str) -> str:
    """Return allele in HLA-X*YY:ZZ format (4-digit resolution).

    arcasHLA may return full 8-digit resolution (A*02:01:01:01); truncate
    to the first two fields which is what MHCflurry expects.
    """
    raw = raw.strip()
    if not raw.startswith("HLA-"):
        raw = f"HLA-{raw}"
    parts = raw.split(":")
    return ":".join(parts[:2])


# ---------------------------------------------------------------------------
# Main aggregation
# ---------------------------------------------------------------------------

def aggregate(
    jsons: list[str],
    genologs: list[str],
    samples_tsv: str,
    loci: list[str],
    min_reads: int,
    fallback: dict[str, str],
    out_alleles: str,
    out_qc: str,
) -> None:
    sample_types = read_samples(samples_tsv)
    log.info("Loaded %d samples from %s", len(sample_types), samples_tsv)

    # Map sample IDs to their json/log paths
    json_by_sample = {}
    log_by_sample = {}
    for jpath, lpath in zip(jsons, genologs):
        sample_id = Path(jpath).name.replace(".genotype.json", "")
        json_by_sample[sample_id] = jpath
        log_by_sample[sample_id] = lpath

    normal_samples = [s for s in json_by_sample if sample_types.get(s) == "Solid Tissue Normal"]
    tumor_samples  = [s for s in json_by_sample if sample_types.get(s) == "Primary Tumor"]
    log.info("Normal samples: %s", normal_samples)
    log.info("Tumor samples:  %s", tumor_samples)

    allele_rows = []
    qc_rows = []

    for locus in loci:
        hla_key = f"HLA-{locus}"

        def _get_calls(sample_list):
            calls = {}
            for sid in sample_list:
                geno = load_genotype(json_by_sample[sid])
                reads = read_locus_reads(log_by_sample[sid], locus)
                alleles = geno.get(hla_key, [])
                if alleles and reads >= min_reads:
                    calls[sid] = ([normalise_allele(a) for a in alleles], reads)
                elif alleles:
                    log.warning(
                        "Sample %s locus %s: only %d reads (min %d), skipping",
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
            log.warning("Locus %s: no confident call, using fallback allele '%s'", locus, fb)

        # Pad to two alleles (homozygous or missing second allele)
        allele1 = chosen_alleles[0] if len(chosen_alleles) > 0 else fallback.get(locus, "")
        allele2 = chosen_alleles[1] if len(chosen_alleles) > 1 else allele1

        allele_rows.append({"locus": locus, "allele1": allele1, "allele2": allele2})

        # QC: check normal/tumor concordance
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
            "locus":       locus,
            "allele1":     allele1,
            "allele2":     allele2,
            "source":      source,
            "reads":       reads,
            "discrepancy": discrepancy,
        })

    # Write outputs
    Path(out_alleles).parent.mkdir(parents=True, exist_ok=True)
    with open(out_alleles, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["locus", "allele1", "allele2"], delimiter="\t")
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
        jsons=snakemake.input.jsons,        # noqa: F821
        genologs=snakemake.input.genologs,  # noqa: F821
        samples_tsv=snakemake.params.samples_tsv,
        loci=snakemake.params.loci,
        min_reads=snakemake.params.min_reads,
        fallback=snakemake.params.fallback,
        out_alleles=snakemake.output.alleles,
        out_qc=snakemake.output.qc,
    )


try:
    snakemake  # noqa: F821
    _snakemake_main()
except NameError:
    pass  # imported directly (e.g. for testing)
