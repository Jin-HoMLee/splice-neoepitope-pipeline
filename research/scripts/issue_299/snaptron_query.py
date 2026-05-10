"""
Issue #299 — Snaptron PSR_GTEx validation for NEJ_GNAS / NEJ_RPL22.

Query Snaptron's GTEx (hg19, ~9166 samples) endpoint for all junctions in
GNAS and RPL22 gene regions. Identify candidate NEJs by molecular signature:
- NEJ_GNAS: A3 acceptor shift, -2 nt loss (frameshift + premature stop)
- NEJ_RPL22: A3 acceptor shift, -6 nt loss (in-frame, loses 2 AAs in alpha-helix)

Per Kwok et al. (Nature 2025): https://doi.org/10.1038/s41586-024-08552-0
"""
import csv
import sys
import urllib.request
from io import StringIO
from pathlib import Path

HERE = Path(__file__).parent

# Gene coords from UCSC hg19 (GRCh37, Ensembl 75 / GENCODE v33)
# Widened ~25 kb on each flank to cover all promoter variants + extended UTRs.
# GNAS has multiple promoters (XLas, NESP55, GNAS-AS1, A/B, 1A); RPL22 has alt 5'UTR.
GENE_REGIONS = {
    "GNAS":  ("chr20", 57380000, 57520000, "+"),
    "RPL22": ("chr1",   6230000,  6300000, "-"),
}

SNAPTRON_GTEX = "https://snaptron.cs.jhu.edu/gtex/snaptron"

def query(gene, chrom, start, end):
    url = f"{SNAPTRON_GTEX}?regions={chrom}:{start}-{end}&header=1"
    print(f"\n=== {gene}: {chrom}:{start}-{end} ===", file=sys.stderr)
    print(f"GET {url}", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=60) as r:
        text = r.read().decode("utf-8")
    rows = list(csv.DictReader(StringIO(text), delimiter="\t"))
    print(f"  -> {len(rows)} junctions", file=sys.stderr)
    return rows

if __name__ == "__main__":
    out = {}
    for gene, (chrom, s, e, strand) in GENE_REGIONS.items():
        out[gene] = query(gene, chrom, s, e)
    # Save raw responses
    for gene, rows in out.items():
        path = HERE / f"snaptron_{gene}_gtex_hg19.tsv"
        with open(path, "w") as f:
            if rows:
                w = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
                w.writeheader()
                for r in rows:
                    w.writerow(r)
        print(f"Wrote {len(rows)} rows -> {path}", file=sys.stderr)
