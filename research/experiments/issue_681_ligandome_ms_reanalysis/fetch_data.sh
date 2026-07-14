#!/usr/bin/env bash
# Rebuild data/ for the Issue 681 ligandome reanalysis. data/ is gitignored (repo-wide `data/` rule);
# outputs/ is the committed artifact. Run from the repo root.
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="$DIR/data"
mkdir -p "$DATA"

# --- HLA Ligand Atlas (benign tissue immunopeptidome, 223,246 peptides)
curl -fsSL --retry 3 -o "$DATA/hla_ligand_atlas_aggregated.tsv.gz" \
  "https://hla-ligand-atlas.org/rel/2020.12/aggregated.tsv.gz"
gzip -t "$DATA/hla_ligand_atlas_aggregated.tsv.gz"   # a truncated gz fails here

# --- SysteMHC Atlas v2.0, non-UniProt peptides (tumor samples; whole table embedded, DataTables paginates it client-side)
curl -fsSL --retry 3 -o "$DATA/systemhc_nonuniprot.html" \
  "https://systemhc.sjtu.edu.cn/Non-UniProt"

# A truncated fetch of this page is silent and costs ~5,000 peptides, so assert the terminator.
# The first pass on 2026-07-13 cut off mid-<td> at 3.2 MB and under-counted 8,911 peptides as 3,870.
if ! tail -c 32 "$DATA/systemhc_nonuniprot.html" | grep -q '</html>'; then
  echo "ERROR: systemhc_nonuniprot.html is truncated (no closing </html>). Re-run." >&2
  exit 1
fi
rows=$(grep -o '<tr' "$DATA/systemhc_nonuniprot.html" | wc -l | tr -d ' ')
echo "systemhc_nonuniprot.html: $rows table rows (expect 20002)"

# --- our predicted panel, from R2
# Creds live in the repo-root .env (note: the var is R2_ENDPOINT, not R2_ENDPOINT_URL), and boto3 is
# installed ONLY in research/.venv - not in the snakemake or tests env. python-dotenv is not installed
# anywhere, so source the .env here and let python read the process environment.
set -a
# shellcheck disable=SC1091
source "$(git rev-parse --show-toplevel)/.env"
set +a

research/.venv/bin/python - "$DATA" <<'PY'
import os, sys, boto3
from botocore.config import Config

dest = sys.argv[1]
s3 = boto3.client(
    "s3",
    endpoint_url=os.environ["R2_ENDPOINT"],
    aws_access_key_id=os.environ["R2_ACCESS_KEY_ID"],
    aws_secret_access_key=os.environ["R2_SECRET_ACCESS_KEY"],
    region_name="auto",
    config=Config(signature_version="s3v4"),
)
bucket = os.environ["R2_BUCKET"]
for p in ("001", "002"):
    key = f"results/patient_{p}/predictions/mhc_presentation.tsv"
    out = os.path.join(dest, f"patient_{p}_mhc.tsv")
    s3.download_file(bucket, key, out)
    print(f"fetched {key} -> {out}")
PY

# --- HGNC gene symbols: the drive list for harvest_caatlas.py. caAtlas is gene-indexed (no bulk
# endpoint), so the harvest is driven by a symbol list.
#
# The filter is load-bearing, not cosmetic: the HGNC complete set carries 45,021 symbols, but the
# harvest that produced the committed outputs used the 19,273 **protein-coding, Approved** ones.
# Fetching the unfiltered set would harvest a completely different gene space. Verified 2026-07-14:
# today's protein-coding+Approved set is 19,296 symbols and is a strict SUPERSET of the original
# 19,273 (HGNC has added 23 since; zero were removed). The extra genes are harmless - harvest_caatlas
# is resumable and simply appends them - but the peptide set will not be bit-identical to the
# committed one if you re-harvest against a newer HGNC release. That is annotation drift, not a bug.
curl -fsSL --retry 3 -o "$DATA/hgnc_complete_set.tsv" \
  "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.tsv"
research/.venv/bin/python - "$DATA" <<'PY'
import csv, os, sys
d = sys.argv[1]
with open(os.path.join(d, "hgnc_complete_set.tsv"), newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
syms = sorted({
    r["symbol"]
    for r in rows
    if r.get("symbol")
    and r.get("locus_group") == "protein-coding gene"
    and r.get("status") == "Approved"
})
with open(os.path.join(d, "hgnc_symbols.txt"), "w") as fh:
    fh.write("\n".join(syms) + "\n")
print(f"hgnc_symbols.txt: {len(syms):,} protein-coding approved symbols (harvest used 19,273)")
PY

# --- SNAF supplementary workbook: NOT fetched here, and deliberately not committed.
#
# Both build_panel.py and recover_presented.py hard-require it (SNAF_XLSX in build_panel.py). It is
# the publisher's 29 MB file, so we neither redistribute it nor scrape it. Obtain it by hand:
#
#   paper : Li et al., SNAF. Sci Transl Med (2024). DOI 10.1126/scitranslmed.ade2886
#   file  : Data S1-S15 (the supplementary workbook, ~29 MB .xlsx)
#   place : ~/Zotero/storage/WK4DHT6M/scitranslmed.ade2886_data_s1_to_s15.xlsx
#           (our Zotero item TZGPRIFK, attachment WK4DHT6M - any path works if you point
#            SNAF_XLSX at it; the Zotero location is just where ours lives)
#
# Without it both scripts raise FileNotFoundError on the first line of the panel build. Flagged in
# review of PR #1172: a third party running this script alone could not reproduce the panel.
if [ ! -f "$HOME/Zotero/storage/WK4DHT6M/scitranslmed.ade2886_data_s1_to_s15.xlsx" ]; then
  echo "NOTE: SNAF workbook not found - see the comment block above for how to obtain it." >&2
  echo "      data/ is otherwise complete; build_panel.py will fail until it is in place." >&2
fi

echo "data/ rebuilt."
