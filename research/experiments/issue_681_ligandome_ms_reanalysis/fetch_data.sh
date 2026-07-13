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

echo "data/ rebuilt."
