"""AC-4 payoff cross-reference: both snowball catalogs vs our n=2156 candidate presenters.

Reads the committed catalog TSVs (peptide column) rather than hardcoding peptide
lists, so the exitron AND ovarian halves are both reproducible and cannot drift
from the artifacts. Presenter peptides are pulled live from R2.

Run: research/.venv/bin/python xref_catalogs_vs_presenters.py   (R2 creds from project-root .env)
"""
import os, io, csv, boto3
from botocore.config import Config

HERE = os.path.dirname(os.path.abspath(__file__))
CATALOGS = {
    "Exitron (Wang 2021, S4+S5)": "exitron_ms_presented_S4_S5.tsv",
    "Ovarian (Zhao 2020, S5)":    "ovarian_tsa_S5.tsv",
}

def load_peptides(fname):
    path = os.path.join(HERE, fname)
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        peps = [(r.get("peptide") or "").strip().upper() for r in rd]
    return [p for p in peps if p]

# --- pull presenter peptides from R2 ---
ep = os.environ["R2_ENDPOINT"]; bk = os.environ["R2_BUCKET"]
s3 = boto3.client("s3", endpoint_url=ep, aws_access_key_id=os.environ["R2_ACCESS_KEY_ID"],
    aws_secret_access_key=os.environ["R2_SECRET_ACCESS_KEY"],
    config=Config(signature_version="s3v4", region_name="auto"))

presenters = []   # (peptide, patient)
for pt in ["patient_001", "patient_002"]:
    key = f"results/{pt}/predictions/mhc_presentation.tsv"
    obj = s3.get_object(Bucket=bk, Key=key)
    rd = csv.DictReader(io.StringIO(obj["Body"].read().decode("utf-8")), delimiter="\t")
    n = 0
    for row in rd:
        pep = (row.get("peptide") or "").strip().upper()
        if pep:
            presenters.append((pep, pt)); n += 1
    print(f"  {pt}: {n} presenter rows")

pep_set = {p[0] for p in presenters}
print(f"TOTAL presenter peptides: {len(presenters)} rows, {len(pep_set)} unique\n")

# --- per-catalog exact + substring (either direction) intersection ---
for label, fname in CATALOGS.items():
    cat = load_peptides(fname)
    exact = [c for c in cat if c in pep_set]
    subs = []
    for c in cat:
        for pep, pt in presenters:
            if pep != c and (pep in c or c in pep):
                subs.append((c, pep, pt))
    print(f"=== {label}: {len(cat)} peptides ===")
    print(f"    EXACT overlap:     {exact if exact else 'NONE'}")
    print(f"    SUBSTRING overlap: {subs if subs else 'NONE'}")
    print()

# --- ceiling control: presenters must actually contain class-I-length peptides ---
lens = sorted({len(p) for p in pep_set})
print(f"presenter peptide lengths present: {lens}  (class-I 8/9/10-mers must appear for the null to be real)")
