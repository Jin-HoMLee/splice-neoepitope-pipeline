#!/usr/bin/env python
"""Harvest the caAtlas immunopeptidome, gene by gene. Issue 681, AC-2.

caAtlas (iScience 2021, doi:10.1016/j.isci.2021.103107) is a **tumor** immunopeptidome built from 43
published datasets. It is the reference most likely to carry real splice-junction hits, and the one
that could give the enrichment result the power it currently lacks.

Two facts about how it is served, both of which cost me time:

1. **The portal is `www.zhang-lab.org/caatlas/`.** `caatlas.omicsbio.info` resolves and answers 200,
   but it is an unrelated lab's site whose `/index.php/gene/<GENE>` route silently returns the
   homepage. The paper's data-availability statement is the authority here.
2. **There is no bulk peptide table.** The download page offers exactly three files (CT antigens,
   cancer-associated antigens, PTM antigens - all gene-level or PTM-scoped). The core peptide set is
   reachable only per-gene, as inline JSON in a bootstrapTable `data: [...]` literal on
   `/index.php/gene/<GENE>`. Hence this harvester.

Each record carries a `source` field (`UniProt` vs otherwise), which is caAtlas's own
canonical/non-canonical split - the non-UniProt rows are where a splice-junction peptide would live.

Harvest over the full HGNC protein-coding symbol list rather than only our panel's genes: restricting
the harvest to SNAF's gene list would make it structurally impossible for a peptide in any other gene
to match, which would bake the answer into the query.

Resumable: appends to outputs/caatlas_raw.jsonl and skips genes already present.
"""

import json
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
RAW = OUT / "caatlas_raw.jsonl"

BASE = "https://www.zhang-lab.org/caatlas/index.php/gene/{}"
# The peptide table is the `data: [...]` literal fed to bootstrapTable for #table_norm.
DATA_RE = re.compile(r"data:\s*(\[\{.*?\}\])", re.S)
WORKERS = 8
TIMEOUT = 30


def parse_gene_page(text):
    """Pull the peptide records out of the inline bootstrapTable data literal."""
    out = []
    for blob in DATA_RE.findall(text):
        try:
            rows = json.loads(blob)
        except json.JSONDecodeError:
            continue
        for r in rows:
            if isinstance(r, dict) and "Peptide_Sequence" in r:
                out.append(r)
    return out


def fetch(gene, session, retries=3):
    for attempt in range(retries):
        try:
            r = session.get(BASE.format(gene), timeout=TIMEOUT)
            if r.status_code == 404:
                return gene, []
            r.raise_for_status()
            return gene, parse_gene_page(r.text)
        except Exception:
            if attempt == retries - 1:
                return gene, None  # None = failed, distinct from [] = no peptides
            time.sleep(1.5 * (attempt + 1))
    return gene, None


def main(symbols_file):
    genes = [g.strip() for g in Path(symbols_file).read_text().split("\n") if g.strip()]

    done = set()
    if RAW.exists():
        for line in RAW.open():
            try:
                done.add(json.loads(line)["gene"])
            except Exception:
                pass
    todo = [g for g in genes if g not in done]
    print(f"{len(genes):,} genes total, {len(done):,} already harvested, {len(todo):,} to go", flush=True)

    OUT.mkdir(exist_ok=True)
    session = requests.Session()
    session.headers["User-Agent"] = "splice-neoepitope-pipeline research (academic reuse)"

    n_pep = n_fail = 0
    with RAW.open("a") as fh, ThreadPoolExecutor(max_workers=WORKERS) as pool:
        for i, (gene, recs) in enumerate(pool.map(lambda g: fetch(g, session), todo), 1):
            if recs is None:
                n_fail += 1
                fh.write(json.dumps({"gene": gene, "error": True}) + "\n")
            else:
                n_pep += len(recs)
                fh.write(json.dumps({"gene": gene, "peptides": recs}) + "\n")
            if i % 250 == 0:
                fh.flush()
                print(f"  {i:,}/{len(todo):,} genes | {n_pep:,} peptide rows | {n_fail} failures", flush=True)

    print(f"done: {n_pep:,} peptide rows, {n_fail} failed genes -> {RAW}")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
