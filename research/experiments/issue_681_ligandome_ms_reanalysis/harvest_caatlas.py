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
    distill()


def distill():
    """caatlas_raw.jsonl (gitignored, ~60 MB) -> outputs/caatlas_peptides.txt.gz (committed, 0.9 MB).

    This is the step that was missing from the first pass: the committed `.gz` is the AC-2
    ceiling-control input (caAtlas returns 0 against SNAF's MS-confirmed splice peptides, one of the
    two legs of the headline finding), so its provenance chain has to be in the repo, not in
    somebody's shell history.

    The `.gz` is committed on purpose (0.9 MB, well under the 10 MB size band, and network-derived
    from a live site, so it is not offline-regenerable). The raw JSONL is not: it is 60 MB and this
    is its only consumer.

    Emits sorted unique peptides in the class-I 8-11mer band, which is exactly what build_panel's
    `load_caatlas` then reads back.
    """
    import gzip

    from build_panel import valid

    peps = []
    with RAW.open() as fh:
        for line in fh:
            rec = json.loads(line)
            for p in rec.get("peptides") or []:
                seq = p.get("Peptide_Sequence")
                if seq:
                    peps.append(seq)

    kept = sorted(valid(peps))
    dest = OUT / "caatlas_peptides.txt.gz"
    # mtime=0 so the gzip header carries no timestamp: re-distilling identical input yields an
    # identical file rather than a spurious diff.
    with gzip.GzipFile(dest, "wb", mtime=0) as gz:
        gz.write(("\n".join(kept) + "\n").encode())
    print(f"distilled {len(peps):,} peptide rows -> {len(kept):,} unique class-I peptides -> {dest}")


if __name__ == "__main__":
    if "--distill-only" in sys.argv:
        # Re-derive the committed .gz from an existing raw harvest, without re-hitting the site.
        sys.exit(distill())
    sys.exit(main(sys.argv[1]))
