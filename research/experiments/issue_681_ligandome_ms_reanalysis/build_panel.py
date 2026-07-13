#!/usr/bin/env python
"""Assemble the AC-1 splice-neoepitope panel and intersect it against public HLA-ligandome atlases.

Issue #681. Run with `research/.venv/bin/python` (needs openpyxl; boto3 only if re-fetching from R2).

The panel has three provenance strata:

  ours          our pipeline's predicted presenters, 2 patients (R2 mhc_presentation.tsv)
  snaf_pred     SNAF-T predicted splice neoantigens (Data S1 TCGA-SKCM, S1 Van Allen, S3 ovarian)
  snaf_ms       SNAF peptides *detected in HLA-ligandome MS* (Data S2, Mel cohort)

`snaf_ms` is the load-bearing stratum. It is not extra coverage - it is the **ceiling control**.
Those peptides are splice-junction-derived AND known to be MS-detectable, because SNAF detected
them in real immunopeptidomics. So if they too are absent from a reference atlas, that atlas
*cannot report this peptide class at all*, and any null we get against it is guaranteed by the
method rather than observed in the biology. See feedback_search_key_must_not_mirror_the_gate.md:
check the CEILING, not just the floor.
"""

import os
import re
import gzip
import csv
import json
import html
import sys
from pathlib import Path

AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "outputs"


def valid(seqs, lo=8, hi=11):
    """Keep canonical-alphabet peptides in the class-I length band."""
    return {s for s in (str(x).strip().upper() for x in seqs) if AA.match(s) and lo <= len(s) <= hi}


# ---------------------------------------------------------------- panel strata

def load_snaf(xlsx):
    """SNAF supplementary workbook -> (predicted, ms_detected).

    Sheet layout is not uniform: the Data S1/S3 sheets carry a title banner, so their real
    header is row 3 and peptides start at row 4. Data S2 has no banner (header row 1).
    """
    import openpyxl

    wb = openpyxl.load_workbook(xlsx, read_only=True)

    def column(sheet, first_data_row, col_idx):
        ws = wb[sheet]
        return [
            r[col_idx]
            for r in ws.iter_rows(min_row=first_data_row, values_only=True)
            if r[col_idx] is not None
        ]

    predicted = valid(
        column("Data S1 - TCGA", 4, 0)
        + column("Data S1 - Van Allen", 4, 0)
        + column("Data S3", 4, 0)
    )
    ms = valid(column("Data S2", 2, 1))
    return predicted, ms


def load_ours(tsvs):
    """Our predicted panel: every peptide we scored, regardless of presentation_class."""
    peps = []
    for t in tsvs:
        with open(t) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("peptide"):
                    peps.append(row["peptide"])
    return valid(peps)


# ------------------------------------------------------------ reference atlases

def load_hla_ligand_atlas(gz):
    """HLA Ligand Atlas aggregated.tsv.gz -> peptides (benign tissue, class I + II)."""
    peps = []
    with gzip.open(gz, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            peps.append(row["peptide_sequence"])
    return valid(peps)


def load_systemhc_nonuniprot(html_path, with_meta=False):
    """SysteMHC v2.0 non-UniProt peptides from the /Non-UniProt page HTML.

    The page embeds the whole table and paginates it client-side with DataTables, so the
    peptides are all in the HTML. Columns are fixed-width (10 per row), peptide at index 2:
    SysteMHC ID | Sample ID | Peptide | Allele | NetMHCpan Rank | Tissue | Cell | Disease | Binder | DB

    Do NOT sweep every <td> indiscriminately - that mixes sample IDs and protein names into the
    peptide set. Parse by column index.

    The served table is exactly 20,000 rows (15,130 unique peptides), a suspiciously round number
    that is very likely a server-side cap, so even this is a lower bound on what is exposed - and
    far below the 78,959 non-UniProt peptides Table 1 of the paper claims. AC-3 records the gap.
    """
    txt = html.unescape(Path(html_path).read_text(errors="ignore"))
    strip = lambda x: re.sub(r"<[^>]+>", "", x).strip()
    recs = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", txt, re.S):
        cells = [strip(c) for c in re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S)]
        if len(cells) >= 10:
            recs.append(
                {
                    "peptide": cells[2].upper(),
                    "allele": cells[3],
                    "tissue": cells[5],
                    "disease": cells[7],
                }
            )
    if with_meta:
        return recs
    return valid(r["peptide"] for r in recs)


def load_caatlas(gz):
    """caAtlas peptides, harvested per-gene by harvest_caatlas.py (see its docstring).

    **Ceiling warning, measured not assumed:** all 557,450 harvested records carry
    `source: UniProt` - there is not one non-UniProt row. caAtlas's gene-indexed surface can
    therefore only report peptides that map to a canonical protein, so a *genuinely novel*
    junction-spanning peptide is unreportable there by construction. A null against caAtlas is
    method-guaranteed for the novel class and must never be read as biological absence.
    """
    with gzip.open(gz, "rt") as fh:
        return valid(line.strip() for line in fh if line.strip())


def main():
    snaf_xlsx = Path(
        os.path.expanduser(
            "~/Zotero/storage/WK4DHT6M/scitranslmed.ade2886_data_s1_to_s15.xlsx"
        )
    )
    snaf_pred, snaf_ms = load_snaf(snaf_xlsx)
    ours = load_ours(sorted(DATA.glob("patient_*_mhc.tsv")))

    atlas = load_hla_ligand_atlas(DATA / "hla_ligand_atlas_aggregated.tsv.gz")
    systemhc = load_systemhc_nonuniprot(DATA / "systemhc_nonuniprot.html")
    caatlas = load_caatlas(OUT / "caatlas_peptides.txt.gz")

    panel = {"ours": ours, "snaf_pred": snaf_pred, "snaf_ms": snaf_ms}
    # Ordered by search database, which is the axis that turns out to matter: the first two are
    # canonical-proteome searches, the third is the only non-canonical one.
    refs = {
        "hla_ligand_atlas": atlas,      # benign tissue, canonical search
        "caatlas": caatlas,             # tumor, 100% UniProt-mapped (canonical by construction)
        "systemhc_nonuniprot": systemhc,  # tumor, NON-canonical by construction
    }

    print("PANEL (unique peptides)")
    for k, v in panel.items():
        print(f"  {k:12s} {len(v):>7,}")
    print("REFERENCE ATLASES (unique peptides)")
    for k, v in refs.items():
        print(f"  {k:20s} {len(v):>7,}")

    # --- control 1: the join works at all (two reference atlases must overlap each other)
    join_ctrl = len(atlas & systemhc)
    print(f"\nCONTROL 1 - join machinery: hla_ligand_atlas n systemhc = {join_ctrl:,}")
    print("  (a nonzero value proves exact string matching works across sources)")

    # --- control 2: the CEILING. Can these atlases report a splice-junction peptide at all?
    print("\nCONTROL 2 - ceiling: SNAF's MS-CONFIRMED splice peptides vs each atlas")
    print("  These are splice-junction peptides that real immunopeptidomics did detect.")
    print("  If an atlas returns ~0 here, it structurally cannot report this peptide class,")
    print("  and any null we score against it is method-guaranteed, not observed.")

    rows = []
    for pname, pset in panel.items():
        for rname, rset in refs.items():
            hits = sorted(pset & rset)
            rows.append(
                {
                    "panel": pname,
                    "reference": rname,
                    "panel_n": len(pset),
                    "reference_n": len(rset),
                    "hits": len(hits),
                    "hit_peptides": ";".join(hits[:50]),
                }
            )
            print(f"  {pname:12s} n {rname:20s} = {len(hits):>5,} / {len(pset):,}")

    OUT.mkdir(exist_ok=True)
    with open(OUT / "intersections.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # The panel itself, one row per unique peptide, carrying provenance and every atlas hit.
    # A peptide can be in more than one stratum (SNAF's MS set is a subset of its predictions),
    # so sources are joined rather than made exclusive.
    with open(OUT / "panel.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["peptide", "length", "sources", "in_hla_ligand_atlas", "in_systemhc_nonuniprot"])
        for pep in sorted(ours | snaf_pred | snaf_ms):
            src = [n for n, s in panel.items() if pep in s]
            w.writerow(
                [
                    pep,
                    len(pep),
                    ",".join(src),
                    int(pep in atlas),
                    int(pep in systemhc),
                ]
            )

    with open(OUT / "panel_summary.json", "w") as fh:
        json.dump(
            {
                "panel_sizes": {k: len(v) for k, v in panel.items()},
                "reference_sizes": {k: len(v) for k, v in refs.items()},
                "control_join_atlas_x_systemhc": join_ctrl,
                "snaf_ms_contained_in_snaf_pred": len(snaf_ms & snaf_pred),
                "snaf_ms_total": len(snaf_ms),
            },
            fh,
            indent=2,
        )
    print(f"\nwrote {OUT/'intersections.tsv'} and {OUT/'panel_summary.json'}")


if __name__ == "__main__":
    sys.exit(main())
