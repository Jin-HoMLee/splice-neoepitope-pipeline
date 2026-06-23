#!/usr/bin/env python3
"""Apply the #566 option-B (MHCflurry-swap) modifications to a clone of bm2-lab/ASNEO.

Option B holds the MHC step constant on MHCflurry, so ASNEO's bundled
NetMHCpan/NetCTLpan/pepmatch binaries are bypassed entirely. This patcher makes
ASNEO stop right after it writes its junction-derived, normal-subtracted candidate
peptides (`putative_peptide.txt`, the `peps - norm_peps` set), and copies that file
to the output dir before the tmp dir is cleaned. Those candidates are then scored
with MHCflurry downstream (outside ASNEO), apples-to-apples with our pipeline.

String-anchored (not line-numbered) so it survives minor upstream drift; each
replacement must match exactly once or the script aborts (fail-loud).

Usage:  python apply_optionB_patch.py /path/to/ASNEO/ASNEO.py
"""
import sys, py_compile

REPLACEMENTS = [
    # 1. Use PATH/conda bedtools -> no need to extract the bundled src/software.tar.gz
    #    (which carries the non-redistributable MHC binaries). Fully open-only.
    (
        "    path['bedtools'] = os.path.join(path['sdir'], 'src/bedtools')",
        "    path['bedtools'] = 'bedtools'  # option B: PATH/conda bedtools; skip bundled src/",
    ),
    # 2. Stop ProcessIsoform after the candidate peptides are written; copy them out;
    #    skip the bundled netMHCpan call + ParseAffit.
    (
        "    logging.info(\"Call netMHCpan\")\n"
        "    RunMHCpan(allele=allele, peptxt=path['pep_txt'], affit=path['affit'])\n"
        "    ParseAffit(protfa=path['prot_fa'], affit=path['affit'], epit=path['epit'], bind=bind)",
        "    # option B (MHCflurry-swap): stop after candidate peptides; skip bundled MHC binaries\n"
        "    shutil.copy(path['pep_txt'], os.path.join(os.path.dirname(path['neo_sorted']), 'putative_peptide.txt'))\n"
        "    logging.info(\"Option B: wrote %s candidate peptides; skipping netMHCpan/netCTLpan\", len(remain_pep))\n"
        "    return",
    ),
    # 3. Skip ProcessNeo in main (it runs the bypassed MHC scoring + xgboost).
    (
        "    ProcessNeo(args.process)",
        "    # ProcessNeo(args.process)  # option B: skipped (MHC scoring bypassed in ProcessIsoform)",
    ),
]


def main(target):
    with open(target) as f:
        src = f.read()
    for old, new in REPLACEMENTS:
        n = src.count(old)
        if n != 1:
            sys.exit(f"ABORT: anchor matched {n} times (expected 1):\n{old[:80]}...")
        src = src.replace(old, new)
    with open(target, "w") as f:
        f.write(src)
    py_compile.compile(target, doraise=True)
    print(f"OK: applied {len(REPLACEMENTS)} option-B replacements to {target}; compiles clean.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
