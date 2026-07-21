#!/usr/bin/env python3
"""make_spectrum.py - generate a theoretical MS2 spectrum for a real junction peptide.

Picks a genuine junction ORF-stretch peptide from the emitted junction FASTA and
writes a synthetic MGF whose precursor and b/y fragment ions are computed from
that peptide's monoisotopic masses. This lets the Sage smoke recover an actual
junction peptide (a check that can fail), not merely prove the FASTA loads.

The target is derived from the *current* FASTA at run time, so the smoke stays
valid as the chr22 junction set evolves - no hand-computed masses, no hardcoded
peptide that could drift out of the database.

Prints one line to stdout: the target peptide sequence (for the caller's
recovery assertion).
"""

import argparse
import sys

# Monoisotopic residue masses (Da).
_RESIDUE = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276,
    "V": 99.06841, "T": 101.04768, "C": 103.00919, "L": 113.08406,
    "I": 113.08406, "N": 114.04293, "D": 115.02694, "Q": 128.05858,
    "K": 128.09496, "E": 129.04259, "M": 131.04049, "H": 137.05891,
    "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}
_WATER = 18.0105646
_PROTON = 1.00727646

_MIN_LEN = 9   # 9-mer -> 2*(9-1)=16 fragment peaks, clears Sage's min_peaks=15
_MAX_LEN = 11  # config peptide_max_len


def _parse_fasta(path):
    header = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def pick_target(fasta_path):
    """Return (header, peptide): the first stretch of length in [_MIN_LEN, _MAX_LEN]
    whose residues are all standard, else a _MIN_LEN-prefix of the first stretch
    at least _MIN_LEN long."""
    fallback = None
    for header, seq in _parse_fasta(fasta_path):
        if not set(seq) <= set(_RESIDUE):
            continue
        if _MIN_LEN <= len(seq) <= _MAX_LEN:
            return header, seq
        if fallback is None and len(seq) >= _MIN_LEN:
            fallback = (header, seq[:_MIN_LEN])
    if fallback is not None:
        return fallback
    raise SystemExit(f"no usable target peptide (>= {_MIN_LEN} standard-AA) in {fasta_path}")


def peptide_mass(pep):
    return sum(_RESIDUE[a] for a in pep) + _WATER


def fragment_ions(pep):
    """Return sorted singly-charged b- and y-ion m/z values."""
    ions = []
    # b-ions: prefix masses + proton
    running = 0.0
    for a in pep[:-1]:
        running += _RESIDUE[a]
        ions.append(running + _PROTON)
    # y-ions: suffix masses + water + proton
    running = 0.0
    for a in reversed(pep[1:]):
        running += _RESIDUE[a]
        ions.append(running + _WATER + _PROTON)
    return sorted(ions)


def write_mgf(path, pep, charge=2):
    m = peptide_mass(pep)
    precursor_mz = (m + charge * _PROTON) / charge
    ions = fragment_ions(pep)
    with open(path, "w") as out:
        out.write("BEGIN IONS\n")
        out.write(f"TITLE=synthetic_junction_peptide_{pep}\n")
        out.write(f"PEPMASS={precursor_mz:.5f} 1000.0\n")
        out.write(f"CHARGE={charge}+\n")
        out.write("RTINSECONDS=300.0\n")
        for i, mz in enumerate(ions):
            # descending, arbitrary but nonuniform intensities
            out.write(f"{mz:.5f} {500.0 + (i % 7) * 90.0:.1f}\n")
        out.write("END IONS\n")


def main():
    p = argparse.ArgumentParser(description="Generate an MGF for a real junction peptide.")
    p.add_argument("--fasta", required=True, help="Emitted junction ORF FASTA")
    p.add_argument("--out", required=True, help="Output MGF path")
    p.add_argument("--charge", type=int, default=2)
    args = p.parse_args()

    header, pep = pick_target(args.fasta)
    write_mgf(args.out, pep, charge=args.charge)
    # stderr: human-readable provenance; stdout: the bare peptide for assertion
    print(f"target junction peptide {pep} from {header} ({len(fragment_ions(pep))} peaks)",
          file=sys.stderr)
    print(pep)


if __name__ == "__main__":
    main()
