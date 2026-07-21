"""ASNEO adapter (#966).

ASNEO (bm2-lab, Apache-2.0) is an end-to-end splice-neoantigen caller. Its
open-only stop-point output (``putative_peptide.txt``, ``ASNEO.py:246-253``)
is a normal-subtracted SET of k-mer peptide strings - one bare peptide per
line. That set operation discards the junction -> peptide linkage, so these
records are peptide-level (``record_level="peptide"``, null junction fields).
Recovering the linkage is tracked in #1258 (a patch to ASNEO), NOT #1100.
"""

from typing import Optional

from common_schema import CommonRecord

CALLER = "asneo"


def parse_asneo(path, genome_build: str = "hg19",
                provenance: Optional[dict] = None) -> list[CommonRecord]:
    """Parse an ASNEO ``putative_peptide.txt`` (one bare peptide per line).

    Blank lines and ``#`` comment lines (our smoke annotations) are skipped.
    Each peptide becomes a peptide-level record with null junction fields.
    """
    base_prov = {"source": CALLER}
    if provenance:
        base_prov.update(provenance)
    records: list[CommonRecord] = []
    with open(path) as fh:
        for line in fh:
            peptide = line.strip()
            if not peptide or peptide.startswith("#"):
                continue
            records.append(
                CommonRecord(
                    caller=CALLER,
                    junction_id=None,
                    genome_build=genome_build,
                    strand=None,
                    peptide=peptide,
                    record_level="peptide",
                    event_type="junction",
                    provenance=dict(base_prov),
                )
            )
    return records
