"""splice2neo adapter (#966).

splice2neo (TRON, MIT) is a junction-to-peptide library, not an end-to-end
caller: it ingests an upstream junction caller's output and emits
``peptide_context`` sequences with no MHC step. Its smoke output
(leaf A, #965) is a TSV:

    junc_id                       tx_id            frame_shift  cts_seq_len  peptide_context
    chr2:152389996-152392205:-    ENST00000409198  TRUE         400          INRHFKYATQLMNEIC

Default build is hg19 (needs liftOver to our hg38 downstream).
"""

import csv

from common_schema import CommonRecord, normalize_junction_id, parse_junction_string

CALLER = "splice2neo"


def _to_bool(value):
    """Map splice2neo's ``frame_shift`` cell to a tri-state bool.

    Blank / ``NA`` -> ``None`` ("unknown"), preserving the schema's
    ``Optional[bool]`` rather than collapsing an absent value to ``False``.
    """
    v = (value or "").strip().upper()
    if v in ("", "NA"):
        return None
    return v == "TRUE"


def parse_splice2neo(path, genome_build: str = "hg19") -> list[CommonRecord]:
    """Parse a splice2neo peptide-context TSV into CommonRecords.

    Rows with no in-frame peptide (``peptide_context`` empty or ``NA``) are
    skipped - they are not neoepitope candidates.
    """
    records: list[CommonRecord] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            peptide = (row.get("peptide_context") or "").strip()
            if not peptide or peptide == "NA":
                continue
            chrom, start, end, strand = parse_junction_string(row["junc_id"])
            records.append(
                CommonRecord(
                    caller=CALLER,
                    junction_id=normalize_junction_id(chrom, start, end, strand),
                    genome_build=genome_build,
                    strand=strand,
                    peptide=peptide,
                    record_level="junction",
                    frame_shift=_to_bool(row.get("frame_shift", "")),
                    event_type="junction",
                    transcript_id=(row.get("tx_id") or None),
                    provenance={
                        "source": CALLER,
                        "cts_seq_len": row.get("cts_seq_len"),
                    },
                )
            )
    return records
