"""Common output schema for the open splice-caller benchmark harness (#966).

The harness runs several open callers on a shared input and collects their
junction + neoepitope-candidate calls into one schema so they can be compared
tool-to-tool (leaf B of the #679 epic). Two design points, both load-bearing:

- ``junction_id`` is normalized to a single canonical string so the same
  underlying junction called by two tools collides on one key. The canonical
  form is ``{chrom}:{start}-{end}:{strand}`` with a UCSC ``chr`` prefix and
  genomic coordinate order (``start < end``); ``strand`` is carried explicitly.
  This is the base case of the #1100 canonical-junction scheme (single
  junctions in one known build); #1100 extends it for the registry's
  multi-junction events and opaque identifiers.
- ``genome_build`` is recorded per record, not assumed: splice2neo and ASNEO
  default to hg19, our pipeline is hg38, so a build column is the honest
  substrate for the liftover the concordance leaves need.
"""

from dataclasses import dataclass, field
from typing import Optional


def normalize_junction_id(chrom: str, start: int, end: int, strand: str) -> str:
    """Return the canonical junction key ``{chrom}:{start}-{end}:{strand}``.

    Adds a UCSC ``chr`` prefix when absent and orders the two coordinates
    genomically so argument order never changes identity.
    """
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    if chrom == "chrMT":
        # UCSC names the mitochondrial contig "chrM"; Ensembl uses "MT". Without
        # this a tool emitting "MT" and one emitting "chrM" never collide.
        chrom = "chrM"
    lo, hi = (start, end) if start <= end else (end, start)
    return f"{chrom}:{lo}-{hi}:{strand}"


def parse_junction_string(s: str) -> tuple[str, int, int, str]:
    """Parse a ``chrom:start-end:strand`` string into its components.

    Handles the splice2neo ``junc_id`` form, e.g. ``chr2:152389996-152392205:-``.
    """
    chrom, coords, strand = s.rsplit(":", 2)
    start_s, end_s = coords.split("-")
    return chrom, int(start_s), int(end_s), strand


@dataclass
class CommonRecord:
    """One harmonized neoepitope-candidate call from a single caller.

    ``presentation_*`` fields are populated by the shared MHCflurry stage that
    scores every tool's peptides identically; they are ``None`` until that
    stage runs.
    """

    caller: str
    junction_id: Optional[str]
    genome_build: str
    strand: Optional[str]
    peptide: str
    record_level: str  # "junction" | "peptide"; orthogonal to event_type
    frame_shift: Optional[bool] = None
    event_type: Optional[str] = None
    transcript_id: Optional[str] = None
    presentation_class: Optional[str] = None
    presentation_score: Optional[float] = None
    presentation_percentile: Optional[float] = None
    provenance: dict = field(default_factory=dict)


def validate(record: "CommonRecord") -> None:
    """Enforce the legal record_level / field combinations.

    A ``junction``-level record must carry junction coordinates; a
    ``peptide``-level record must not (its junction linkage is unavailable by
    construction - see #1258). Raise loudly on any illegal combination so an
    adapter bug cannot smuggle a mislabeled record downstream.
    """
    if record.record_level == "junction":
        if record.junction_id is None or record.strand is None:
            raise ValueError(
                f"junction-level record for {record.caller!r} needs junction_id "
                f"and strand; got junction_id={record.junction_id!r}, strand={record.strand!r}"
            )
    elif record.record_level == "peptide":
        if record.junction_id is not None or record.strand is not None:
            raise ValueError(
                f"peptide-level record for {record.caller!r} must have null junction_id "
                f"and strand; got junction_id={record.junction_id!r}, strand={record.strand!r}"
            )
    else:
        raise ValueError(
            f"unknown record_level {record.record_level!r} for {record.caller!r}; "
            f"expected 'junction' or 'peptide'"
        )
