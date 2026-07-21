"""Tests for the ASNEO adapter (#966, AC-2/AC-5).

ASNEO's open-only stop-point output is one bare candidate peptide per line
(putative_peptide.txt); records are peptide-level with null junction fields.
"""

from adapters.asneo import parse_asneo


def test_parses_bare_peptides_as_peptide_level(tmp_path):
    f = tmp_path / "putative_peptide.txt"
    f.write_text(
        "# ASNEO chr22 smoke - first N candidate peptides\n"
        "LEQGTHPKFQ\n"
        "\n"
        "ILQPKPVD\n"
    )
    records = parse_asneo(f, genome_build="hg19",
                          provenance={"sj_tab_digest": "abc123", "thresholds": "reads2 psi0.05"})
    assert [r.peptide for r in records] == ["LEQGTHPKFQ", "ILQPKPVD"]
    r0 = records[0]
    assert r0.caller == "asneo"
    assert r0.record_level == "peptide"
    assert r0.junction_id is None and r0.strand is None
    assert r0.genome_build == "hg19"
    assert r0.event_type == "junction"
    assert r0.provenance["source"] == "asneo"
    assert r0.provenance["sj_tab_digest"] == "abc123"
