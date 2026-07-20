"""Tests for the splice2neo adapter (#966, AC-2/AC-5).

Format mirrors the committed leaf-A smoke output
``research/experiments/issue_965_open_caller_audit/outputs/splice2neo_smoke_out.tsv``:
    junc_id  tx_id  frame_shift  cts_seq_len  peptide_context
"""

from adapters.splice2neo import parse_splice2neo


def test_maps_rows_to_common_records(tmp_path):
    tsv = tmp_path / "s2n.tsv"
    tsv.write_text(
        "junc_id\ttx_id\tframe_shift\tcts_seq_len\tpeptide_context\n"
        "chr2:152389996-152392205:-\tENST00000409198\tTRUE\t400\tINRHFKYATQLMNEIC\n"
        "2:100-200:+\tENST00000397345\tFALSE\t1945\tDMLTALYNSHMWSQVMSDGM\n"
    )
    records = parse_splice2neo(tsv, genome_build="hg19")

    assert len(records) == 2
    r0 = records[0]
    assert r0.caller == "splice2neo"
    assert r0.junction_id == "chr2:152389996-152392205:-"
    assert r0.genome_build == "hg19"
    assert r0.strand == "-"
    assert r0.peptide == "INRHFKYATQLMNEIC"
    assert r0.frame_shift is True
    assert r0.transcript_id == "ENST00000409198"


def test_normalizes_ensembl_contig_and_reads_plus_strand(tmp_path):
    tsv = tmp_path / "s2n.tsv"
    tsv.write_text(
        "junc_id\ttx_id\tframe_shift\tcts_seq_len\tpeptide_context\n"
        "2:100-200:+\tENST00000397345\tFALSE\t1945\tDMLTALYNSHMWSQVMSDGM\n"
    )
    (r,) = parse_splice2neo(tsv, genome_build="hg19")
    assert r.junction_id == "chr2:100-200:+"
    assert r.strand == "+"
    assert r.frame_shift is False


def test_skips_rows_with_no_peptide(tmp_path):
    # splice2neo emits NA peptide_context for junctions with no in-frame peptide;
    # those are not neoepitope candidates and must not enter the schema.
    tsv = tmp_path / "s2n.tsv"
    tsv.write_text(
        "junc_id\ttx_id\tframe_shift\tcts_seq_len\tpeptide_context\n"
        "chr2:152389996-152392205:-\tENST00000409198\tTRUE\t400\tNA\n"
        "chr3:1-2:+\tENST00000000001\tFALSE\t10\tMEIC\n"
    )
    records = parse_splice2neo(tsv, genome_build="hg19")
    assert [r.peptide for r in records] == ["MEIC"]
