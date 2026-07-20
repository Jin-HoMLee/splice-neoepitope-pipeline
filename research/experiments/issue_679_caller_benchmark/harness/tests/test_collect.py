"""Tests for the harness collect CLI (#966, AC-1/AC-3/AC-5).

The collect step is the harness entry point: it takes one or more caller
outputs, dispatches each to its registered adapter, merges the results into
one common schema, and writes a unified TSV. Adding a caller is a one-line
registry entry (the AC-5 extension point).
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from collect import (
    ADAPTERS,
    _cell,
    _parse_input_spec,
    collect,
    get_adapter,
    records_to_tsv,
)

HARNESS_DIR = Path(__file__).resolve().parent.parent
SPLICE2NEO_TSV = (
    "junc_id\ttx_id\tframe_shift\tcts_seq_len\tpeptide_context\n"
    "chr2:152389996-152392205:-\tENST00000409198\tTRUE\t400\tINRHFKYATQLMNEIC\n"
)


def test_registry_maps_splice2neo_to_its_adapter():
    from adapters.splice2neo import parse_splice2neo

    assert ADAPTERS["splice2neo"] is parse_splice2neo


def test_get_adapter_rejects_unknown_caller_and_names_the_known_ones():
    with pytest.raises(ValueError) as exc:
        get_adapter("netmhcpan")
    msg = str(exc.value)
    assert "netmhcpan" in msg
    assert "splice2neo" in msg  # error lists the registered callers


def test_collect_dispatches_input_to_the_registered_adapter(tmp_path):
    tsv = tmp_path / "s2n.tsv"
    tsv.write_text(SPLICE2NEO_TSV)

    records = collect([("splice2neo", tsv)])

    assert len(records) == 1
    assert records[0].caller == "splice2neo"
    assert records[0].peptide == "INRHFKYATQLMNEIC"


def test_collect_merges_multiple_inputs(tmp_path):
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    a.write_text(SPLICE2NEO_TSV)
    b.write_text(SPLICE2NEO_TSV)

    records = collect([("splice2neo", a), ("splice2neo", b)])

    assert len(records) == 2


def test_records_to_tsv_writes_header_and_json_encoded_provenance(tmp_path):
    tsv = tmp_path / "s2n.tsv"
    tsv.write_text(SPLICE2NEO_TSV)
    records = collect([("splice2neo", tsv)])

    out = tmp_path / "unified.tsv"
    records_to_tsv(records, out)

    lines = out.read_text().splitlines()
    header = lines[0].split("\t")
    assert header[0] == "caller"
    assert "junction_id" in header
    assert "peptide" in header
    assert "provenance" in header

    row = dict(zip(header, lines[1].split("\t")))
    assert row["caller"] == "splice2neo"
    assert row["junction_id"] == "chr2:152389996-152392205:-"
    # provenance is a dict on the record; on disk it must round-trip as JSON.
    assert json.loads(row["provenance"])["source"] == "splice2neo"


def test_unified_tsv_round_trips_through_a_tab_csv_reader(tmp_path):
    # The provenance-as-JSON encoding is justified by the TSV reading back with
    # a tab reader (no csv-quoting surprises). Pin that load-bearing claim: a
    # real csv.DictReader must recover the JSON verbatim.
    tsv = tmp_path / "s2n.tsv"
    tsv.write_text(SPLICE2NEO_TSV)
    out = tmp_path / "unified.tsv"
    records_to_tsv(collect([("splice2neo", tsv)]), out)

    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["peptide"] == "INRHFKYATQLMNEIC"
    assert json.loads(rows[0]["provenance"])["source"] == "splice2neo"


def test_cell_rejects_a_value_with_an_embedded_tab():
    # A tab inside a value would silently shift every downstream column, so the
    # writer must refuse it loudly rather than corrupt the TSV.
    with pytest.raises(ValueError):
        _cell("pep\ttide")


def test_parse_input_spec_rejects_a_token_without_a_colon():
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_input_spec("splice2neo")  # no ":path"


def test_main_writes_a_unified_tsv_end_to_end(tmp_path):
    tsv = tmp_path / "s2n.tsv"
    tsv.write_text(SPLICE2NEO_TSV)
    out = tmp_path / "unified.tsv"

    result = subprocess.run(
        [
            sys.executable,
            str(HARNESS_DIR / "collect.py"),
            "--input",
            f"splice2neo:{tsv}",
            "--out",
            str(out),
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert out.exists()
    lines = out.read_text().splitlines()
    assert lines[0].startswith("caller\t")
    assert len(lines) == 2  # header + one record
