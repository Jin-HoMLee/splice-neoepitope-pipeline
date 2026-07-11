"""Unit tests for the UCSC RepeatMasker fetcher (Issue #919).

The load-bearing tests here are the truncation guards. A capped UCSC response
yields a valid, well-formed, *short* BED - it looks exactly like a successful
fetch, and would silently understate repeat overlap in the #919 analysis. So
truncation must raise, not warn.

Network is never touched: every test drives the pure functions off a synthetic
payload shaped like a real UCSC ``/getData/track`` response.
"""
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))

from fetch_rmsk import (  # noqa: E402
    BED_COLUMNS,
    RmskTruncatedError,
    build_track_url,
    extract_items,
    items_to_bed_rows,
    write_bed,
)


def _item(start, end, name="AluSx1", cls="SINE", family="Alu", strand="+", score=2021, chrom="chr22"):
    """One rmsk record, shaped like the real UCSC API payload."""
    return {
        "genoName": chrom,
        "genoStart": start,
        "genoEnd": end,
        "repName": name,
        "swScore": score,
        "strand": strand,
        "repClass": cls,
        "repFamily": family,
    }


def _payload(items, **extra):
    payload = {"chrom": "chr22", "rmsk": items, "itemsReturned": len(items)}
    payload.update(extra)
    return payload


# -- URL ----------------------------------------------------------------------

def test_build_track_url_uses_semicolon_separators():
    """UCSC's API wants `;`-separated params; `&` silently returns the wrong slice."""
    url = build_track_url("hg38", "chr22")
    assert url == (
        "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=rmsk;chrom=chr22"
    )
    assert "&" not in url


# -- truncation guards --------------------------------------------------------

def test_truncated_response_raises():
    """`maxItemsLimit` is UCSC's truncation marker - it must be fatal, not a warning."""
    payload = _payload([_item(100, 200)], maxItemsLimit=1)
    with pytest.raises(RmskTruncatedError, match="capped"):
        extract_items(payload, "chr22")


def test_complete_response_does_not_raise():
    """The marker is absent on a complete response - that is the only passing shape."""
    payload = _payload([_item(100, 200), _item(300, 400)])
    assert len(extract_items(payload, "chr22")) == 2


def test_item_count_disagreeing_with_itemsReturned_raises():
    """A payload we do not understand is not something to paper over."""
    payload = _payload([_item(100, 200)])
    payload["itemsReturned"] = 99
    with pytest.raises(RmskTruncatedError, match="itemsReturned"):
        extract_items(payload, "chr22")


def test_missing_track_key_raises():
    payload = {"chrom": "chr22", "itemsReturned": 0}
    with pytest.raises(KeyError):
        extract_items(payload, "chr22")


# -- BED conversion -----------------------------------------------------------

def test_coordinates_pass_through_unchanged():
    """UCSC genoStart/genoEnd are already 0-based half-open (BED convention).

    Regression guard against importing the regtools BED12 anchor-outer trap
    (Issue #370): there is no coordinate arithmetic to do here, and adding some
    would shift every repeat interval.
    """
    rows = items_to_bed_rows([_item(10510227, 10510528)], "chr22")
    assert rows[0][1] == "10510227"
    assert rows[0][2] == "10510528"


def test_row_has_bed6_plus_class_and_family():
    rows = items_to_bed_rows([_item(100, 200, name="L1MC5a", cls="LINE", family="L1", strand="-")], "chr22")
    assert len(rows[0]) == len(BED_COLUMNS) == 8
    assert rows[0] == ["chr22", "100", "200", "L1MC5a", "2021", "-", "LINE", "L1"]


def test_rows_are_coordinate_sorted():
    """Downstream overlap uses a sorted sweep, so order is load-bearing, not cosmetic."""
    rows = items_to_bed_rows(
        [_item(900, 950), _item(100, 200), _item(500, 600), _item(100, 150)], "chr22"
    )
    assert [int(r[1]) for r in rows] == [100, 100, 500, 900]
    assert [int(r[2]) for r in rows] == [150, 200, 600, 950]


def test_foreign_chromosome_records_dropped():
    """UCSC filters server-side; this only fires if the response ever widens."""
    rows = items_to_bed_rows([_item(100, 200), _item(300, 400, chrom="chr21")], "chr22")
    assert len(rows) == 1
    assert rows[0][0] == "chr22"


def test_missing_optional_fields_default_to_placeholder():
    """A record without repClass must not crash the fetch or shift the columns."""
    rows = items_to_bed_rows([{"genoName": "chr22", "genoStart": 1, "genoEnd": 2}], "chr22")
    assert rows[0] == ["chr22", "1", "2", ".", "0", ".", ".", "."]


# -- output -------------------------------------------------------------------

def test_write_bed_creates_parent_directory(tmp_path):
    out = tmp_path / "references" / "rmsk" / "hg38" / "rmsk.chr22.bed"
    write_bed(items_to_bed_rows([_item(100, 200)], "chr22"), out)
    assert out.read_text() == "chr22\t100\t200\tAluSx1\t2021\t+\tSINE\tAlu\n"
