"""Unit tests for the junction/repeat overlap analysis (Issue #919).

The interval logic is where this can go quietly wrong: a half-open boundary
error or a bisect that stops at the first hit would shift the repeat-overlap
rate a few percent in a way no assertion on the report would catch, and the
whole point of the analysis is to compare two rates that differ by a few
percent. So the boundary and nesting cases are tested directly.
"""
import pytest

# conftest.py puts the experiment dir (which holds junction_repeat_overlap.py) on
# sys.path - same pattern as the issue_547 experiment.
from junction_repeat_overlap import (  # noqa: E402
    Junction,
    Repeat,
    RepeatIndex,
    classify,
    compare,
    format_report,
    is_annotated,
    load_junctions,
    load_rmsk,
    overlaps_repeat,
    parse_junction_id,
    repeat_classes,
    repeat_rate,
)


def _index(*intervals):
    """RepeatIndex over (start, end[, class]) triples on chr22."""
    repeats = [
        Repeat(start=s, end=e, name="R", rep_class=(c if len(rest) else "SINE"), rep_family="Alu")
        for s, e, *rest in intervals
        for c in [rest[0] if rest else "SINE"]
    ]
    return RepeatIndex({"chr22": repeats})


# -- junction id parsing ------------------------------------------------------

def test_parse_junction_id_converts_donor_to_zero_based():
    """donor is 1-based in the id; the intron start is donor - 1.

    Mirrors filter_junctions._parse_junction_id. An off-by-one here silently
    shifts every splice site by a base, which at repeat boundaries flips the
    overlap call.
    """
    j = parse_junction_id("chr22:10711002:10744610:+", reads=3)
    assert (j.chrom, j.start, j.end, j.strand, j.reads) == ("chr22", 10711001, 10744610, "+", 3)


def test_junction_id_round_trips():
    assert parse_junction_id("chr22:100:200:-").junction_id == "chr22:100:200:-"


def test_splice_sites_are_the_intron_terminal_bases():
    j = parse_junction_id("chr22:101:200:+")   # intron = [100, 200)
    assert j.donor_site == 100      # first intronic base
    assert j.acceptor_site == 199   # last intronic base, not 200


@pytest.mark.parametrize("bad", ["chr22:100:200", "chr22", "", "chr22:1:2:3:4"])
def test_malformed_junction_id_raises(bad):
    with pytest.raises(ValueError):
        parse_junction_id(bad)


# -- interval lookup ----------------------------------------------------------

def test_position_inside_repeat_hits():
    assert _index((100, 200)).at("chr22", 150)


def test_interval_is_half_open():
    """BED is half-open: start is inside, end is not."""
    index = _index((100, 200))
    assert index.at("chr22", 100)          # start: inclusive
    assert index.at("chr22", 199)          # last base: inside
    assert not index.at("chr22", 200)      # end: exclusive
    assert not index.at("chr22", 99)


def test_nested_repeats_are_all_found():
    """Repeats nest (a SINE inside an LTR). A bisect that stops at the first hit
    would report only one, understating the class breakdown."""
    hits = _index((100, 500, "LTR"), (200, 300, "SINE")).at("chr22", 250)
    assert {r.rep_class for r in hits} == {"LTR", "SINE"}


def test_long_repeat_spanning_query_is_found_despite_later_starts():
    """The scan-back window must be the longest repeat, not a fixed guess: a long
    LINE starting well before the query is still spanning it."""
    hits = _index((1_000, 7_000, "LINE"), (5_900, 5_950, "SINE"), (5_960, 5_980, "SINE")).at("chr22", 6_000)
    assert [r.rep_class for r in hits] == ["LINE"]


def test_unknown_chromosome_returns_no_hits():
    assert _index((100, 200)).at("chr1", 150) == []


def test_empty_index_returns_no_hits():
    assert RepeatIndex({}).at("chr22", 150) == []


# -- junction classification --------------------------------------------------

def test_overlap_when_only_donor_in_repeat():
    j = Junction("chr22", 150, 400, "+", 1)      # donor 150 in [100,200), acceptor 399 outside
    assert overlaps_repeat(j, _index((100, 200)))


def test_overlap_when_only_acceptor_in_repeat():
    j = Junction("chr22", 50, 150, "+", 1)       # acceptor 149 in [100,200)
    assert overlaps_repeat(j, _index((100, 200)))


def test_no_overlap_when_repeat_only_spans_the_intron_body():
    """A repeat in the middle of an intron is unremarkable - introns are
    repeat-dense. Only a splice site *inside* a repeat is the signal."""
    j = Junction("chr22", 50, 400, "+", 1)       # sites at 50 and 399; repeat is [100,200)
    assert not overlaps_repeat(j, _index((100, 200)))


def test_repeat_classes_deduplicated_and_sorted():
    j = Junction("chr22", 150, 250, "+", 1)      # both sites land in repeats
    classes = repeat_classes(j, _index((100, 200, "SINE"), (200, 300, "SINE")))
    assert classes == ["SINE"]


def test_classify_returns_donor_and_acceptor_hits_separately():
    j = Junction("chr22", 150, 250, "+", 1)
    donor, acceptor = classify(j, _index((100, 200, "LINE"), (200, 300, "SINE")))
    assert [r.rep_class for r in donor] == ["LINE"]
    assert [r.rep_class for r in acceptor] == ["SINE"]


# -- set comparison -----------------------------------------------------------

def _jset(*ids):
    return {i: parse_junction_id(i, reads=1) for i in ids}


def test_compare_splits_lost_retained_gained():
    off = _jset("chr22:100:200:+", "chr22:300:400:+")
    on = _jset("chr22:100:200:+", "chr22:500:600:+")
    c = compare(off, on)
    assert [j.junction_id for j in c.lost] == ["chr22:300:400:+"]
    assert [j.junction_id for j in c.retained] == ["chr22:100:200:+"]
    assert [j.junction_id for j in c.gained] == ["chr22:500:600:+"]


def test_gained_is_empty_for_a_pure_subset():
    """The filter only removes alignments, so a real run must gain nothing."""
    off = _jset("chr22:100:200:+", "chr22:300:400:+")
    on = _jset("chr22:100:200:+")
    assert compare(off, on).gained == []


def test_repeat_rate_counts_and_fraction():
    index = _index((100, 200))
    junctions = [Junction("chr22", 150, 900, "+", 1), Junction("chr22", 800, 900, "+", 1)]
    hits, total, frac = repeat_rate(junctions, index)
    assert (hits, total, frac) == (1, 2, 0.5)


def test_repeat_rate_on_empty_set_is_zero_not_a_crash():
    assert repeat_rate([], _index((100, 200))) == (0, 0, 0.0)


# -- IO -----------------------------------------------------------------------

def test_load_rmsk_reads_bed6_plus_two(tmp_path):
    bed = tmp_path / "rmsk.bed"
    bed.write_text("chr22\t100\t200\tAluSx1\t2021\t+\tSINE\tAlu\n")
    index = load_rmsk(bed)
    assert len(index) == 1
    assert index.at("chr22", 150)[0].rep_family == "Alu"


def test_load_rmsk_rejects_a_short_bed(tmp_path):
    """A plain BED6 would silently read repClass out of a missing column."""
    bed = tmp_path / "short.bed"
    bed.write_text("chr22\t100\t200\tAluSx1\t2021\t+\n")
    with pytest.raises(ValueError, match="BED6"):
        load_rmsk(bed)


def test_load_junctions_reads_id_and_reads(tmp_path):
    tsv = tmp_path / "raw_junctions.tsv"
    tsv.write_text("chr22:100:200:+\t7\nchr22:300:400:-\t1\n")
    junctions = load_junctions(tsv)
    assert set(junctions) == {"chr22:100:200:+", "chr22:300:400:-"}
    assert junctions["chr22:100:200:+"].reads == 7


def test_report_names_the_enrichment(tmp_path):
    off = _jset("chr22:151:900:+", "chr22:801:900:+")
    on = _jset("chr22:801:900:+")
    report = format_report(compare(off, on), _index((100, 200)), "test")
    assert "enrichment" in report.lower()
    assert "lost to the filter" in report


# -- the composition-artifact guard -------------------------------------------
#
# The unstratified enrichment is the number that made the #919 chr22 result look
# like a mild win (1.13x) when the true within-pool enrichment was 0.98x. This
# script is explicitly reusable for the whole-genome re-run (#1095), so it must
# not hand the next operator that same confounded headline unqualified.

def test_report_without_annotation_warns_that_enrichment_is_confounded():
    off = _jset("chr22:151:900:+", "chr22:801:900:+")
    on = _jset("chr22:801:900:+")
    report = format_report(compare(off, on), _index((100, 200)), "test", annotated=None)
    assert "composition-confounded" in report
    assert "--annotated-bed" in report


def test_report_with_annotation_stratifies_and_flags_the_naive_number():
    off = _jset("chr22:151:900:+", "chr22:801:900:+")
    on = _jset("chr22:801:900:+")
    annotated = frozenset({("chr22", 800, 900, "+")})
    report = format_report(compare(off, on), _index((100, 200)), "test", annotated=annotated)
    assert "CONFOUNDED" in report
    assert "Stratified by annotation status" in report
    assert "Enrichment within the unannotated pool" in report
    assert "lost, unannotated" in report


def test_stratified_enrichment_dissolves_a_pure_composition_artifact():
    """The regression this whole guard exists for.

    Construct the #919 shape: every junction, lost or retained, is repeat-overlapping
    when unannotated and repeat-free when annotated. The filter is doing *nothing*
    repeat-selective - but the lost set is annotation-poor, so the naive ratio still
    reads as an enrichment. Stratified, it must come out at exactly 1.00x (none).
    """
    index = _index((100, 200), (300, 400))
    # lost: 1 unannotated (in repeat). retained: 1 unannotated (in repeat) + 1 annotated (not).
    off = _jset("chr22:151:900:+", "chr22:351:950:+", "chr22:501:990:+")
    on = _jset("chr22:351:950:+", "chr22:501:990:+")
    annotated = frozenset({("chr22", 500, 990, "+")})
    report = format_report(compare(off, on), index, "test", annotated=annotated)
    # Naive: lost 100% vs retained 50% -> 2.00x, a pure artifact of annotated content.
    assert "**2.00x** (CONFOUNDED" in report
    # Stratified within the unannotated pool: 100% vs 100% -> no enrichment.
    assert "Enrichment within the unannotated pool: 1.00x" in report


def test_enrichment_is_not_inf_on_a_degenerate_set():
    """A retained set with zero repeat overlap must read n/a, not `inf`."""
    off = _jset("chr22:151:900:+", "chr22:801:900:+")
    on = _jset("chr22:801:900:+")
    report = format_report(compare(off, on), _index((100, 200)), "test")
    assert "inf" not in report


def test_is_annotated_matches_on_the_intron_interval():
    j = parse_junction_id("chr22:101:200:+")
    assert is_annotated(j, frozenset({("chr22", 100, 200, "+")}))
    assert not is_annotated(j, frozenset({("chr22", 100, 200, "-")}))
