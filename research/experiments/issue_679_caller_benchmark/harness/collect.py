"""Harness collect entry point for the open splice-caller benchmark (#966).

The collect step is the reusable substrate leaf B provides to the rest of the
epic: it takes one or more caller outputs, dispatches each to its registered
adapter, merges the results into the common schema, and writes a unified TSV
that the downstream concordance leaves (detection, burden) plug into.

Extension point (AC-5): adding a caller is a one-line ``ADAPTERS`` entry
pointing at an ``adapters/<caller>.py`` parse function that returns
``list[CommonRecord]``. No change to this module's logic is needed.

Invocation:

    python collect.py --input splice2neo:path/to/splice2neo_out.tsv --out unified.tsv

``--input`` is repeatable (``caller:path``); every caller must be registered
in ``ADAPTERS``.
"""

import argparse
import dataclasses
import json
import sys
import tempfile

from adapters.asneo import parse_asneo
from adapters.splice2neo import parse_splice2neo
from common_schema import CommonRecord, validate
from runners.asneo import run_asneo

# The adapter registry: caller name -> a parse function returning CommonRecords.
# One line per caller is the whole ingest-side extension point (AC-5).
ADAPTERS = {
    "splice2neo": parse_splice2neo,
    "asneo": parse_asneo,
}

# The runner registry: caller name -> a function (sj_tab, workdir) -> output path.
# One line per caller alongside its adapter extends AC-5 to the run step.
RUNNERS = {
    "asneo": run_asneo,
}


def get_adapter(caller):
    """Return the parse function registered for ``caller``.

    Raises ``ValueError`` naming the unknown caller and the known ones so a
    typo or an unregistered adapter fails loudly at dispatch, not silently.
    """
    try:
        return ADAPTERS[caller]
    except KeyError:
        known = ", ".join(sorted(ADAPTERS))
        raise ValueError(
            f"unknown caller {caller!r}; registered callers: {known}"
        )


def collect(specs):
    """Dispatch each ``(caller, path)`` spec to its adapter and merge results.

    Every record is validated against the legal record_level / field
    combinations before it is admitted, so an adapter bug fails loudly here
    rather than writing a mislabeled row. Returns one flat
    ``list[CommonRecord]`` across all inputs.
    """
    records = []
    for caller, path in specs:
        for rec in get_adapter(caller)(path):
            validate(rec)
            records.append(rec)
    return records


def _field_names():
    return [f.name for f in dataclasses.fields(CommonRecord)]


def _cell(value):
    """Render one field as a plain tab-delimited cell.

    ``dict`` fields (``provenance``) are JSON-encoded verbatim so the column
    round-trips with a naive ``split("\\t")`` - the format stays greppable and
    ``pandas.read_csv(sep="\\t")``-friendly rather than csv-quoted. A tab inside
    a value would corrupt the columns, so it is rejected loudly.
    """
    if value is None:
        return ""
    if isinstance(value, dict):
        value = json.dumps(value, sort_keys=True)
    text = str(value)
    if "\t" in text or "\n" in text:
        raise ValueError(f"field value contains a tab/newline: {text!r}")
    return text


def records_to_tsv(records, path):
    """Write records to a unified, plainly tab-delimited TSV.

    Columns are the CommonRecord fields in declaration order; ``provenance`` is
    JSON-encoded (see ``_cell``).
    """
    fields = _field_names()
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(fields) + "\n")
        for r in records:
            fh.write("\t".join(_cell(getattr(r, name)) for name in fields) + "\n")


def _parse_input_spec(spec):
    """Split a ``caller:path`` CLI token, tolerating colons in the path."""
    caller, sep, path = spec.partition(":")
    if not sep or not caller or not path:
        raise argparse.ArgumentTypeError(
            f"--input must be caller:path, got {spec!r}"
        )
    return caller, path


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Collect open splice-caller outputs into the common schema."
    )
    parser.add_argument(
        "--input",
        dest="inputs",
        action="append",
        default=[],
        type=_parse_input_spec,
        metavar="CALLER:PATH",
        help="a caller name and its pre-run output path to ingest (repeatable)",
    )
    parser.add_argument(
        "--run",
        dest="runs",
        action="append",
        default=[],
        type=_parse_input_spec,
        metavar="CALLER:SJ_TAB",
        help="run a caller on an SJ.out.tab, then ingest its output (repeatable)",
    )
    parser.add_argument("--out", required=True, help="unified TSV output path")
    args = parser.parse_args(argv)

    if not args.inputs and not args.runs:
        parser.error("provide at least one --input or --run")

    specs = list(args.inputs)
    for caller, sj_tab in args.runs:
        try:
            runner = RUNNERS[caller]
        except KeyError:
            known = ", ".join(sorted(RUNNERS))
            raise ValueError(f"no runner for {caller!r}; registered runners: {known}")
        out_path = runner(sj_tab, tempfile.mkdtemp(prefix=f"{caller}_run_"))
        specs.append((caller, out_path))

    records = collect(specs)
    records_to_tsv(records, args.out)
    print(f"collected {len(records)} records from {len(specs)} input(s) -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
