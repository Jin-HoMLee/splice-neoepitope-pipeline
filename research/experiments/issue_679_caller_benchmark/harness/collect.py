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

from adapters.splice2neo import parse_splice2neo
from common_schema import CommonRecord

# The caller registry: caller name -> a parse function returning CommonRecords.
# One line per caller is the whole extension point (AC-5).
ADAPTERS = {
    "splice2neo": parse_splice2neo,
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

    Returns one flat ``list[CommonRecord]`` across all inputs.
    """
    records = []
    for caller, path in specs:
        records.extend(get_adapter(caller)(path))
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
        required=True,
        type=_parse_input_spec,
        metavar="CALLER:PATH",
        help="a caller name and its output path (repeatable)",
    )
    parser.add_argument("--out", required=True, help="unified TSV output path")
    args = parser.parse_args(argv)

    records = collect(args.inputs)
    records_to_tsv(records, args.out)
    print(f"collected {len(records)} records from {len(args.inputs)} input(s) -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
