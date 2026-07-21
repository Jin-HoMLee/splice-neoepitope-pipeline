"""Caller adapters: each parses one tool's native output into CommonRecords.

Adding a caller = drop a module here exposing a ``parse_<caller>(path, ...)``
function that returns ``list[CommonRecord]``, then add one line registering it
in the ``ADAPTERS`` dict in ``collect.py`` (#966, AC-5 extension point; see the
experiment README's "Adding a caller").
"""
