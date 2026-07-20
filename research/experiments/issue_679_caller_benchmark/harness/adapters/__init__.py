"""Caller adapters: each parses one tool's native output into CommonRecords.

Adding a caller = drop a module here exposing a ``parse_<caller>(path, ...)``
function that returns ``list[CommonRecord]``, then register it in
``callers.yaml`` (#966, AC-5 extension point).
"""
