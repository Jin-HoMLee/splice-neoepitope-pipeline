import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / ".agents" / "hooks"))
import graphql_meter as gm  # noqa: E402


def test_fragment_and_probe_are_literal_strings():
    assert gm.RATE_LIMIT_FRAGMENT == "rateLimit { cost remaining }"
    assert gm.RATE_LIMIT_PROBE_QUERY == "query { rateLimit { cost remaining } }"


def test_extract_rate_limit_reads_data_ratelimit():
    resp = {"data": {"rateLimit": {"cost": 2, "remaining": 4998}}}
    assert gm.extract_rate_limit(resp) == (2, 4998)


def test_extract_rate_limit_missing_is_none_none():
    assert gm.extract_rate_limit({"data": {}}) == (None, None)
    assert gm.extract_rate_limit({}) == (None, None)
    assert gm.extract_rate_limit(None) == (None, None)
    assert gm.extract_rate_limit("not a dict") == (None, None)


def test_log_graphql_spend_appends_one_line(tmp_path, monkeypatch):
    log = tmp_path / "graphql_spend.jsonl"
    monkeypatch.setattr(gm, "SPEND_LOG_PATH", log)
    resp = {"data": {"rateLimit": {"cost": 26, "remaining": 4700}}}
    gm.log_graphql_spend("board_open_items", resp, query_name="board")
    lines = log.read_text().splitlines()
    assert len(lines) == 1
    rec = json.loads(lines[0])
    assert rec["consumer"] == "board_open_items"
    assert rec["cost"] == 26
    assert rec["remaining"] == 4700
    assert rec["query_name"] == "board"
    assert rec["ts"].endswith("Z")
    assert list(rec.keys()) == ["ts", "consumer", "cost", "remaining", "query_name"]


def test_log_graphql_probe_nulls_cost(tmp_path, monkeypatch):
    log = tmp_path / "graphql_spend.jsonl"
    monkeypatch.setattr(gm, "SPEND_LOG_PATH", log)
    probe = {"data": {"rateLimit": {"cost": 0, "remaining": 4321}}}
    gm.log_graphql_probe("post_gh_pr_create", probe, query_name="status_mutation")
    rec = json.loads(log.read_text().splitlines()[0])
    assert rec["cost"] is None
    assert rec["remaining"] == 4321


def test_log_is_fail_open_on_unwritable_path(monkeypatch):
    # A directory that cannot be created / written must not raise.
    monkeypatch.setattr(gm, "SPEND_LOG_PATH", Path("/proc/nonexistent/graphql_spend.jsonl"))
    gm.log_graphql_spend("x", {"data": {"rateLimit": {"cost": 1, "remaining": 1}}})  # no raise


def test_log_is_fail_open_on_malformed_response(tmp_path, monkeypatch):
    log = tmp_path / "graphql_spend.jsonl"
    monkeypatch.setattr(gm, "SPEND_LOG_PATH", log)
    gm.log_graphql_spend("x", {"garbage": True})  # no rateLimit; must not raise
    rec = json.loads(log.read_text().splitlines()[0])
    assert rec["cost"] is None and rec["remaining"] is None


def test_board_query_carries_ratelimit_fragment():
    import importlib.util
    root = Path(__file__).resolve().parents[2]
    spec = importlib.util.spec_from_file_location("boi", root / "scripts" / "board_open_items.py")
    boi = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(boi)
    assert "rateLimit { cost remaining }" in boi.QUERY
