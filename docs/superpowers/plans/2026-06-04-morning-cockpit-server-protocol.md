# Morning Cockpit — Server + Protocol (Phase 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the purpose-built local server + the push/pull protocol that backs the visual morning-routine cockpit — token-gated, server-stamped event sequencing, atomic state writes — fully unit-tested in CI, with no browser, no `gh`, and no skill yet.

**Architecture:** A stdlib `http.server` (no third-party deps) serves a `state.json` the session writes atomically and accepts token-authenticated `POST /event`s that it stamps with a monotonic server-side `seq` and appends to `events.jsonl`. A re-arming Bash watcher (driven by the session in later phases) exits when `seq > cursor`. This phase delivers the server, the `EventStore`/runtime helpers, their tests, the CI wiring, and the Phase-0 real-browser smoke harness.

**Tech Stack:** Python 3.11 stdlib only (`http.server`, `json`, `secrets`, `threading`, `tempfile`); pytest (CI-tools convention); Bash for the watcher/smoke harness.

**Spec:** [`docs/superpowers/specs/2026-06-04-visual-morning-routine-design.md`](../specs/2026-06-04-visual-morning-routine-design.md) — this plan implements §3, §4 (schemas), §8.1 (watcher), §10 (idle watchdog), §12 (tests), §14 (security controls), and §18 build Phase 0 + Phase 1. Tracking: [Issue #655](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/655).

---

## File Structure

All code co-locates with its tests in `tools/morning/`, mirroring `tools/news/` and `tools/ci/` (flat `test_*.py`, per-dir `pytest.ini`):

| File | Responsibility |
|---|---|
| `tools/morning/runtime.py` | Atomic JSON/text writes (`temp + os.replace`), tolerant reads. Pure functions, no I/O policy. |
| `tools/morning/events.py` | `EventStore` — server-side monotonic `seq` assignment, append, `read_since(cursor)` with partial-line safety, `truncate`. |
| `tools/morning/server.py` | `MorningServer` — HTTP handler (`/`, `/health`, `/state.json`, `POST /event`), token/ctype/size gates, log redaction, idle watchdog, `launch()`. |
| `tools/morning/cockpit_min.html` | Minimal token-embedded button page for the Phase-0 smoke (the real cockpit is Phase 2). |
| `tools/morning/watcher.sh` | The re-arming watcher (exits on `seq > cursor`); used by the smoke harness now, by the engine in Phase 3. |
| `tools/morning/test_runtime.py` · `test_events.py` · `test_server.py` | Unit + protocol tests (CI-safe, no network beyond loopback). |
| `tools/morning/conftest.py` | Adds the dir to `sys.path` for flat intra-module imports. |
| `tools/morning/pytest.ini` | `live` marker registration (mirror `tools/news/pytest.ini`). |
| `tools/morning/README.md` | Interpreter (CI Python 3.11), how to run tests, the Phase-0 smoke procedure. |
| `.github/workflows/tests.yml` | Add `pytest tools/morning/ -v` to the `ci-tools-pytest` run block. |

---

### Task 1: Scaffold `tools/morning/` + wire into CI

**Files:**
- Create: `tools/morning/conftest.py`, `tools/morning/pytest.ini`, `tools/morning/test_smoke.py`
- Modify: `.github/workflows/tests.yml:77-81`

- [ ] **Step 1: Create `conftest.py`** (puts the dir on `sys.path` so `import server` works in flat layout)

```python
# tools/morning/conftest.py
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
```

- [ ] **Step 2: Create `pytest.ini`** (mirror `tools/news/pytest.ini`)

```ini
[pytest]
markers =
    live: tests requiring a real browser / harness (opt-in via -m live)
```

- [ ] **Step 3: Create a trivial collection-sanity test**

```python
# tools/morning/test_smoke.py
def test_collection_works():
    assert True
```

- [ ] **Step 4: Run it to confirm the dir collects**

Run: `python -m pytest tools/morning/ -v`
Expected: PASS (1 passed)

- [ ] **Step 5: Add `tools/morning/` to the `ci-tools-pytest` job**

In `.github/workflows/tests.yml`, change the `Run CI-tooling tests` run block from:

```yaml
        run: |
          pytest tools/ci/ -v
          pytest tools/news/ -v
```

to:

```yaml
        run: |
          pytest tools/ci/ -v
          pytest tools/news/ -v
          pytest tools/morning/ -v
```

- [ ] **Step 6: Commit**

```bash
git add tools/morning/conftest.py tools/morning/pytest.ini tools/morning/test_smoke.py .github/workflows/tests.yml
git commit -m "feat(morning): scaffold tools/morning/ + wire into ci-tools-pytest"
```

---

### Task 2: `runtime.py` — atomic writes + tolerant reads

**Files:**
- Create: `tools/morning/runtime.py`, `tools/morning/test_runtime.py`

- [ ] **Step 1: Write the failing tests**

```python
# tools/morning/test_runtime.py
import json
import os
import runtime

def test_atomic_write_json_roundtrips(tmp_path):
    p = str(tmp_path / "state.json")
    runtime.atomic_write_json(p, {"a": 1, "b": [2, 3]})
    with open(p) as f:
        assert json.load(f) == {"a": 1, "b": [2, 3]}

def test_atomic_write_leaves_no_tmp_files(tmp_path):
    p = str(tmp_path / "state.json")
    runtime.atomic_write_json(p, {"x": 1})
    leftovers = [n for n in os.listdir(tmp_path) if n != "state.json"]
    assert leftovers == []

def test_atomic_write_overwrites_in_place(tmp_path):
    p = str(tmp_path / "state.json")
    runtime.atomic_write_json(p, {"v": 1})
    runtime.atomic_write_json(p, {"v": 2})
    assert runtime.read_json(p) == {"v": 2}

def test_read_json_missing_returns_default(tmp_path):
    assert runtime.read_json(str(tmp_path / "nope.json"), {"boot": True}) == {"boot": True}

def test_read_json_corrupt_returns_default(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text("{not json")
    assert runtime.read_json(str(p), {"boot": True}) == {"boot": True}

def test_atomic_write_text_roundtrips(tmp_path):
    p = str(tmp_path / "port")
    runtime.atomic_write_text(p, "50824")
    with open(p) as f:
        assert f.read() == "50824"
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tools/morning/test_runtime.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'runtime'`

- [ ] **Step 3: Implement `runtime.py`**

```python
# tools/morning/runtime.py
"""Atomic file writes (temp + os.replace) and tolerant reads.

state.json is polled ~1s by the cockpit while the session rewrites it; a
non-atomic write guarantees torn reads (spec §4.1, review finding). os.replace
is atomic on POSIX same-filesystem.
"""
import json
import os
import tempfile


def _atomic_write(path, data):
    d = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".tmp-morning-")
    try:
        with os.fdopen(fd, "w") as f:
            f.write(data)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def atomic_write_json(path, obj):
    _atomic_write(path, json.dumps(obj))


def atomic_write_text(path, text):
    _atomic_write(path, str(text))


def read_json(path, default=None):
    try:
        with open(path) as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return default
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tools/morning/test_runtime.py -v`
Expected: PASS (6 passed)

- [ ] **Step 5: Commit**

```bash
git add tools/morning/runtime.py tools/morning/test_runtime.py
git commit -m "feat(morning): atomic json/text writes + tolerant reads (runtime.py)"
```

---

### Task 3: `events.py` — `EventStore` append with server-assigned `seq`

**Files:**
- Create: `tools/morning/events.py`, `tools/morning/test_events.py`

- [ ] **Step 1: Write the failing tests**

```python
# tools/morning/test_events.py
import json
import events

def test_append_assigns_monotonic_seq(tmp_path):
    store = events.EventStore(str(tmp_path / "events.jsonl"))
    assert store.append({"event_kind": "action", "id": "a"}, now_ms=1) == 1
    assert store.append({"event_kind": "nav", "id": "next"}, now_ms=2) == 2
    assert store.max_seq() == 2

def test_append_stamps_seq_and_recv_at_on_disk(tmp_path):
    p = tmp_path / "events.jsonl"
    store = events.EventStore(str(p))
    store.append({"event_kind": "action", "id": "x"}, now_ms=123)
    rec = json.loads(p.read_text().splitlines()[0])
    assert rec["seq"] == 1 and rec["recv_at"] == 123 and rec["id"] == "x"

def test_seq_continues_from_existing_file(tmp_path):
    p = tmp_path / "events.jsonl"
    p.write_text(json.dumps({"seq": 5, "id": "old"}) + "\n")
    store = events.EventStore(str(p))
    assert store.append({"id": "new"}, now_ms=1) == 6

def test_truncate_resets_seq_and_file(tmp_path):
    p = tmp_path / "events.jsonl"
    store = events.EventStore(str(p))
    store.append({"id": "a"}, now_ms=1)
    store.truncate()
    assert store.max_seq() == 0
    assert p.read_text() == ""
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tools/morning/test_events.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'events'`

- [ ] **Step 3: Implement the append/seq half of `events.py`**

```python
# tools/morning/events.py
"""Server-mediated event queue.

The server assigns a strictly-increasing seq to every event (spec §4.2); seq —
NOT the browser wall-clock — is the only dedup/progress key. Single-process
appends under a lock; readers ignore a trailing partial line.
"""
import json
import threading
import time


class EventStore:
    def __init__(self, path):
        self.path = path
        self._lock = threading.Lock()
        self._seq = self._scan_max_seq()

    def _scan_max_seq(self):
        mx = 0
        try:
            with open(self.path) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        mx = max(mx, int(json.loads(line).get("seq", 0)))
                    except (json.JSONDecodeError, ValueError, TypeError):
                        continue
        except FileNotFoundError:
            pass
        return mx

    def append(self, event, now_ms=None):
        with self._lock:
            self._seq += 1
            rec = dict(event)
            rec["seq"] = self._seq
            rec["recv_at"] = now_ms if now_ms is not None else int(time.time() * 1000)
            with open(self.path, "a") as f:
                f.write(json.dumps(rec) + "\n")
            return self._seq

    def max_seq(self):
        with self._lock:
            return self._seq

    def truncate(self):
        with self._lock:
            open(self.path, "w").close()
            self._seq = 0
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tools/morning/test_events.py -v`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add tools/morning/events.py tools/morning/test_events.py
git commit -m "feat(morning): EventStore append with server-assigned seq (events.py)"
```

---

### Task 4: `events.py` — `read_since(cursor)` with partial-line safety

**Files:**
- Modify: `tools/morning/events.py`, `tools/morning/test_events.py`

- [ ] **Step 1: Add the failing tests**

```python
# append to tools/morning/test_events.py
def test_read_since_returns_only_newer(tmp_path):
    store = events.EventStore(str(tmp_path / "events.jsonl"))
    store.append({"id": "a"}, now_ms=1)
    store.append({"id": "b"}, now_ms=2)
    store.append({"id": "c"}, now_ms=3)
    seqs = [r["seq"] for r in store.read_since(1)]
    assert seqs == [2, 3]

def test_read_since_ignores_trailing_partial_line(tmp_path):
    p = tmp_path / "events.jsonl"
    p.write_text(json.dumps({"seq": 1, "id": "a"}) + "\n" + '{"seq": 2, "id": "par')  # no newline
    store = events.EventStore(str(p))
    seqs = [r["seq"] for r in store.read_since(0)]
    assert seqs == [1]  # the half-written line is skipped

def test_read_since_skips_garbled_lines(tmp_path):
    p = tmp_path / "events.jsonl"
    p.write_text(json.dumps({"seq": 1, "id": "a"}) + "\n" + "GARBAGE\n" + json.dumps({"seq": 2, "id": "b"}) + "\n")
    store = events.EventStore(str(p))
    seqs = [r["seq"] for r in store.read_since(0)]
    assert seqs == [1, 2]

def test_read_since_missing_file_returns_empty(tmp_path):
    store = events.EventStore(str(tmp_path / "nope.jsonl"))
    assert store.read_since(0) == []
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tools/morning/test_events.py -k read_since -v`
Expected: FAIL with `AttributeError: 'EventStore' object has no attribute 'read_since'`

- [ ] **Step 3: Add `read_since` to `EventStore`**

```python
# add method to EventStore in tools/morning/events.py
    def read_since(self, cursor):
        try:
            with open(self.path) as f:
                data = f.read()
        except FileNotFoundError:
            return []
        nl = data.rfind("\n")
        if nl == -1:
            return []  # only a partial line, nothing complete yet
        out = []
        for line in data[:nl].split("\n"):
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue  # garbled line advances past, never crashes (spec §14)
            if int(rec.get("seq", 0)) > cursor:
                out.append(rec)
        return out
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tools/morning/test_events.py -v`
Expected: PASS (8 passed)

- [ ] **Step 5: Commit**

```bash
git add tools/morning/events.py tools/morning/test_events.py
git commit -m "feat(morning): EventStore.read_since with partial-line + garble safety"
```

---

### Task 5: `server.py` — `MorningServer` skeleton, `/health` + `/state.json`

**Files:**
- Create: `tools/morning/server.py`, `tools/morning/test_server.py`

- [ ] **Step 1: Write the failing tests** (a fixture starts a real server on an ephemeral port)

```python
# tools/morning/test_server.py
import json
import threading
import urllib.request
import urllib.error
import pytest
import server as server_mod


@pytest.fixture
def live_server(tmp_path):
    srv = server_mod.MorningServer(str(tmp_path), idle_timeout=9999)
    httpd = srv.start(port=0)
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    base = "http://127.0.0.1:%d" % httpd.server_address[1]
    yield srv, base
    httpd.shutdown()


def _get(base, path):
    with urllib.request.urlopen(base + path, timeout=5) as r:
        return r.status, r.read().decode()


def test_health_returns_nonce_and_seq(live_server):
    srv, base = live_server
    status, body = _get(base, "/health")
    assert status == 200
    data = json.loads(body)
    assert data["nonce"] == srv.nonce and data["seq"] == 0


def test_state_json_serves_boot_when_absent(live_server):
    srv, base = live_server
    status, body = _get(base, "/state.json")
    assert status == 200 and json.loads(body) == {"boot": True}


def test_state_json_serves_written_state(live_server, tmp_path):
    srv, base = live_server
    import runtime
    runtime.atomic_write_json(srv.state_path, {"state_seq": 7, "role": "pm"})
    status, body = _get(base, "/state.json")
    assert json.loads(body)["state_seq"] == 7
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tools/morning/test_server.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'server'`

- [ ] **Step 3: Implement the server skeleton**

```python
# tools/morning/server.py
"""Purpose-built morning-cockpit server (stdlib only).

Push: serves state.json (session writes atomically). Pull: token-gated
POST /event -> EventStore (server-stamped seq). See spec §3, §4, §14.
"""
import json
import os
import secrets
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import runtime
from events import EventStore

MAX_BODY = 4096


class MorningServer:
    def __init__(self, runtime_dir, idle_timeout=1800, idle_check=None, token=None, nonce=None):
        self.dir = runtime_dir
        os.makedirs(runtime_dir, exist_ok=True)
        self.events = EventStore(os.path.join(runtime_dir, "events.jsonl"))
        self.state_path = os.path.join(runtime_dir, "state.json")
        self.token = token or secrets.token_hex(16)
        self.nonce = nonce or secrets.token_hex(8)
        self.idle_timeout = idle_timeout
        self.idle_check = idle_check if idle_check is not None else min(30, idle_timeout)
        self._last_activity = time.monotonic()
        self._httpd = None

    def touch(self):
        self._last_activity = time.monotonic()

    def page(self):
        path = os.path.join(os.path.dirname(__file__), "cockpit_min.html")
        with open(path) as f:
            return f.read().replace("__TOKEN__", self.token)

    def _handler_class(self):
        srv = self

        class H(BaseHTTPRequestHandler):
            def log_message(self, fmt, *args):
                try:
                    status = args[1] if len(args) > 1 else "?"
                    sys.stderr.write("%s %s %s\n" % (self.command, self.path.split("?")[0], status))
                except Exception:
                    pass

            def _send(self, code, body, ctype="application/json"):
                b = body.encode() if isinstance(body, str) else body
                self.send_response(code)
                self.send_header("Content-Type", ctype)
                self.send_header("Content-Length", str(len(b)))
                self.end_headers()
                self.wfile.write(b)

            def do_GET(self):
                srv.touch()
                p = self.path.split("?")[0]
                if p == "/":
                    return self._send(200, srv.page(), "text/html")
                if p == "/health":
                    return self._send(200, json.dumps({"nonce": srv.nonce, "seq": srv.events.max_seq()}))
                if p == "/state.json":
                    return self._send(200, json.dumps(runtime.read_json(srv.state_path, {"boot": True})))
                return self._send(404, "{}")

        return H

    def start(self, port=0):
        self._httpd = ThreadingHTTPServer(("127.0.0.1", port), self._handler_class())
        return self._httpd
```

Note: `page()` reads `cockpit_min.html`; create a stub now so `/` doesn't crash if hit (full page in Task 6).

- [ ] **Step 4: Create the minimal page stub**

```html
<!-- tools/morning/cockpit_min.html -->
<!doctype html><html><head><meta charset="utf-8"><title>morning</title></head>
<body><button onclick="fire()">POST /event</button>
<script>
const TOKEN="__TOKEN__";
async function fire(){
  await fetch('/event',{method:'POST',headers:{'Content-Type':'application/json','X-Morning-Token':TOKEN},body:JSON.stringify({event_kind:'action',id:'manual_click'})});
}
</script></body></html>
```

- [ ] **Step 5: Run to verify pass**

Run: `python -m pytest tools/morning/test_server.py -v`
Expected: PASS (3 passed)

- [ ] **Step 6: Commit**

```bash
git add tools/morning/server.py tools/morning/cockpit_min.html tools/morning/test_server.py
git commit -m "feat(morning): MorningServer skeleton + /health + /state.json"
```

---

### Task 6: `server.py` — `POST /event` happy path

**Files:**
- Modify: `tools/morning/server.py`, `tools/morning/test_server.py`

- [ ] **Step 1: Add a POST helper + the failing test**

```python
# append to tools/morning/test_server.py
def _post(base, path, body, headers):
    req = urllib.request.Request(base + path, data=body.encode(), method="POST")
    for k, v in headers.items():
        req.add_header(k, v)
    try:
        with urllib.request.urlopen(req, timeout=5) as r:
            return r.status, r.read().decode()
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode()


def test_post_event_appends_and_returns_seq(live_server):
    srv, base = live_server
    status, body = _post(base, "/event", '{"event_kind":"action","id":"commit_644"}',
                         {"Content-Type": "application/json", "X-Morning-Token": srv.token})
    assert status == 200 and json.loads(body) == {"ok": True, "seq": 1}
    assert srv.events.max_seq() == 1
    recs = srv.events.read_since(0)
    assert recs[0]["id"] == "commit_644" and recs[0]["event_kind"] == "action"
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tools/morning/test_server.py::test_post_event_appends_and_returns_seq -v`
Expected: FAIL — server returns 501/404 (no `do_POST` yet)

- [ ] **Step 3: Add `do_POST` to the handler class** (inside `_handler_class`, after `do_GET`)

```python
            def do_POST(self):
                srv.touch()
                if self.path.split("?")[0] != "/event":
                    return self._send(404, "{}")
                if self.headers.get("X-Morning-Token") != srv.token:
                    return self._send(403, json.dumps({"error": "token"}))
                if self.headers.get("Content-Type", "").split(";")[0].strip() != "application/json":
                    return self._send(415, json.dumps({"error": "ctype"}))
                n = int(self.headers.get("Content-Length", 0) or 0)
                if n > MAX_BODY:
                    return self._send(413, json.dumps({"error": "size"}))
                try:
                    ev = json.loads(self.rfile.read(n))
                except Exception:
                    return self._send(400, json.dumps({"error": "json"}))
                if not isinstance(ev, dict):
                    return self._send(400, json.dumps({"error": "shape"}))
                seq = srv.events.append({"event_kind": ev.get("event_kind"), "id": ev.get("id")})
                return self._send(200, json.dumps({"ok": True, "seq": seq}))
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tools/morning/test_server.py -v`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add tools/morning/server.py tools/morning/test_server.py
git commit -m "feat(morning): POST /event happy path (server-stamped seq)"
```

---

### Task 7: `server.py` — security gates (token / ctype / size; never append on reject)

**Files:**
- Modify: `tools/morning/test_server.py`

The gate logic already exists from Task 6; this task locks it with adversarial tests (spec §14 controls 1–2, 4).

- [ ] **Step 1: Add the failing/locking tests**

```python
# append to tools/morning/test_server.py
def test_post_no_token_is_403_and_no_append(live_server):
    srv, base = live_server
    status, _ = _post(base, "/event", '{"id":"x"}', {"Content-Type": "application/json"})
    assert status == 403 and srv.events.max_seq() == 0

def test_post_wrong_token_is_403(live_server):
    srv, base = live_server
    status, _ = _post(base, "/event", '{"id":"x"}',
                     {"Content-Type": "application/json", "X-Morning-Token": "wrong"})
    assert status == 403 and srv.events.max_seq() == 0

def test_post_bad_ctype_is_415_and_no_append(live_server):
    srv, base = live_server
    status, _ = _post(base, "/event", "x",
                     {"Content-Type": "text/plain", "X-Morning-Token": srv.token})
    assert status == 415 and srv.events.max_seq() == 0

def test_post_oversize_is_413_and_no_append(live_server):
    srv, base = live_server
    big = '{"id":"' + "z" * 5000 + '"}'
    status, _ = _post(base, "/event", big,
                     {"Content-Type": "application/json", "X-Morning-Token": srv.token})
    assert status == 413 and srv.events.max_seq() == 0
```

- [ ] **Step 2: Run to verify pass** (gates implemented in Task 6 — these lock them)

Run: `python -m pytest tools/morning/test_server.py -k "403 or 415 or 413 or token" -v`
Expected: PASS (4 passed)

- [ ] **Step 3: Commit**

```bash
git add tools/morning/test_server.py
git commit -m "test(morning): lock POST /event security gates (token/ctype/size)"
```

---

### Task 8: `server.py` — log redaction

**Files:**
- Modify: `tools/morning/test_server.py` (the `log_message` override exists from Task 5)

- [ ] **Step 1: Add the failing/locking test** (capture stderr, assert no body/query/token leaks)

```python
# append to tools/morning/test_server.py
def test_log_redacts_query_and_body(live_server, capfd):
    srv, base = live_server
    _post(base, "/event?token=SECRETLEAK", '{"id":"hush"}',
          {"Content-Type": "application/json", "X-Morning-Token": srv.token})
    _get(base, "/health")
    err = capfd.readouterr().err
    assert "SECRETLEAK" not in err   # query string never logged
    assert "hush" not in err          # body never logged
    assert srv.token not in err       # token never logged
    assert "POST /event" in err       # method + path-without-query IS logged
```

- [ ] **Step 2: Run to verify pass** (redaction implemented in Task 5's `log_message`)

Run: `python -m pytest tools/morning/test_server.py::test_log_redacts_query_and_body -v`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add tools/morning/test_server.py
git commit -m "test(morning): lock log redaction (no query/body/token in logs)"
```

---

### Task 9: `server.py` — idle watchdog + `launch()` (server-info, nonce, port, atomic)

**Files:**
- Modify: `tools/morning/server.py`, `tools/morning/test_server.py`

- [ ] **Step 1: Write the failing tests**

```python
# append to tools/morning/test_server.py
import time

def test_idle_watchdog_shuts_down_when_idle(tmp_path):
    srv = server_mod.MorningServer(str(tmp_path), idle_timeout=0.5, idle_check=0.1)
    httpd = srv.start(port=0)
    t = threading.Thread(target=httpd.serve_forever, daemon=True)
    t.start()
    srv.start_idle_watchdog()
    t.join(timeout=3.0)
    assert not t.is_alive()  # serve_forever returned => watchdog shut it down

def test_launch_writes_server_info_atomically(tmp_path):
    srv = server_mod.MorningServer(str(tmp_path), idle_timeout=9999)
    httpd = srv.start(port=0)
    srv.write_server_info(httpd.server_address[1])
    import runtime
    info = runtime.read_json(str(tmp_path / "server-info"))
    assert info["token"] == srv.token and info["nonce"] == srv.nonce
    assert info["pid"] and info["port"] == httpd.server_address[1]
    with open(tmp_path / "port") as f:
        assert f.read() == str(httpd.server_address[1])
    httpd.shutdown()
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tools/morning/test_server.py -k "watchdog or server_info" -v`
Expected: FAIL — `AttributeError` (`start_idle_watchdog` / `write_server_info` missing)

- [ ] **Step 3: Add the watchdog + launch helpers to `MorningServer`**

```python
# add methods to MorningServer in tools/morning/server.py
    def start_idle_watchdog(self):
        def loop():
            while True:
                time.sleep(self.idle_check)
                if self._httpd is None:
                    return
                if time.monotonic() - self._last_activity > self.idle_timeout:
                    self._httpd.shutdown()
                    return
        threading.Thread(target=loop, daemon=True).start()

    def write_server_info(self, port):
        runtime.atomic_write_json(os.path.join(self.dir, "server-info"),
                                  {"token": self.token, "nonce": self.nonce,
                                   "pid": os.getpid(), "port": port})
        runtime.atomic_write_text(os.path.join(self.dir, "port"), str(port))

    def launch(self, port=0):
        httpd = self.start(port)
        self.write_server_info(httpd.server_address[1])
        self.start_idle_watchdog()
        sys.stderr.write("morning server listening on %d\n" % httpd.server_address[1])
        return httpd
```

- [ ] **Step 4: Add the CLI entrypoint at the bottom of `server.py`**

```python
# bottom of tools/morning/server.py
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("runtime_dir")
    ap.add_argument("--port", type=int, default=0)
    ap.add_argument("--idle-timeout", type=int, default=1800)
    a = ap.parse_args()
    s = MorningServer(a.runtime_dir, idle_timeout=a.idle_timeout)
    httpd = s.launch(a.port)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        httpd.shutdown()
```

- [ ] **Step 5: Run to verify pass**

Run: `python -m pytest tools/morning/test_server.py -v`
Expected: PASS (all)

- [ ] **Step 6: Commit**

```bash
git add tools/morning/server.py tools/morning/test_server.py
git commit -m "feat(morning): idle watchdog + launch (server-info/nonce/port, atomic)"
```

---

### Task 10: Phase-0 smoke harness — re-arming watcher + real-browser re-verify

This productionises the Phase-0 spike (§18 Phase 0) against the real server, with a real browser (the spike used `curl`). Not CI — a documented manual smoke.

**Files:**
- Create: `tools/morning/watcher.sh`, `tools/morning/README.md`

- [ ] **Step 1: Create the re-arming watcher script**

```bash
# tools/morning/watcher.sh
# Exits 0 when events.jsonl has a seq greater than $1 (the cursor); else times out.
# The session re-arms this each wake (spec §8.1). seq comes from the server, not the clock.
set -u
DIR="${MORNING_DIR:?set MORNING_DIR}"
CURSOR="${1:-0}"
EV="$DIR/events.jsonl"
for _ in $(seq 1 1800); do   # 1800 * 0.5s = 15 min, then re-arm
  MAX=$(grep -o '"seq": *[0-9]*' "$EV" 2>/dev/null | grep -o '[0-9]*' | sort -n | tail -1)
  MAX=${MAX:-0}
  if [ "$MAX" -gt "$CURSOR" ]; then
    echo "WAKE maxseq=$MAX cursor=$CURSOR"
    exit 0
  fi
  sleep 0.5
done
echo "TIMEOUT cursor=$CURSOR"
exit 1
```

- [ ] **Step 2: Create the README with the smoke procedure**

````markdown
# tools/morning — cockpit server + protocol

Interpreter: **Python 3.11** (CI `ci-tools-pytest`). Stdlib only; no third-party deps.

## Tests
```bash
python -m pytest tools/morning/ -v
```

## Phase-0 real-browser smoke (manual — not CI)
Re-verifies the live bridge with a *real browser* (the design spike used a curl driver).

```bash
export MORNING_DIR=$(mktemp -d)
python tools/morning/server.py "$MORNING_DIR" --idle-timeout 1800 &   # note the printed port
# open http://127.0.0.1:<port>/ in a browser
```
Then, from a Claude Code session, drive ≥5 cycles:
1. `bash tools/morning/watcher.sh <cursor>` with `MORNING_DIR` exported, `run_in_background: true`
2. Click the page button once → the watcher exits → the session is re-invoked
3. Advance the cursor to the new max seq, re-arm (step 1), repeat
Assert: every click wakes the session; `server-info`/`port` present; clean teardown (`kill %1`) frees the port.
````

- [ ] **Step 3: Make the watcher executable + verify the full stack runs once headless**

```bash
chmod +x tools/morning/watcher.sh
export MORNING_DIR=$(mktemp -d)
python tools/morning/server.py "$MORNING_DIR" --port 0 &
SRVPID=$!; sleep 1; PORT=$(cat "$MORNING_DIR/port")
TOKEN=$(python3 -c "import json;print(json.load(open('$MORNING_DIR/server-info'))['token'])")
( MORNING_DIR="$MORNING_DIR" bash tools/morning/watcher.sh 0 & echo $! > /tmp/wpid )
curl -s -X POST "http://127.0.0.1:$PORT/event" -H "Content-Type: application/json" -H "X-Morning-Token: $TOKEN" -d '{"event_kind":"action","id":"smoke"}'
wait "$(cat /tmp/wpid)"   # watcher should exit 0
kill "$SRVPID"
```
Expected: the `curl` returns `{"ok": true, "seq": 1}`, the watcher prints `WAKE maxseq=1 cursor=0` and exits 0.

- [ ] **Step 4: Commit**

```bash
git add tools/morning/watcher.sh tools/morning/README.md
git commit -m "feat(morning): re-arming watcher + Phase-0 real-browser smoke harness"
```

---

## Self-Review

**Spec coverage (Phase 1 scope):**
- §4.1 atomic `state.json` writes → Task 2 (`atomic_write_json`) + Task 5 (`/state.json` serves it).
- §4.2 server-assigned `seq`, cursor dedup, partial-line/garble safety → Tasks 3–4.
- §3 server components (`/`, `/health`, `/state.json`, `POST /event`) → Tasks 5–6.
- §14 controls: token (Task 6/7), ctype (7), size (7), log redaction (8), `events.jsonl` garble-tolerance (4). _Origin-check + 0600 perms + token-regeneration-on-reuse are deferred to Phase 5 (launch/lifecycle) — noted so they aren't dropped._
- §10 idle watchdog (server-owned) → Task 9.
- §8.1 re-arming watcher (`seq > cursor`) → Task 10.
- §18 Phase 0 real-browser re-verify → Task 10 smoke harness. §12 CI wiring → Task 1.
- **Out of Phase-1 scope (later plans):** `boot_nonce`/`state_seq`/`schema_version` *producers* (the session writes these in Phase 3; the server only serves whatever JSON it's given, which these tests honor), `bridge`/heartbeat (Phase 3), the diff-render cockpit (Phase 2), idempotent-singleton 3-part liveness reuse (Phase 5).

**Placeholder scan:** none — every step carries real code/commands + expected output.

**Type consistency:** `MorningServer(runtime_dir, idle_timeout, idle_check, token, nonce)`, `.start(port)→httpd`, `.launch(port)→httpd`, `.start_idle_watchdog()`, `.write_server_info(port)`, `.state_path`, `.events`, `.token`, `.nonce`, `.page()` — used consistently across Tasks 5–10. `EventStore(path)` with `.append(event, now_ms)→int`, `.max_seq()`, `.read_since(cursor)→list`, `.truncate()` — consistent across Tasks 3–10. `runtime.atomic_write_json/atomic_write_text/read_json` — consistent across Tasks 2, 5, 9.

---

## Next plans (carved as sub-issues of #655 at execution)
2. Cockpit page (`cockpit.html/js/css`) — diff-render, link allowlist, bridge-state footer.
3. Engine + re-arming watcher loop + hybrid confirm dispatch + `state.json` producer (`schema_version`/`state_seq`/`boot_nonce`/`bridge`).
4. PM role config + builders (wrap `check_ready_queue.sh`/`check_milestone_health.sh`/closure-audit).
5. `SKILL.md` + role-detection allowlist + idempotent-singleton launch + crash-safe teardown + remaining §14 (origin, 0600, token-regen).
6. Dev/Sci configs — prove corrected D6 (per-role formatters + new block types).
