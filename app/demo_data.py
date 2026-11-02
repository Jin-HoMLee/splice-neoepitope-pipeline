"""Schema-checked access to the bundled, read-only synthetic demo."""

import hashlib
import json
from pathlib import Path
from typing import Any, Dict

DEMO_PATH = Path(__file__).parent / "demo" / "demo.json"
DEMO_MANIFEST_PATH = DEMO_PATH.parent / "manifest.json"


class DemoFixtureError(ValueError):
    pass


def load_demo(path: Path = DEMO_PATH) -> Dict[str, Any]:
    raw = path.read_bytes()
    if path == DEMO_PATH:
        manifest = json.loads(DEMO_MANIFEST_PATH.read_text(encoding="utf-8"))
        expected = manifest["files"][path.name]
        observed = hashlib.sha256(raw).hexdigest()
        if observed != expected:
            raise DemoFixtureError("demo fixture checksum mismatch")
    payload = json.loads(raw)
    metadata = payload.get("metadata", {})
    required_flags = {
        "demo_fixture": True,
        "synthetic": True,
        "clinical_use": False,
        "network_required": False,
    }
    mismatches = [key for key, value in required_flags.items() if metadata.get(key) is not value]
    if mismatches:
        raise DemoFixtureError("invalid demo safety metadata: " + ", ".join(mismatches))
    for key in ("funnel", "origins", "candidates", "provenance"):
        if key not in payload:
            raise DemoFixtureError(f"demo fixture is missing {key}")
    return payload
