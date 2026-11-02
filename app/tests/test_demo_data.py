import json
from pathlib import Path

import pytest

from app.demo_data import DEMO_PATH, DemoFixtureError, load_demo


def test_bundled_demo_is_explicitly_synthetic_and_offline():
    demo = load_demo()
    metadata = demo["metadata"]
    assert metadata["demo_fixture"] is True
    assert metadata["synthetic"] is True
    assert metadata["network_required"] is False
    assert metadata["clinical_use"] is False
    assert "not a MHCflurry prediction" in metadata["prediction_status"]
    assert demo["funnel"][0]["count"] == 124


def test_unsafe_demo_metadata_is_rejected(tmp_path: Path):
    payload = json.loads(DEMO_PATH.read_text())
    payload["metadata"]["synthetic"] = False
    path = tmp_path / "demo.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(DemoFixtureError, match="synthetic"):
        load_demo(path)
