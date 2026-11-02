from app.dag_component import EDGES, NODES, render_dag


class FakeComponents:
    def __init__(self):
        self.document = None
        self.options = None

    def html(self, document, **options):
        self.document = document
        self.options = options


def test_graph_contract_has_valid_endpoints_and_semantic_layers():
    ids = {node["id"] for node in NODES}
    assert len(ids) == len(NODES)
    assert all(source in ids and target in ids for source, target in EDGES)
    assert {node["layer"] for node in NODES} == {"input", "intermediate", "output"}
    assert any(node["id"] == "report" and node["layer"] == "output" for node in NODES)


def test_offline_component_has_no_remote_assets_or_arrowheads():
    fake = FakeComponents()
    render_dag(fake)
    assert "http://" not in fake.document
    assert "https://" not in fake.document
    assert "marker-end" not in fake.document
    assert "data-source" in fake.document
    assert "data-target" in fake.document
    assert fake.options == {"height": 510, "scrolling": False}
