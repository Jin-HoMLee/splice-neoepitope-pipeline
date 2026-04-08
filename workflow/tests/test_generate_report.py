"""Tests for generate_report.py — HTML generation helpers."""

import pandas as pd

from generate_report import _df_to_html


class TestDfToHtml:
    def test_empty_dataframe_returns_no_data_message(self):
        df = pd.DataFrame(columns=["a", "b"])
        assert _df_to_html(df) == "<p><em>No data.</em></p>"

    def test_small_dataframe_returns_html_table(self):
        df = pd.DataFrame({"sample": ["s1"], "count": [5]})
        html = _df_to_html(df)
        assert "<table" in html
        assert "s1" in html

    def test_large_dataframe_shows_truncation_message(self):
        df = pd.DataFrame({"x": range(200)})
        html = _df_to_html(df, max_rows=10)
        assert "Showing 10 of 200 rows" in html

    def test_large_dataframe_respects_max_rows(self):
        df = pd.DataFrame({"x": range(50)})
        html = _df_to_html(df, max_rows=10)
        # Only 10 data rows should be in the table (row 0–9, not 10–49)
        assert "49" not in html
