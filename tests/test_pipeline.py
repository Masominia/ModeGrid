import json

import pytest

from modegrid import engine
from modegrid.pipeline import run


@pytest.mark.skipif(not engine.available(), reason="C++ solver not built")
def test_fast_pipeline(tmp_path):
    m = run(out=tmp_path, fast=True)
    assert not [c for c in m["checks"] if c["status"] == "FAIL"]
    assert m["structures"] >= 25 and m["largest_differences"]["neff"] < 1e-9
    assert "width_sweep" in m["studies"] and "couplers" in m["studies"]
    assert (tmp_path / "report.md").exists() and (tmp_path / "figures" / "overview.png").exists()
    records = list((tmp_path / "cases").glob("*.json"))
    assert len(records) == m["structures"] == len(list((tmp_path / "figures" / "cases").glob("*.png")))
    assert json.loads(records[0].read_text())["against_reference"]
