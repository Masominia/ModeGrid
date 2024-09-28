import json

import pytest

from modegrid import cases, engine
from modegrid.engine import ROOT

ALL = cases.discover(ROOT / "cases")


def test_library_is_present_and_described():
    assert len(ALL) >= 90
    for path in ALL:
        assert len(cases.notes(path)) >= 2, f"{path.name} has no description"
        assert "expect " in path.read_text(), f"{path.name} checks nothing"


def test_families():
    fam = {"strips": ["strip_"], "grid": ["conv_", "cut_"]}
    assert cases.family_of("strip_te_w450", fam) == "strips"
    assert cases.family_of("cut_tm_d125", fam) == "grid"
    assert cases.family_of("slot_gap60", fam) == "other"


@pytest.mark.skipif(not engine.available(), reason="C++ solver not built")
@pytest.mark.parametrize("path", ALL, ids=lambda p: p.stem)
def test_every_structure_meets_its_expectations(path):
    result = engine.run(path, "check")
    assert result["checks"], "no expectations evaluated"
    failed = [c for c in result["checks"] if c["status"] != "PASS"]
    assert not failed, failed
    assert result["metrics"]["converged"] == 1


@pytest.mark.skipif(not engine.available(), reason="C++ solver not built")
def test_record_round_trip(tmp_path):
    record, result = cases.run_case(ROOT / "cases" / "strip_te_w450.wg")
    back = json.loads(cases.write_record(record, tmp_path).read_text())
    assert back["case"] == "strip_te_w450" and 2.3 < back["metrics"]["neff"] < 2.5
    assert "field" not in back and result["field"].shape == (111, 60) == result["eps"].shape
