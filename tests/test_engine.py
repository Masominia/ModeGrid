import numpy as np
import pytest

from modegrid import engine

pytestmark = pytest.mark.skipif(not engine.available(), reason="C++ solver not built")

HEAD = "grid dx=0.05 dy=0.04\nwindow width=2.4 ymin=-0.8 ymax=1.04\n"


def write(tmp_path, text):
    path = tmp_path / "structure.wg"
    path.write_text(text)
    return path


def test_run_returns_grids(tmp_path):
    r = engine.run(write(tmp_path, "title Test\n" + HEAD + "rect si -0.25 0.25 0 0.24\nsymmetry even\n"))
    assert r["title"] == "Test"
    assert r["field"].shape == (46, 24) == r["eps"].shape
    assert abs(np.sum(r["field"] ** 2) - 1) < 1e-4 and r["field"].max() > 0
    assert abs(r["eps"].max() - 3.476**2) < 1e-4 and abs(r["eps"].min() - 1.444**2) < 1e-4
    assert r["metrics"]["guided"] == 1 and 2.3 < r["metrics"]["neff"] < 2.7


def test_check_reports_pass_and_fail(tmp_path):
    text = HEAD + "rect si -0.25 0.25 0 0.24\nsymmetry even\nexpect guided yes\nexpect neff 1.0 0.1\n"
    r = engine.run(write(tmp_path, text), "check")
    assert [c["status"] for c in r["checks"]] == ["PASS", "FAIL"]
    assert "field" not in r


def test_bad_file_raises_with_the_line(tmp_path):
    with pytest.raises(RuntimeError, match="line 2"):
        engine.run(write(tmp_path, "wavelength 1.55\nnonsense 1\n"))
    with pytest.raises(RuntimeError, match="whole number of cells"):
        engine.run(write(tmp_path, "window width=3.01 ymin=-1 ymax=1\n"))
    with pytest.raises(RuntimeError, match="cannot read"):
        engine.run(tmp_path / "missing.wg")


def test_polarisations_are_ordered(tmp_path):
    n = {}
    for field in ("te", "tm", "scalar"):
        text = HEAD + f"rect si -0.25 0.25 0 0.24\nsymmetry even\nfield {field}\n"
        n[field] = engine.run(write(tmp_path, text), "check")["metrics"]["neff"]
    assert n["tm"] < n["te"] < n["scalar"]
