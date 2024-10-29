import numpy as np
import pytest

from modegrid import engine, reference
from modegrid.compare import compare
from modegrid.engine import ROOT

HEAD = "grid dx=0.05 dy=0.04\nwindow width=2.4 ymin=-0.8 ymax=1.04\n"


def structure(tmp_path, text):
    path = tmp_path / "s.wg"
    path.write_text(text)
    return reference.read_structure(path)


def test_reader_understands_every_statement(tmp_path):
    s = structure(
        tmp_path,
        "# note\ntitle A B\nwavelength 1.31\ngrid dx=0.025 dy=0.01\nwindow width=4 ymin=-1 ymax=1.2\nbackground air\n"
        "layer sio2 -1 0\nrect 3.5 -0.25 0.25 0 0.22 # core\nsymmetry odd\nfield tm\ncoupler\nexpect neff 2 1\n",
    )
    assert (s["title"], s["wavelength"], s["dx"], s["dy"], s["width"]) == ("A B", 1.31, 0.025, 0.01, 4.0)
    assert (s["background"], s["symmetry"], s["field"], s["coupler"]) == ("air", "odd", "tm", True)
    assert s["rects"][0][0] == "sio2" and s["rects"][0][1] < -1e8
    assert s["rects"][1] == ("3.5", -0.25, 0.25, 0.0, 0.22)


def test_painting_weights_cut_cells_by_area(tmp_path):
    s = structure(tmp_path, "grid dx=1 dy=1\nwindow width=4 ymin=0 ymax=4\nbackground 1\nrect 2 -0.5 0 0.5 1\n")
    eps, x0, y0 = reference.permittivity(s, 1.55, half=False)
    assert (x0, y0, eps.shape) == (-2.0, 0.0, (4, 4))
    assert abs(eps[0, 1] - (0.75 * 1 + 0.25 * 4)) < 1e-12 and eps.sum() == pytest.approx(16 + 0.75)


def test_operator_of_a_uniform_box_has_the_analytic_spectrum(tmp_path):
    s = structure(tmp_path, "grid dx=0.05 dy=0.04\nwindow width=1.2 ymin=0 ymax=0.64\nbackground 2\n")
    eps, _, _ = reference.permittivity(s, 1.55, half=False)
    k = 2 * np.pi / 1.55
    lx, ly = (2 / 0.05 * np.sin(np.pi / 48)) ** 2, (2 / 0.04 * np.sin(np.pi / 32)) ** 2
    for field in ("te", "tm", "scalar"):
        n, e = reference.solve(eps, 0.05, 0.04, 1.55, field, "none")
        assert abs((k * n) ** 2 - (4 * k * k - lx - ly)) < 1e-8
        assert e.min() > 0


def test_te_operator_is_not_symmetric_at_an_interface(tmp_path):
    s = structure(tmp_path, HEAD + "rect si -0.25 0.25 0 0.24\n")
    eps, _, _ = reference.permittivity(s, 1.55, half=False)
    te = reference.operator(eps, 0.05, 0.04, 1.55, "te", "none")
    scalar = reference.operator(eps, 0.05, 0.04, 1.55, "scalar", "none")
    assert abs(te - te.T).max() > 1 and abs(scalar - scalar.T).max() < 1e-12


def test_slab_index_and_effective_index_method():
    n = reference.slab_index(1.55, 0.22, 3.476, 1.444, 1.444, "te")
    k = 2 * np.pi / 1.55
    kappa, gamma = k * np.sqrt(3.476**2 - n**2), k * np.sqrt(n**2 - 1.444**2)
    assert abs(np.tan(kappa * 0.11) - gamma / kappa) < 1e-9
    assert reference.slab_index(1.55, 0.01, 3.476, 1.0, 1.444, "te") is None
    eim = reference.effective_index_method(1.55, 0.45, 0.22, 3.476, 1.444, "te")
    assert 1.444 < eim < n
    assert (
        reference.effective_index_method(1.55, 0.45, 0.02, 1.5, 1.0, "te") is not None
    )  # symmetric slabs never cut off


def test_column_index_matches_the_exact_slab():
    column = np.full(260, 1.444**2)
    column[120:142] = 3.476**2  # 220 nm of silicon on a 10 nm grid
    for pol, tol in (("te", 2e-3), ("tm", 3e-3)):
        exact = reference.slab_index(1.55, 0.22, 3.476, 1.444, 1.444, pol)
        assert abs(reference.column_index(column, 0.01, 1.55, pol) - exact) < tol
    assert reference.column_index(np.full(50, 4.0), 0.01, 1.55, "te") == 2.0


@pytest.mark.skipif(not engine.available(), reason="C++ solver not built")
@pytest.mark.parametrize(
    "name", ["strip_te_w450", "strip_tm_w550", "rib_slab90", "coupler_te_gap200", "slot_gap100", "odd_te_w400"]
)
def test_solver_and_reference_agree(name):
    path = ROOT / "cases" / f"{name}.wg"
    cmp = compare(path, engine.run(path))
    assert not cmp["outside"], cmp["outside"]
    assert cmp["differences"]["neff"] < 1e-9 and cmp["field"] < 1e-4
