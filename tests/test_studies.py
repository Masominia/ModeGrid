import numpy as np

from modegrid import studies


def test_convergence_extrapolates_a_first_order_error():
    metrics = {f"conv_te_d{h:02d}": {"dx": h / 1e3, "neff": 2.5 + 0.0004 * h} for h in (10, 12, 15, 20, 30, 40)}
    metrics["strip_te_w450"] = {"dx": 0.025, "neff": 2.4}
    t = studies.convergence(metrics, "te")
    assert len(t) == 6 and abs(t.attrs["limit"] - 2.5) < 1e-12
    assert np.allclose(t["above the extrapolated value"], 0.0004 * t["cell, nm"])
    assert studies.convergence(metrics, "tm").empty


def test_width_sweep_adds_the_effective_index_method():
    metrics = {
        f"strip_te_w{w}": {"neff": 2.0 + w / 1e4, "group_index": 4.0, "confinement": 0.6} for w in (1000, 400, 500)
    }
    t = studies.width_sweep(metrics)
    assert list(t["width, nm"]) == [400, 500, 1000]
    assert t["effective index method"].is_monotonic_increasing and "difference" in t


def test_coupler_fit_recovers_the_decay_length():
    metrics = {
        f"coupler_te_gap{g}": {"splitting": 0.2 * np.exp(-g / 110.0), "coupling_length": 1.0}
        for g in (100, 150, 200, 300)
    }
    t = studies.couplers(metrics)
    assert abs(t.attrs["decay_nm"] - 110.0) < 1e-6 and t.attrs["largest_residual"] < 1e-9
    assert studies.couplers({"coupler_te_gap100": {"splitting": 0.1, "coupling_length": 1.0}}).empty
