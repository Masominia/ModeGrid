"""Tables drawn from several structures of the library at once."""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import reference


def convergence(metrics: dict[str, dict], polarisation: str) -> pd.DataFrame:
    """Index against cell size for the aligned-grid family, with the first-order extrapolation to zero cell size."""
    rows = [
        (round(m["dx"] * 1e3), m["neff"]) for name, m in metrics.items() if name.startswith(f"conv_{polarisation}_")
    ]
    if len(rows) < 3:
        return pd.DataFrame()
    rows.sort()
    (h1, n1), (h2, n2) = rows[0], rows[1]
    limit = n1 - (n2 - n1) * h1 / (h2 - h1)
    table = pd.DataFrame(rows, columns=["cell, nm", "n_eff"])
    table["above the extrapolated value"] = table["n_eff"] - limit
    table.attrs["limit"] = float(limit)
    return table


def width_sweep(metrics: dict[str, dict], lam: float = 1.55) -> pd.DataFrame:
    """Quasi-TE index of the silicon strips against width, beside the effective index method."""
    rows = []
    for name, m in sorted(metrics.items()):
        if not name.startswith("strip_te_w"):
            continue
        width = int(name.removeprefix("strip_te_w")) / 1e3
        eim = reference.effective_index_method(lam, width, 0.22, 3.476, 1.444, "te")
        rows.append(
            {
                "width, nm": round(1e3 * width),
                "n_eff": m["neff"],
                "group index": m["group_index"],
                "confinement": m["confinement"],
                "effective index method": eim,
            }
        )
    table = pd.DataFrame(rows)
    if len(table):
        table = table.sort_values("width, nm").reset_index(drop=True)
        table["difference"] = table["effective index method"] - table["n_eff"]
    return table


def couplers(metrics: dict[str, dict], prefix: str = "coupler_te_gap") -> pd.DataFrame:
    """Coupling length against gap, and the decay length of the splitting from a straight-line fit of its logarithm."""
    rows = [
        (int(name.removeprefix(prefix)), m["splitting"], m["coupling_length"])
        for name, m in metrics.items()
        if name.startswith(prefix)
    ]
    if len(rows) < 3:
        return pd.DataFrame()
    table = pd.DataFrame(sorted(rows), columns=["gap, nm", "n_even - n_odd", "coupling length, um"])
    fit = np.polyfit(table["gap, nm"], np.log(table["n_even - n_odd"]), 1)
    residual = np.log(table["n_even - n_odd"]) - np.polyval(fit, table["gap, nm"])
    table.attrs["decay_nm"] = float(-1.0 / fit[0])
    table.attrs["largest_residual"] = float(np.max(np.abs(residual)))
    return table
