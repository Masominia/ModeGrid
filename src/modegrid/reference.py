"""An independent solution of every structure with SciPy's sparse eigen-solver, and analytic slab indices.

It shares no code with the C++ solver: its own structure reader, its own painting of the permittivity grid, a
vectorised assembly of the operator as a sparse matrix, and ARPACK in shift-invert mode instead of banded LU
with inverse iteration.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.sparse.linalg import eigs

MATERIALS = {
    "si": (3.476, -0.0865),
    "silicon": (3.476, -0.0865),
    "sio2": (1.444, -0.0120),
    "silica": (1.444, -0.0120),
    "oxide": (1.444, -0.0120),
    "sin": (1.996, -0.0240),
    "nitride": (1.996, -0.0240),
    "air": (1.0, 0.0),
}
GROUP_STEP = 0.005


def index(material: str, lam: float) -> float:
    if material in MATERIALS:
        n0, slope = MATERIALS[material]
        return n0 + slope * (lam - 1.55)
    return float(material)


def read_structure(path: Path | str) -> dict:
    """Reads a structure file into a dictionary, ignoring its expectations."""
    s = {
        "title": "",
        "wavelength": 1.55,
        "dx": 0.02,
        "dy": 0.02,
        "width": 3.0,
        "ymin": -1.0,
        "ymax": 1.22,
        "background": "sio2",
        "rects": [],
        "symmetry": "none",
        "field": "te",
        "coupler": False,
    }
    for raw in Path(path).read_text().splitlines():
        w = raw.split("#")[0].split()
        if not w or w[0] == "expect":
            continue
        if w[0] == "title":
            s["title"] = " ".join(w[1:])
        elif w[0] == "wavelength":
            s["wavelength"] = float(w[1])
        elif w[0] in ("grid", "window"):
            s.update({k: float(v) for k, v in (x.split("=") for x in w[1:])})
        elif w[0] in ("background", "symmetry", "field"):
            s[w[0]] = w[1]
        elif w[0] == "layer":
            s["rects"].append((w[1], -1e9, 1e9, float(w[2]), float(w[3])))
        elif w[0] == "rect":
            s["rects"].append((w[1], *(float(v) for v in w[2:6])))
        elif w[0] == "coupler":
            s["coupler"] = True
    return s


def permittivity(s: dict, lam: float, half: bool) -> tuple[np.ndarray, float, float]:
    """The grid of cell permittivities (rows bottom to top) and the coordinates of its lower left corner."""
    span = 0.5 * s["width"] if half else s["width"]
    nx, ny = round(span / s["dx"]), round((s["ymax"] - s["ymin"]) / s["dy"])
    x0 = 0.0 if half else -0.5 * s["width"]
    xe, ye = x0 + s["dx"] * np.arange(nx + 1), s["ymin"] + s["dy"] * np.arange(ny + 1)
    eps = np.full((ny, nx), index(s["background"], lam) ** 2)
    for material, xa, xb, ya, yb in s["rects"]:
        fx = np.clip(np.minimum(xe[1:], xb) - np.maximum(xe[:-1], xa), 0, None) / s["dx"]
        fy = np.clip(np.minimum(ye[1:], yb) - np.maximum(ye[:-1], ya), 0, None) / s["dy"]
        f = np.outer(fy, fx)
        f[f > 1 - 1e-9] = 1.0
        f[f < 1e-9] = 0.0
        eps = eps * (1 - f) + f * index(material, lam) ** 2
    return eps, x0, s["ymin"]


def operator(eps: np.ndarray, dx: float, dy: float, lam: float, field: str, symmetry: str) -> sp.csc_matrix:
    """The semi-vectorial operator: its eigenvalues are (k n_eff)^2."""
    ny, nx = eps.shape
    k2 = (2 * np.pi / lam) ** 2
    number = np.arange(nx * ny).reshape(ny, nx)
    diagonal = k2 * eps.copy()
    rows, cols, vals = [], [], []
    for axis, h, jump in ((1, dx, field == "te"), (0, dy, field == "tm")):
        for direction in (-1, 1):
            q = np.roll(eps, -direction, axis=axis)  # the neighbour's permittivity
            edge = [slice(None), slice(None)]
            edge[axis] = 0 if direction == -1 else -1
            edge = tuple(edge)
            q[edge] = eps[edge]  # beyond the window: an image of the cell itself
            off = (2 * q / (eps + q) if jump else np.ones_like(eps)) / h**2
            diagonal -= (2 * eps / (eps + q) if jump else np.ones_like(eps)) / h**2
            inside = np.ones_like(eps, dtype=bool)
            inside[edge] = False
            rows.append(number[inside])
            cols.append(np.roll(number, -direction, axis=axis)[inside])
            vals.append(off[inside])
            mirror = axis == 1 and direction == -1 and symmetry == "even"
            diagonal[edge] += off[edge] if mirror else -off[edge]
    rows.append(number.ravel())
    cols.append(number.ravel())
    vals.append(diagonal.ravel())
    return sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nx * ny, nx * ny))


def solve(eps: np.ndarray, dx: float, dy: float, lam: float, field: str, symmetry: str) -> tuple[float, np.ndarray]:
    """Effective index and field (unit norm, peak positive) of the mode of highest index."""
    k = 2 * np.pi / lam
    a = operator(eps, dx, dy, lam, field, symmetry)
    start = np.ones(a.shape[0])
    value, vector = eigs(a, k=1, sigma=1.0005 * k * k * eps.max(), v0=start, tol=1e-13, maxiter=5000)
    e = vector[:, 0].real
    e /= np.linalg.norm(e)
    if e[np.argmax(np.abs(e))] < 0:
        e = -e
    return float(np.sqrt(value[0].real) / k), e.reshape(eps.shape)


def column_index(column: np.ndarray, dy: float, lam: float, field: str) -> float:
    """Index of the slab mode of one grid column, or the material index of a uniform column."""
    if np.ptp(column) < 1e-12:
        return float(np.sqrt(column.max()))
    k2 = (2 * np.pi / lam) ** 2
    n = len(column)
    a = np.diag(k2 * column)
    for i in range(n):
        for j in (i - 1, i + 1):
            p, q = column[i], column[j] if 0 <= j < n else column[i]
            off = (2 * q / (p + q) if field == "tm" else 1.0) / dy**2
            a[i, i] -= (2 * p / (p + q) if field == "tm" else 1.0) / dy**2
            if 0 <= j < n:
                a[i, j] = off
            else:
                a[i, i] -= off
    return float(np.sqrt(np.linalg.eigvals(a).real.max() / k2))


def cladding_index(eps: np.ndarray, dy: float, lam: float, field: str, half: bool) -> float:
    """The index a mode must exceed to be confined by the structure, not by the walls of the window."""
    n = max(np.sqrt(max(eps[0, :].max(), eps[-1, :].max())), column_index(eps[:, -1], dy, lam, field))
    return float(n if half else max(n, column_index(eps[:, 0], dy, lam, field)))


def evaluate(s: dict) -> dict:
    """Every quantity the solver reports for a structure, computed here from scratch."""
    lam = s["wavelength"]
    symmetry = "even" if s["coupler"] else s["symmetry"]
    half = symmetry != "none"
    eps, x0, y0 = permittivity(s, lam, half)
    n, e = solve(eps, s["dx"], s["dy"], lam, s["field"], symmetry)
    ny, nx = eps.shape
    x, y = x0 + (np.arange(nx) + 0.5) * s["dx"], y0 + (np.arange(ny) + 0.5) * s["dy"]
    p = e**2
    total = p.sum()
    cx = x0 if half else float((p.sum(axis=0) * x).sum() / total)
    cy = float((p.sum(axis=1) * y).sum() / total)
    out = {
        "neff": n,
        "nx": nx,
        "ny": ny,
        "n_cladding": cladding_index(eps, s["dy"], lam, s["field"], half),
        "confinement": float(p[eps >= eps.max() * (1 - 1e-9)].sum() / total),
        "area": float((2 if half else 1) * total**2 / (p**2).sum() * s["dx"] * s["dy"]),
        "rms_x": float(np.sqrt((p.sum(axis=0) * (x - cx) ** 2).sum() / total)),
        "rms_y": float(np.sqrt((p.sum(axis=1) * (y - cy) ** 2).sum() / total)),
        "centre_y": cy,
        "field": e,
    }
    out["guided"] = int(n > out["n_cladding"] + 1e-9)
    side = []
    for v in (lam - GROUP_STEP, lam + GROUP_STEP):
        side.append(solve(permittivity(s, v, half)[0], s["dx"], s["dy"], v, s["field"], symmetry)[0])
    out["group_index"] = n - lam * (side[1] - side[0]) / (2 * GROUP_STEP)
    if s["coupler"]:
        n_odd = solve(eps, s["dx"], s["dy"], lam, s["field"], "odd")[0]
        out.update(n_even=n, n_odd=n_odd, splitting=n - n_odd, coupling_length=lam / (2 * abs(n - n_odd)))
    return out


def slab_index(lam: float, thickness: float, n_core: float, n_clad: float, n_sub: float, pol: str, order: int = 0):
    """Exact effective index of a three-layer slab mode, or None below cut-off."""
    k0 = 2 * np.pi / lam
    c_clad, c_sub = ((n_core / n_clad) ** 2, (n_core / n_sub) ** 2) if pol == "tm" else (1.0, 1.0)

    def mismatch(n):
        kappa = k0 * np.sqrt(n_core**2 - n**2)
        gamma, delta = k0 * np.sqrt(max(n**2 - n_clad**2, 0.0)), k0 * np.sqrt(max(n**2 - n_sub**2, 0.0))
        return kappa * thickness - np.arctan2(c_clad * gamma, kappa) - np.arctan2(c_sub * delta, kappa) - order * np.pi

    lo, hi = max(n_clad, n_sub) + 1e-12, n_core - 1e-12
    if mismatch(lo) <= 0:
        return None
    return brentq(mismatch, lo, hi, xtol=1e-14)


def effective_index_method(lam: float, width: float, thickness: float, n_core: float, n_clad: float, pol: str):
    """The classic estimate for a buried strip: a slab solve across the height, then one across the width."""
    first = slab_index(lam, thickness, n_core, n_clad, n_clad, pol)
    if first is None:
        return None
    return slab_index(lam, width, first, n_clad, n_clad, "tm" if pol == "te" else "te")
