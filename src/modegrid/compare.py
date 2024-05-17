"""Compares what the solver reports with the independent SciPy solution."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import reference

# quantity: absolute tolerance (relative for the area and the coupler quantities marked below)
TOLERANCES = {
    "neff": 1e-9,
    "group_index": 1e-6,
    "n_cladding": 1e-7,
    "guided": 0.0,
    "confinement": 1e-7,
    "area": 2e-5,
    "rms_x": 2e-5,
    "rms_y": 1e-6,
    "centre_y": 1e-7,
    "n_even": 1e-9,
    "n_odd": 1e-9,
    "splitting": 2e-9,
    "coupling_length": 1e-5,
}
RELATIVE = {"area", "coupling_length"}


def compare(path: Path | str, result: dict) -> dict:
    """Solves one structure again. Returns the quantities out of tolerance and the largest differences."""
    ref = reference.evaluate(reference.read_structure(path))
    got = result["metrics"]
    differences, outside = {}, []
    for name, tol in TOLERANCES.items():
        if name not in ref:
            continue
        diff = abs(got[name] - ref[name]) / (abs(ref[name]) if name in RELATIVE else 1.0)
        differences[name] = float(diff)
        if diff > tol:
            outside.append(f"{name}: {got[name]:.10g} against {ref[name]:.10g}")
    field = float(np.max(np.abs(result["field"] - ref["field"])) / np.max(np.abs(ref["field"])))
    if field > 1e-4:
        outside.append(f"field differs by {field:.2e} of its peak")
    return {"outside": outside, "differences": differences, "field": field}
