"""Python side of the C++ solver: runs the binary on a structure file and parses its output."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / "cpp" / "build" / "modegrid"


def available() -> bool:
    return BINARY.exists()


def run(structure: Path | str, command: str = "run") -> dict:
    """Solve a structure. Returns the reported quantities, the checks and, for `run`, the grids.

    Raises RuntimeError with the solver's message if the file is invalid.
    """
    proc = subprocess.run([str(BINARY), command, str(structure)], capture_output=True, text=True)
    if proc.returncode == 2:
        raise RuntimeError(proc.stderr.strip())
    res: dict = {"metrics": {}, "checks": [], "title": ""}
    eps, field = [], []
    for line in proc.stdout.splitlines():
        key, _, rest = line.partition(" ")
        if key == "eps":
            eps.append(rest)
        elif key == "field":
            field.append(rest)
        elif key == "check":
            status, _, tail = rest.partition(" ")
            name, _, detail = tail.partition(" :: ")
            res["checks"].append({"check": name, "status": status, "detail": detail})
        elif key == "title":
            res["title"] = rest
        else:
            res["metrics"][key] = float(rest)
    if field:
        res["eps"] = np.array([r.split() for r in eps], dtype=float)  # rows bottom to top
        res["field"] = np.array([r.split() for r in field], dtype=float)
    return res
