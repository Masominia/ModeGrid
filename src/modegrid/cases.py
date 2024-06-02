"""The waveguide library: discovery, description and per-structure results."""

from __future__ import annotations

import json
from pathlib import Path

from . import engine


def discover(directory: Path) -> list[Path]:
    return sorted(directory.glob("*.wg"))


def notes(path: Path) -> list[str]:
    """The comment lines at the top of a structure: what it shows and what it is checked against."""
    out = []
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        out.append(line[1:].strip())
    return out


def family_of(name: str, families: dict[str, list[str]]) -> str:
    for family, prefixes in families.items():
        if any(name.startswith(p) for p in prefixes):
            return family
    return "other"


def run_case(path: Path) -> tuple[dict, dict]:
    """Solve one structure. Returns a record without the grids, and the raw result for plotting."""
    res = engine.run(path)
    record = {
        "case": path.stem,
        "title": res["title"],
        "notes": notes(path),
        "metrics": res["metrics"],
        "checks": res["checks"],
    }
    return record, res


def write_record(record: dict, directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{record['case']}.json"
    path.write_text(json.dumps(record, indent=1) + "\n")
    return path
