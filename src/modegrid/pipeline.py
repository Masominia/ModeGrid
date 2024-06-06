"""Solves the whole waveguide library, checks it and writes the report."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

from . import cases, engine, plotting, studies
from .compare import compare
from .engine import ROOT
from .report import md_table


def run(cfg_path=None, out=None, fast: bool = False) -> dict:
    cfg = yaml.safe_load(open(cfg_path or ROOT / "configs" / "default.yaml"))
    out = Path(out) if out else ROOT / "reports"
    (out / "figures" / "cases").mkdir(parents=True, exist_ok=True)
    if not engine.available():
        raise SystemExit("C++ solver not built: run `make build` first")
    paths = cases.discover(ROOT / cfg["cases_dir"])
    if fast:
        keep = tuple(cfg["fast_families"])
        paths = [
            p
            for i, p in enumerate(paths)
            if i % cfg["fast_every"] == 0 or p.stem in cfg["gallery"] or p.stem.startswith(keep)
        ]

    checks, rows, gallery, worst, all_metrics = [], [], {}, {}, {}
    for path in paths:
        record, result = cases.run_case(path)
        family = cases.family_of(path.stem, cfg["families"])
        cmp = compare(path, result)
        record["against_reference"] = {"largest_differences": cmp["differences"], "field": cmp["field"]}
        for name, diff in cmp["differences"].items():
            worst[name] = max(worst.get(name, 0.0), diff)
        worst["field"] = max(worst.get("field", 0.0), cmp["field"])
        checks.append(
            {
                "area": "independent solution",
                "check": f"{path.stem}: every reported quantity and the field solved again with SciPy",
                "status": "FAIL" if cmp["outside"] else "PASS",
                "detail": "; ".join(cmp["outside"])
                or f"{len(cmp['differences'])} quantities agree, field within {cmp['field']:.1e} of its peak",
            }
        )
        m = result["metrics"]
        checks.append(
            {
                "area": "solver",
                "check": f"{path.stem}: the iteration converged",
                "status": "PASS" if m["converged"] == 1 else "FAIL",
                "detail": f"{int(m['iterations'])} iterations on {int(m['nx'] * m['ny'])} cells",
            }
        )
        cases.write_record(record, out / "cases")
        plotting.case_figure(result, record["title"], out / "figures" / "cases" / f"{path.stem}.png")
        for c in record["checks"]:
            checks.append(
                {"area": family, "check": f"{path.stem}: {c['check']}", "status": c["status"], "detail": c["detail"]}
            )
        all_metrics[path.stem] = m
        rows.append(
            {
                "structure": path.stem,
                "family": family,
                "cells": int(m["nx"] * m["ny"]),
                "n_eff": m["neff"],
                "group index": m["group_index"],
                "guided": "yes" if m["guided"] else "no",
                "confinement": m["confinement"],
                "area, um2": m["area"],
                "checks": len(record["checks"]),
                "passed": sum(c["status"] == "PASS" for c in record["checks"]),
            }
        )
        if path.stem in cfg["gallery"]:
            gallery[path.stem] = (result, record["title"])

    table = pd.DataFrame(rows)
    fam = table.groupby("family", sort=False).agg(
        structures=("structure", "size"), checks=("checks", "sum"), passed=("passed", "sum")
    )
    if gallery:
        plotting.overview({k: gallery[k] for k in cfg["gallery"] if k in gallery}, out / "figures" / "overview.png")

    # studies across structures
    conv = {pol: studies.convergence(all_metrics, pol) for pol in ("te", "tm")}
    sweep = studies.width_sweep(all_metrics)
    pairs = studies.couplers(all_metrics)
    study_lines, study_metrics = [], {}
    for pol, t in conv.items():
        if len(t):
            values = t["n_eff"].to_numpy()
            monotone = bool((values[1:] > values[:-1]).all())
            at20 = float(t.loc[(t["cell, nm"] - 20).abs().idxmin(), "above the extrapolated value"])
            checks.append(
                {
                    "area": "grid convergence",
                    "check": f"quasi-{pol.upper()} index falls steadily as the grid is refined",
                    "status": "PASS" if monotone else "FAIL",
                    "detail": f"extrapolated {t.attrs['limit']:.5f}; a 20 nm grid is {at20:+.4f} above it",
                }
            )
            study_lines += [
                f"### Grid convergence, quasi-{pol.upper()} (480 x 240 nm silicon strip)",
                "",
                f"First-order extrapolation to zero cell size: {t.attrs['limit']:.5f}.",
                "",
                md_table(t, 5, index=False),
                "",
            ]
            study_metrics[f"convergence_{pol}"] = {"limit": t.attrs["limit"], "rows": t.to_dict("records")}
    if len(sweep) >= 3:
        wide = sweep[sweep["width, nm"] >= 400]
        worst_eim = float(wide["difference"].abs().max()) if len(wide) else float("nan")
        checks.append(
            {
                "area": "physical consistency",
                "check": "effective index method agrees with the solver for strips of 400 nm and wider",
                "status": "PASS" if worst_eim < 0.06 else "FAIL",
                "detail": f"largest difference {worst_eim:.4f}",
            }
        )
        rising = bool((sweep["n_eff"].to_numpy()[1:] > sweep["n_eff"].to_numpy()[:-1]).all())
        checks.append(
            {
                "area": "physical consistency",
                "check": "index rises with strip width",
                "status": "PASS" if rising else "FAIL",
                "detail": f"{len(sweep)} widths",
            }
        )
        study_lines += [
            "### Silicon strips against width, quasi-TE (220 nm thick, 1.55 um)",
            "",
            md_table(sweep, 4, index=False),
            "",
        ]
        study_metrics["width_sweep"] = sweep.to_dict("records")
    if len(pairs):
        checks.append(
            {
                "area": "physical consistency",
                "check": "supermode splitting of two strips falls exponentially with the gap",
                "status": "PASS" if pairs.attrs["largest_residual"] < 0.08 else "FAIL",
                "detail": f"decay length {pairs.attrs['decay_nm']:.1f} nm, largest deviation from the fit {pairs.attrs['largest_residual']:.3f} in the logarithm",
            }
        )
        study_lines += [
            "### Directional coupler, quasi-TE (two 450 x 220 nm strips)",
            "",
            f"Splitting decays by 1/e every {pairs.attrs['decay_nm']:.1f} nm of gap.",
            "",
            md_table(pairs, 5, index=False),
            "",
        ]
        study_metrics["couplers"] = {"decay_nm": pairs.attrs["decay_nm"], "rows": pairs.to_dict("records")}
    if "mirror_full" in all_metrics and "mirror_half" in all_metrics:
        d = abs(all_metrics["mirror_full"]["neff"] - all_metrics["mirror_half"]["neff"])
        checks.append(
            {
                "area": "physical consistency",
                "check": "whole structure and mirrored half give the same index",
                "status": "PASS" if d < 1e-9 else "FAIL",
                "detail": f"difference {d:.1e}",
            }
        )

    val = pd.DataFrame(checks)
    counts = val["status"].value_counts().to_dict()
    agreement = pd.DataFrame({"largest difference over the library": worst}).rename_axis("quantity")
    lines = [
        "# ModeGrid report",
        "",
        f"{len(table)} structures, {int(table['cells'].sum())} cells solved, {len(val)} checks: {counts.get('PASS', 0)} PASS, {counts.get('FAIL', 0)} FAIL.",
        "",
        "Every structure is a text file in `cases/`. Its expectations come from an independent solution of the same grid or from a physical bound, stated in the comments at the top of the file.",
        "Per-structure results are in `reports/cases/`, mode profiles in `reports/figures/cases/`.",
        "",
        "## Families",
        "",
        md_table(fam, 0),
        "",
        "## Agreement with the independent solution",
        "",
        "Each structure is read, painted, assembled and solved again by `reference.py` with SciPy's sparse eigen-solver. Differences are absolute, except the area and the coupling length (relative) and the field (relative to its peak).",
        "",
        md_table(agreement, 12),
        "",
        "## Studies",
        "",
        *study_lines,
        "## All structures",
        "",
        md_table(table.drop(columns=["family"]), 4, index=False),
        "",
    ]
    (out / "report.md").write_text("\n".join(lines))
    (out / "validation.md").write_text("# Validation checks\n\n" + md_table(val, index=False) + "\n")
    metrics = {
        "structures": len(table),
        "cells": int(table["cells"].sum()),
        "families": fam.to_dict("index"),
        "largest_differences": worst,
        "studies": study_metrics,
        "structures_table": table.set_index("structure").to_dict("index"),
        "checks": checks,
    }
    (out / "metrics.json").write_text(json.dumps(metrics, indent=1, default=float))
    return metrics
