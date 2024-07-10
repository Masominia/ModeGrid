"""Static figures for the report (matplotlib, PNG for the README)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Validated categorical order (fixed, never cycled) + recessive ink.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
NEUTRAL = "#8a8985"
CRITICAL = "#c4271a"

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID,
        "axes.labelcolor": INK2,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "legend.frameon": False,
        "lines.linewidth": 1.6,
        "figure.dpi": 150,
        "axes.axisbelow": True,
    }
)


def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _unfold(result):
    """Field and permittivity on the whole cross-section, with the extent for imshow."""
    m = result["metrics"]
    field, eps = result["field"], result["eps"]
    x0, x1 = m["x0"], m["x0"] + m["nx"] * m["dx"]
    if m["x0"] == 0.0:  # right half only: mirror it
        sign = -1.0 if result.get("odd") else 1.0
        field = np.hstack([sign * field[:, ::-1], field])
        eps = np.hstack([eps[:, ::-1], eps])
        x0 = -x1
    return field, eps, (x0, x1, m["y0"], m["y0"] + m["ny"] * m["dy"])


def _profile(ax, result, zoom=True):
    field, eps, extent = _unfold(result)
    ax.imshow(np.abs(field), origin="lower", extent=extent, cmap="Blues", aspect="equal", interpolation="bilinear")
    xs = np.linspace(extent[0], extent[1], eps.shape[1] + 1)[:-1] + 0.5 * (extent[1] - extent[0]) / eps.shape[1]
    ys = np.linspace(extent[2], extent[3], eps.shape[0] + 1)[:-1] + 0.5 * (extent[3] - extent[2]) / eps.shape[0]
    if np.ptp(eps) > 0:
        ax.contour(xs, ys, eps, levels=[0.5 * (eps.min() + eps.max())], colors=[INK2], linewidths=0.8)
    if zoom:
        m = result["metrics"]
        half = min(extent[1], max(0.6, 4.0 * m["rms_x"]))
        tall = max(0.5, 4.0 * m["rms_y"])
        ax.set_xlim(-half, half)
        ax.set_ylim(max(extent[2], m["centre_y"] - tall), min(extent[3], m["centre_y"] + tall))
    ax.grid(False)
    ax.set(xlabel="x, um", ylabel="y, um")


def _cut(ax, result):
    field, eps, extent = _unfold(result)
    row = int(np.argmax(np.abs(field).max(axis=1)))
    xs = np.linspace(extent[0], extent[1], field.shape[1] + 1)[:-1] + 0.5 * (extent[1] - extent[0]) / field.shape[1]
    peak = np.abs(field[row]).max()
    ax.fill_between(xs, 0, np.where(eps[row] > 0.5 * (eps.min() + eps.max()), 1.1, 0.0), color=GRID, lw=0, label="core")
    ax.plot(xs, field[row] / peak, color=SERIES[0], label="field")
    half = min(extent[1], max(0.6, 4.0 * result["metrics"]["rms_x"]))
    ax.set_xlim(-half, half)
    ax.set_ylim(min(0.0, (field[row] / peak).min() * 1.1), 1.1)
    ax.set(xlabel="x, um", ylabel="field through its peak, normalised")
    ax.legend(loc="upper right", fontsize=7.5)


def case_figure(result, title, path, odd=None):
    result = dict(result, odd=("antisymmetric" in title) if odd is None else odd)
    fig, (left, right) = plt.subplots(1, 2, figsize=(8.4, 3.3))
    _profile(left, result)
    _cut(right, result)
    m = result["metrics"]
    fig.suptitle(f"{title}: n_eff {m['neff']:.4f}, group index {m['group_index']:.3f}", fontsize=9, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=105)
    plt.close(fig)


def overview(gallery, path):
    names = list(gallery)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.2))
    for ax, name in zip(axes.ravel(), names, strict=False):
        result, title = gallery[name]
        _profile(ax, dict(result, odd="antisymmetric" in title))
        short = title if len(title) < 50 else title[:47] + "..."
        ax.set_title(f"{short}\nn_eff {result['metrics']['neff']:.4f}", fontsize=8.5)
    for ax in axes.ravel()[len(names) :]:
        ax.axis("off")
    _save(fig, path)
