"""Figures for seed-connectivity results (fsLR surface grids, annotated heatmaps)."""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap

from .seedfc import cortex_to_surface

# Diverging blue <-> red with a neutral gray midpoint (dataviz reference palette)
DIVERGING = LinearSegmentedColormap.from_list("ocs_div", ["#0d366b", "#2a78d6", "#f0efec", "#e34948", "#8f1d1d"])
VIEWS = [("L", "lateral"), ("L", "medial"), ("R", "medial"), ("R", "lateral")]


def fslr_surfaces(kind: str = "veryinflated") -> dict[str, str]:
    import templateflow.api as tflow
    return {h: str(tflow.get("fsLR", density="32k", hemi=h, suffix=kind, extension=".surf.gii")) for h in "LR"}


def surface_grid(maps: np.ndarray, row_labels: list[str], bm, vmax: float, title: str, cbar_label: str,
                 threshold: float | None = None, path=None, surfaces: dict | None = None):
    """Rows = maps [R x G] on the fsLR surface (L lat, L med, R med, R lat), shared symmetric color scale.

    threshold: values with |v| < threshold are not drawn (surface shown in gray).
    """
    from nilearn import plotting as nplot
    surfaces = surfaces or fslr_surfaces()
    maps = np.atleast_2d(maps)
    fig, axs = plt.subplots(len(maps), 4, figsize=(12, 2.3 * len(maps)), subplot_kw={"projection": "3d"},
                            gridspec_kw=dict(wspace=0, hspace=0), squeeze=False)
    fig.subplots_adjust(top=1 - 0.35 / len(maps), bottom=0.5 / (2.3 * len(maps) + 1), left=0.14, right=0.99)
    for i, m in enumerate(maps):
        for j, (h, view) in enumerate(VIEWS):
            nplot.plot_surf_stat_map(surfaces[h], cortex_to_surface(m, bm, h), hemi={"L": "left", "R": "right"}[h],
                                     view=view, cmap=DIVERGING, vmin=-vmax, vmax=vmax, symmetric_cbar=True,
                                     threshold=threshold, colorbar=False, bg_on_data=False, axes=axs[i, j], figure=fig)
        axs[i, 0].text2D(-0.05, 0.5, row_labels[i], transform=axs[i, 0].transAxes, ha="right", va="center", fontsize=9)
    for j, (h, view) in enumerate(VIEWS):
        axs[0, j].set_title(f"{h} {view}", fontsize=9)
    cax = fig.add_axes([0.38, 0.25 / (2.3 * len(maps) + 1), 0.3, 0.01])
    fig.colorbar(ScalarMappable(cmap=DIVERGING, norm=plt.Normalize(-vmax, vmax)), cax=cax, orientation="horizontal",
                 label=cbar_label)
    fig.suptitle(title, y=0.995, fontsize=11)
    if path is not None:
        fig.savefig(path, dpi=200, bbox_inches="tight")
    return fig


def annotated_heatmap(values, row_labels, col_labels, title: str, cbar_label: str, marks=None,
                      lim: float | None = None, fmt: str = "{:.2f}", path=None):
    """Diverging heatmap with the value in each cell; ``marks`` (bool array) adds an asterisk."""
    values = np.asarray(values, dtype=float)
    lim = lim or float(np.nanmax(np.abs(values))) or 1.0
    fig, ax = plt.subplots(figsize=(1 + 0.75 * values.shape[1], 0.9 + 0.45 * values.shape[0]), layout="constrained")
    im = ax.imshow(values, cmap=DIVERGING, vmin=-lim, vmax=lim, aspect="auto")
    for (r, c), v in np.ndenumerate(values):
        if np.isfinite(v):
            star = "*" if marks is not None and marks[r, c] else ""
            ax.text(c, r, fmt.format(v) + star, ha="center", va="center", fontsize=8,
                    color="white" if abs(v) > 0.6 * lim else "#0b0b0b", fontweight="bold" if star else None)
    ax.set_xticks(range(values.shape[1]), col_labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=9)
    ax.set_yticks(range(values.shape[0]), row_labels, fontsize=9)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.colorbar(im, ax=ax, shrink=0.8, label=cbar_label)
    ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
    if path is not None:
        fig.savefig(path, dpi=200, bbox_inches="tight")
    return fig
