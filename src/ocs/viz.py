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


# ---------------------------------------------------------------------------
# Sample counts (inventory and after motion QC)
# ---------------------------------------------------------------------------
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]   # timepoint slots, fixed order
C_MUTED, C_GRID, C_AX, C_MISSING = "#8a897f", "#e1e0d9", "#c3c2b7", "#e34948"
MOTION_FILES = ("dcan_hdf5", "motion_tsv")


def _style(ax):
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(axis="x", color=C_GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(C_AX)


def sample_counts(runs, motion_qc=None, session_qc=None):
    """Found vs passing counts per level (runs, sessions, subjects) and timepoint.

    Without motion_qc/session_qc, "passing" = has both motion files (DCAN HDF5 + motion.tsv); a session/subject
    passes if at least one of its runs does. With them, passing = run_ok / session_ok from motion QC.
    Returns a DataFrame: level, timepoint ('all' for subject rows), found, passing.
    """
    import pandas as pd
    if motion_qc is None:
        r = runs.assign(ok=runs[list(MOTION_FILES)].notna().all(axis=1))
        s = r.groupby(["subject", "session", "timepoint"], as_index=False).ok.any()
    else:
        r = motion_qc.rename(columns={"run_ok": "ok"})
        s = session_qc.rename(columns={"session_ok": "ok"})
    tps = sorted(s.timepoint.dropna().astype(int).unique())
    rows = []
    for tp in tps:
        d = r[r.timepoint == tp]; rows.append(("Runs", f"T{tp}", len(d), int(d.ok.sum())))
    for tp in tps:
        d = s[s.timepoint == tp]; rows.append(("Sessions", f"T{tp}", len(d), int(d.ok.sum())))
    by_sub = s.pivot_table(index="subject", columns="timepoint", values="ok", aggfunc="any").reindex(columns=tps)
    rows.append(("Subjects, any timepoint", "all", len(by_sub), int(by_sub.fillna(False).astype(bool).any(axis=1).sum())))
    if len(tps) > 1:
        rows.append((f"Subjects, all {len(tps)} timepoints", "all", int(by_sub.notna().all(axis=1).sum()),
                     int(by_sub.fillna(False).astype(bool).all(axis=1).sum())))
    return pd.DataFrame(rows, columns=["level", "timepoint", "found", "passing"])


def sample_counts_figure(runs, motion_qc=None, session_qc=None, min_run_s=None, max_mean_fd=None, path=None):
    """Two panels. A: found (light) vs passing (solid) runs / sessions / subjects per timepoint.
    B: inventory -> runs missing each file type (motion files highlighted);
       motion QC -> runs failing each inclusion criterion per timepoint.
    Returns (fig, counts DataFrame)."""
    from matplotlib.patches import Patch
    qc = motion_qc is not None
    counts = sample_counts(runs, motion_qc, session_qc)
    tps = sorted(runs.timepoint.dropna().astype(int).unique())
    col = {f"T{tp}": CATEGORICAL[i % len(CATEGORICAL)] for i, tp in enumerate(tps)} | {"all": "#55534c"}
    pass_lbl = "passing motion QC" if qc else "with motion files"

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 0.5 * max(len(counts), 7) + 1.6), layout="constrained",
                                   gridspec_kw=dict(width_ratios=[1.25, 1]))
    # A - found vs passing
    y = np.arange(len(counts))[::-1]
    xmax = counts.found.max()
    for yi, row in zip(y, counts.itertuples()):
        c = col[row.timepoint]
        axA.barh(yi, row.found, height=0.7, color=c, alpha=0.22, lw=0)
        axA.barh(yi, row.passing, height=0.7, color=c, lw=0)
        pct = f" ({100 * row.passing / row.found:.0f}%)" if row.found else ""
        axA.text(row.found + xmax * 0.01, yi, f"{row.passing} / {row.found}{pct}", va="center", fontsize=9, color="#2b2a26")
    axA.set_yticks(y, [f"{l}  {t}" if t != "all" else l for l, t in zip(counts.level, counts.timepoint)])
    axA.set_xlim(0, xmax * 1.25); axA.set_xlabel("count")
    axA.set_title(f"A. Found vs {pass_lbl}", loc="left", fontweight="bold")
    _style(axA)
    axA.legend(handles=[Patch(color=C_MUTED, alpha=0.22, label="found"), Patch(color=C_MUTED, label=pass_lbl)]
               + [Patch(color=col[f"T{tp}"], label=f"T{tp}") for tp in tps],
               frameon=False, loc="lower right", fontsize=8)

    # B - why runs are missing / failing
    if not qc:
        file_cols = [c for c in runs.columns if c not in ("subject", "session", "task", "run", "timepoint")]
        file_cols = [c for c in MOTION_FILES if c in file_cols] + [c for c in file_cols if c not in MOTION_FILES]
        miss = runs[file_cols].isna().sum()
        yb = np.arange(len(file_cols))[::-1]
        for yi, c in zip(yb, file_cols):
            m = int(miss[c])
            axB.barh(yi, m, height=0.6, color=C_MISSING if c in MOTION_FILES else C_MUTED, lw=0)
            axB.text(m + 0.02 * max(1, miss.max()), yi, f"{m} / {len(runs)}", va="center", fontsize=9, color="#2b2a26")
        axB.set_yticks(yb, [f"{c}  (motion)" if c in MOTION_FILES else c for c in file_cols])
        axB.set_xlim(0, max(1, miss.max()) * 1.3); axB.set_xlabel("runs missing the file")
        axB.set_title("B. Runs missing each file", loc="left", fontweight="bold")
    else:
        mq = motion_qc
        reasons = {"no dtseries": ~mq.has_dtseries,
                   "no motion files": mq.seconds_remaining.isna(),
                   f"< {min_run_s:g} s retained": mq.seconds_remaining.notna() & ~mq.ok_seconds,
                   f"mean FD > {max_mean_fd:g} mm": mq.mean_fd_remaining.notna() & ~mq.ok_mean_fd,
                   "any (run excluded)": ~mq.run_ok}
        h = 0.8 / len(tps)
        yb = np.arange(len(reasons))[::-1]
        top = 1
        for i, tp in enumerate(tps):
            sel = (mq.timepoint == tp).to_numpy()
            for yi, (lab, m) in zip(yb, reasons.items()):
                n = int(m.to_numpy()[sel].sum()); top = max(top, n)
                yy = yi + 0.4 - h * (i + 0.5)
                axB.barh(yy, n, height=h * 0.9, color=col[f"T{tp}"], lw=0)
                axB.text(n, yy, f" {n}", va="center", fontsize=8, color="#2b2a26")
        axB.set_yticks(yb, list(reasons)); axB.set_xlim(0, top * 1.25)
        axB.set_xlabel("runs (criteria can overlap)")
        axB.set_title("B. Runs failing each criterion", loc="left", fontweight="bold")
    _style(axB)
    stage = "after motion QC" if qc else "inventory"
    fig.suptitle(f"Sample {stage}: {runs.subject.nunique()} subjects, "
                 f"{len(runs[['subject', 'session']].drop_duplicates())} sessions, {len(runs)} runs found", fontsize=13)
    if path is not None:
        fig.savefig(path, dpi=200, bbox_inches="tight")
    return fig, counts
