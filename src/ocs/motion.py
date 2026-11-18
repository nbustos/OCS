"""Motion QC from XCP-D DCAN HDF5 files and motion.tsv."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from rsx_analyzer import read_motion_metrics


def dcan_thresholds(hdf5_path: str | Path) -> list[float]:
    """FD thresholds present in a DCAN file (groups ``dcan_motion/fd_<thr>``)."""
    import h5py
    with h5py.File(hdf5_path, "r") as h:
        if "dcan_motion" not in h:
            raise KeyError(f"No 'dcan_motion' group in {hdf5_path}")
        return sorted(float(k[3:]) for k in h["dcan_motion"].keys() if k.startswith("fd_"))


def read_dcan_all(hdf5_path: str | Path) -> pd.DataFrame:
    """One row per FD threshold with frames/seconds remaining and mean FD.

    Same columns and checks as ``rsx_analyzer.read_motion_metrics``, but reads every
    ``dcan_motion/fd_<thr>`` group by its actual name in one pass. (XCP-D names groups by
    ``str(float)``, e.g. ``fd_0.0`` / ``fd_1.0``, which the ``f"fd_{thr:g}"`` lookup misses.)
    """
    import h5py
    rows = []
    with h5py.File(hdf5_path, "r") as h:
        if "dcan_motion" not in h:
            raise KeyError(f"No 'dcan_motion' group in {hdf5_path}")
        for key, g in h["dcan_motion"].items():
            if not key.startswith("fd_"):
                continue
            v = {f: float(np.asarray(g[f][()]).reshape(-1)[0]) for f in
                 ("remaining_frame_mean_FD", "remaining_seconds", "remaining_total_frame_count",
                  "skip", "threshold", "total_frame_count")}
            if not np.isclose(v["threshold"], float(key[3:])):
                raise ValueError(f"Stored threshold {v['threshold']} != group {key} in {hdf5_path}")
            total, remaining = v["total_frame_count"], v["remaining_total_frame_count"]
            if total < 0 or remaining < 0 or remaining > total:
                raise ValueError(f"Invalid frame counts in {hdf5_path}:{key}")
            rows.append({"fd_threshold": v["threshold"], "mean_fd_remaining": v["remaining_frame_mean_FD"],
                         "seconds_remaining": v["remaining_seconds"], "frames_remaining": remaining,
                         "frames_total": total, "frames_censored": total - remaining,
                         "fraction_remaining": remaining / total if total else np.nan, "skip": v["skip"]})
    return pd.DataFrame(rows).sort_values("fd_threshold").reset_index(drop=True)


def fd_summary(motion_tsv: str | Path, thresholds=(0.2, 0.3, 0.5)) -> dict:
    """Uncensored FD summary from XCP-D motion.tsv (first frame FD is 0/NaN by design)."""
    fd_full = pd.read_csv(motion_tsv, sep="\t")["framewise_displacement"].to_numpy(float)
    fd = fd_full[1:]
    out = {"n_frames_motion_tsv": len(fd_full), "mean_fd_all": np.nanmean(fd),
           "median_fd_all": np.nanmedian(fd), "max_fd_all": np.nanmax(fd)}
    for t in thresholds:
        out[f"pct_fd_gt_{t:g}"] = 100 * np.nanmean(fd > t)
    return out


def matching_thresholds(dcan: pd.DataFrame, n_frames_dtseries: int) -> list[float]:
    """FD thresholds whose DCAN remaining-frame count equals the dtseries length.

    XCP-D <0.7 removes censored volumes, so the threshold XCP-D applied is among these.
    Low-motion runs match many thresholds (nothing censored at any of them); intersect
    across runs with ``consensus_threshold`` to pin down the one actually used.
    """
    return sorted(float(t) for t in dcan.loc[dcan["frames_remaining"] == n_frames_dtseries, "fd_threshold"])


def consensus_threshold(per_run: list[list[float]]) -> float | None:
    """Threshold present in every run's match list (the highest if several remain)."""
    sets = [set(m) for m in per_run if m]
    common = set.intersection(*sets) if sets else set()
    return max(common) if common else None
