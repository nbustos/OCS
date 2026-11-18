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
    """One row per FD threshold with frames/seconds remaining and mean FD (via rsx)."""
    return pd.DataFrame([read_motion_metrics(hdf5_path, t) for t in dcan_thresholds(hdf5_path)])


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
