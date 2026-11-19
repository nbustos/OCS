"""Seed-based functional connectivity on fsLR 91k CIFTI data."""
from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np


def seed_corr(target: np.ndarray, seeds: np.ndarray) -> np.ndarray:
    """Fisher-z correlation of each seed time series with every target grayordinate.

    target: [T x G] (e.g. smoothed denoised dtseries); seeds: [T x S] seed mean time series.
    Returns [S x G] float32. Grayordinates with zero variance or non-finite data get NaN.
    """
    X = np.asarray(target, dtype=np.float64)
    S = np.asarray(seeds, dtype=np.float64).reshape(len(X), -1)
    if len(S) != len(X):
        raise ValueError(f"seed length {len(S)} != target length {len(X)}")
    Xc = X - X.mean(axis=0)
    Sc = S - S.mean(axis=0)
    xs, ss = np.sqrt((Xc ** 2).sum(axis=0)), np.sqrt((Sc ** 2).sum(axis=0))
    with np.errstate(divide="ignore", invalid="ignore"):
        r = (Sc.T @ Xc) / np.outer(ss, xs)
    r[:, ~np.isfinite(xs) | (xs == 0)] = np.nan
    return np.arctanh(np.clip(r, -0.999999, 0.999999)).astype(np.float32)


def write_dscalar(maps: np.ndarray, names: list[str], bm, path: str | Path) -> Path:
    """Write [S x G] maps as a dscalar.nii (one named map per row) on brain model ``bm``."""
    maps = np.asarray(maps, dtype=np.float32)
    if maps.shape != (len(names), len(bm)):
        raise ValueError(f"maps {maps.shape} != ({len(names)}, {len(bm)})")
    path = Path(path)
    nib.Cifti2Image(maps, (nib.cifti2.ScalarAxis(names), bm)).to_filename(str(path))
    return path


def read_dscalar(path: str | Path) -> tuple[np.ndarray, list[str], object]:
    """Return (maps [S x G], names, brain model) from a dscalar.nii."""
    img = nib.load(str(path))
    return np.asarray(img.dataobj, dtype=np.float32), [str(n) for n in img.header.get_axis(0).name], img.header.get_axis(1)


def cortex_to_surface(values: np.ndarray, bm, hemi: str) -> np.ndarray:
    """Map a 91k grayordinate vector onto the full fsLR 32k vertex array of one hemisphere (NaN = medial wall)."""
    struct = {"L": "CIFTI_STRUCTURE_CORTEX_LEFT", "R": "CIFTI_STRUCTURE_CORTEX_RIGHT"}[hemi]
    for name, slc, sub in bm.iter_structures():
        if name == struct:
            out = np.full(sub.nvertices[name], np.nan, dtype=np.float32)
            out[sub.vertex] = np.asarray(values)[slc]
            return out
    raise ValueError(f"{struct} not in brain model")
