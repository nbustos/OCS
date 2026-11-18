"""Map volumetric ROIs onto CIFTI 91k grayordinates and compute ROI-level QC."""
from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd


def load_dtseries(path: str | Path):
    """Return (data [T x G] float32, BrainModelAxis, TR seconds)."""
    img = nib.load(str(path))
    series, bm = img.header.get_axis(0), img.header.get_axis(1)
    return np.asarray(img.get_fdata(dtype=np.float32)), bm, float(series.step)


def volume_grayordinates(bm) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Indices, MNI coordinates (mm) and structure names of the volumetric grayordinates."""
    idx = np.where(bm.volume_mask)[0]
    ijk = bm.voxel[idx]
    xyz = nib.affines.apply_affine(bm.affine, ijk)
    return idx, xyz, np.asarray(bm.name)[idx]


def mask_from_nifti(img_or_path, bm, threshold: float = 0.5, labels=None,
                    min_coverage: float | None = None, subsamples: int = 5) -> np.ndarray:
    """Grayordinate indices whose MNI position falls in a NIfTI ROI.

    Probabilistic/binary map: value >= threshold. Label map: value in ``labels``.
    The NIfTI may be on any grid as long as it is in the same MNI space as the CIFTI
    volume (NLin6Asym).

    min_coverage=None (default): nearest-voxel sampling at each grayordinate centre.
    min_coverage=f: each 2 mm grayordinate voxel is sampled on a ``subsamples``^3 grid and kept if
    at least fraction f of it lies in the ROI. Use for masks much finer than 2 mm (e.g. 1 mm
    NM-MRI masks), where centre sampling keeps only a handful of grayordinates.
    """
    img = nib.load(str(img_or_path)) if isinstance(img_or_path, (str, Path)) else img_or_path
    data = np.asarray(img.dataobj)
    if data.ndim == 4:
        raise ValueError("Pass a 3D image (index a single volume of 4D atlases first)")
    idx = np.where(bm.volume_mask)[0]
    ijk = bm.voxel[idx].astype(float)
    if min_coverage is None:
        offsets = np.zeros((1, 3))
    else:
        o = (np.arange(subsamples) + 0.5) / subsamples - 0.5
        offsets = np.stack(np.meshgrid(o, o, o, indexing="ij"), -1).reshape(-1, 3)
    inv = np.linalg.inv(img.affine) @ bm.affine                       # CIFTI voxel -> NIfTI voxel
    hits = np.zeros(len(idx))
    for off in offsets:
        vox = np.rint(nib.affines.apply_affine(inv, ijk + off)).astype(int)
        inside = np.all((vox >= 0) & (vox < np.array(data.shape)), axis=1)
        vals = np.zeros(len(idx))
        vals[inside] = data[tuple(vox[inside].T)]
        hits += np.isin(vals, list(labels)) if labels is not None else vals >= threshold
    frac = hits / len(offsets)
    keep = frac >= min_coverage if min_coverage is not None else frac > 0
    return idx[keep]


def split_hemi(idx: np.ndarray, bm) -> dict[str, np.ndarray]:
    """Split volumetric grayordinate indices into L (x<0) / R (x>0) by MNI x."""
    x = nib.affines.apply_affine(bm.affine, bm.voxel[idx])[:, 0]
    return {"L": idx[x < 0], "R": idx[x > 0]}


def roi_qc(data: np.ndarray, idx: np.ndarray, raw: np.ndarray | None = None) -> dict:
    """QC metrics for one ROI in one run.

    data: denoised dtseries [T x G]; raw: optional un-denoised (fMRIPrep) dtseries for tSNR.
    - frac_valid: ROI grayordinates with finite, non-constant signal (coverage / dropout)
    - sd_ratio: median temporal SD in ROI / median SD over all valid grayordinates
    - homogeneity: mean pairwise correlation between ROI grayordinates
    - pc1_var: variance explained by the first principal component of the ROI
    - r_global: correlation of ROI mean with the global (all-grayordinate) mean signal
    - tsnr_median: median temporal SNR in the ROI from the raw data (if given)
    """
    out = {"n_grayordinates": int(len(idx))}
    if len(idx) == 0:
        return out
    X = data[:, idx]
    sd = X.std(0)
    valid = np.isfinite(X).all(0) & (sd > 0)
    allsd = data.std(0)
    allvalid = np.isfinite(data).all(0) & (allsd > 0)
    out.update(n_valid=int(valid.sum()), frac_valid=float(valid.mean()))
    if valid.sum() == 0:
        return out
    Xv = X[:, valid]
    out["sd_ratio"] = float(np.median(sd[valid]) / np.median(allsd[allvalid]))
    if valid.sum() >= 2:
        Z = (Xv - Xv.mean(0)) / Xv.std(0)
        C = (Z.T @ Z) / len(Z)
        out["homogeneity"] = float(C[np.triu_indices_from(C, 1)].mean())
        ev = np.linalg.eigvalsh(C)
        out["pc1_var"] = float(ev[-1] / ev.sum())
    roi_mean = Xv.mean(1)
    gs = data[:, allvalid].mean(1)
    out["r_global"] = float(np.corrcoef(roi_mean, gs)[0, 1])
    if raw is not None:
        R = raw[:, idx][:, valid]
        rsd = R.std(0)
        out["tsnr_median"] = float(np.median(np.where(rsd > 0, R.mean(0) / rsd, np.nan)))
    return out


def roi_timeseries(data: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Mean time series over valid (finite, non-constant) ROI grayordinates."""
    X = data[:, idx]
    valid = np.isfinite(X).all(0) & (X.std(0) > 0)
    return X[:, valid].mean(1) if valid.any() else np.full(len(data), np.nan)
