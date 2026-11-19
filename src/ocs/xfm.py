"""Apply templateflow's MNI152NLin2009cAsym <-> MNI152NLin6Asym ANTs transforms without ANTs.

The templateflow files (``tpl-MNI152NLin6Asym_from-MNI152NLin2009cAsym_mode-image_xfm.h5`` and the
reverse-named one) are ITK composite transforms: an affine plus a displacement field stored on the
2009cAsym 1 mm grid. Checked empirically against the template T1s, both map *2009c points to NLin6
points* (q = T(p)). Pulling a 2009c image onto an NLin6 grid therefore needs the inverse,
computed here by fixed-point iteration. Applied naively in the direction the file name suggests
(e.g. via nitransforms), the alignment gets worse (T1 r 0.69 -> 0.47) instead of better (-> 0.82).
"""
from __future__ import annotations

import numpy as np
import nibabel as nib
from scipy.ndimage import map_coordinates

_LPS = np.array([-1.0, -1.0, 1.0])   # ITK works in LPS; NIfTI world is RAS


class ITKComposite:
    """ITK composite (AffineTransform + DisplacementFieldTransform) from an ANTs .h5 file.

    ITK applies the queue [affine, field] last-in-first: field first, then affine.
    """

    def __init__(self, h5_path):
        import h5py
        with h5py.File(str(h5_path), "r") as f:
            types = [f[f"TransformGroup/{i}/TransformType"][0].decode() for i in (1, 2)]
            if not (types[0].startswith("AffineTransform") and types[1].startswith("DisplacementFieldTransform")):
                raise ValueError(f"unexpected transform types {types}")
            p = f["TransformGroup/1/TranformParameters"][()]
            self.center = f["TransformGroup/1/TranformFixedParameters"][()]
            fp = f["TransformGroup/2/TranformFixedParameters"][()]
            field = f["TransformGroup/2/TranformParameters"][()]
        self.M, self.t = p[:9].reshape(3, 3), p[9:12]
        size = fp[:3].astype(int)
        self.origin, spacing, direction = fp[3:6], fp[6:9], fp[9:18].reshape(3, 3)
        self.index_from_lps = np.linalg.inv(direction * spacing)
        self.D = field.reshape(size[2], size[1], size[0], 3).transpose(2, 1, 0, 3)   # x, y, z, vector

    def _field(self, x):
        idx = self.index_from_lps @ (x - self.origin).T
        return np.stack([map_coordinates(self.D[..., k], idx, order=1, mode="nearest") for k in range(3)], 1)

    def forward(self, ras: np.ndarray) -> np.ndarray:
        """Map RAS points [N x 3] with the stored transform (2009c -> NLin6 for the templateflow files)."""
        x = np.asarray(ras, float) * _LPS
        x = x + self._field(x)
        return ((x - self.center) @ self.M.T + self.t + self.center) * _LPS

    def inverse(self, ras: np.ndarray, iterations: int = 30, tol: float = 1e-4) -> np.ndarray:
        """Inverse mapping by fixed-point iteration (NLin6 -> 2009c for the templateflow files)."""
        q = np.asarray(ras, float)
        p = q.copy()
        for _ in range(iterations):
            err = self.forward(p) - q
            p -= err
            if np.abs(err).max() < tol:
                break
        return p


def resample(img: nib.Nifti1Image, ref: nib.Nifti1Image, ref_to_img, order: int = 1,
             chunk: int = 2_000_000) -> np.ndarray:
    """Pull ``img`` onto the grid of ``ref``; ``ref_to_img`` maps ref RAS points to img RAS points."""
    data = np.asarray(img.dataobj, dtype=float)
    inv = np.linalg.inv(img.affine)
    ijk = np.indices(ref.shape[:3]).reshape(3, -1).T
    out = np.empty(len(ijk))
    for s in range(0, len(ijk), chunk):
        src = nib.affines.apply_affine(inv, ref_to_img(nib.affines.apply_affine(ref.affine, ijk[s:s + chunk])))
        out[s:s + chunk] = map_coordinates(data, src.T, order=order, mode="constant")
    return out.reshape(ref.shape[:3])
