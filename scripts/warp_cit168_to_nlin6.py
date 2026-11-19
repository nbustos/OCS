#!/usr/bin/env python
"""Warp the CIT168 deterministic atlas from MNI152NLin2009cAsym into MNI152NLin6Asym (the CIFTI 91k
subcortical space) with the templateflow transform, and verify the direction on the template T1s.

Writes resources/atlases/CIT168_Pauli2018/CIT168_det_MNI152NLin6Asym_res-0.7mm.nii.gz
(same 0.7 mm box as the 2009c file, nearest-neighbour labels).

    python scripts/warp_cit168_to_nlin6.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import nibabel as nib
import templateflow.api as tflow

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ocs.paths import RESOURCES          # noqa: E402
from ocs.xfm import ITKComposite, resample   # noqa: E402

CIT = RESOURCES / "atlases" / "CIT168_Pauli2018"
SRC = CIT / "CIT168_det_MNI152NLin2009cAsym_res-0.7mm.nii.gz"
OUT = CIT / "CIT168_det_MNI152NLin6Asym_res-0.7mm.nii.gz"


def main():
    h5 = tflow.get("MNI152NLin6Asym", **{"from": "MNI152NLin2009cAsym"}, mode="image", suffix="xfm",
                   extension=".h5")
    X = ITKComposite(h5)

    # direction check: 2009c T1 pulled onto the NLin6 grid must match the NLin6 T1 better than no warp
    n6 = nib.load(str(tflow.get("MNI152NLin6Asym", resolution=1, suffix="T1w", desc=None, extension=".nii.gz")))
    c9 = nib.load(str(tflow.get("MNI152NLin2009cAsym", resolution=1, suffix="T1w", desc=None, extension=".nii.gz")))
    r = n6.get_fdata(); m = r > np.percentile(r, 60)
    corr = lambda a: float(np.corrcoef(a[m], r[m])[0, 1])
    r_id, r_inv = corr(resample(c9, n6, lambda p: p)), corr(resample(c9, n6, X.inverse))
    print(f"T1 match to NLin6: no warp r={r_id:.3f}  warped r={r_inv:.3f}")
    if r_inv <= r_id + 0.05:
        sys.exit("warp did not improve the template match - transform direction is wrong, nothing written")

    det = nib.load(str(SRC))
    ref = nib.Nifti1Image(np.zeros(det.shape, np.int16), det.affine)
    lab = np.rint(resample(det, ref, X.inverse, order=0)).astype(np.int16)
    out = nib.Nifti1Image(lab, det.affine)
    out.header.set_xyzt_units("mm"); out.set_qform(det.affine, 4); out.set_sform(det.affine, 4)   # 4 = MNI
    out.header["descrip"] = b"CIT168 det warped 2009cAsym->NLin6Asym (templateflow xfm)"
    out.to_filename(str(OUT))
    a, b = np.asarray(det.dataobj), lab
    print(f"wrote {OUT}")
    print("label  voxels 2009c -> NLin6")
    for k in range(1, 17):
        print(f"  {k:2d}   {int((a == k).sum()):6d} -> {int((b == k).sum()):6d}")


if __name__ == "__main__":
    main()
