"""Machine-specific paths. Override with environment variables instead of editing code."""

import os
from pathlib import Path

XCPD_DIR = Path(
    os.environ.get(
        "OCS_XCPD_DIR",
        "/Volumes/Synology/1_current_projects/OCS_7006/5_Data/"
        "marsonc_scan_data/marsonc/derivatives/xcp_d",
    )
)
FMRIPREP_DIR = Path(os.environ.get("OCS_FMRIPREP_DIR", XCPD_DIR.parent / "fmriprep"))

REPO_DIR = Path(__file__).resolve().parents[2]
RESOURCES = REPO_DIR / "resources"

# NM-MRI template masks (shared on request, never in this repo). Order: $OCS_NM_MASKS, the Synology
# folder on cortex if mounted, else the git-ignored local copy in resources/private/.
NM_MASKS_CORTEX = Path(
    "/Volumes/Synology/1_current_projects/AKILI_OCD_AAAU8631/5_Data/"
    "akili_scan_data/akili/scripts/NM/Masks"
)
NM_MASKS_DIR = Path(
    os.environ.get("OCS_NM_MASKS")
    or (NM_MASKS_CORTEX if NM_MASKS_CORTEX.is_dir() else RESOURCES / "private" / "akili_NM_masks")
)


def nm_mask(name: str) -> Path:
    """Path of an NM mask by stem (e.g. 'SN_mask'), accepting .nii.gz or .nii."""
    for ext in (".nii.gz", ".nii"):
        p = NM_MASKS_DIR / f"{name}{ext}"
        if p.exists():
            return p
    return NM_MASKS_DIR / f"{name}.nii.gz"   # missing: reported by CELL 1
