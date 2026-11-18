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
