"""Find XCP-D (and matching fMRIPrep) run files and return one row per run."""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

RUN_RE = re.compile(r"sub-(?P<subject>[^_]+)_ses-(?P<session>[^_]+)_task-(?P<task>[^_]+)_run-(?P<run>[^_]+)")

# suffix (after the run prefix) -> column name
XCPD_FILES = {
    "space-fsLR_den-91k_desc-denoised_bold.dtseries.nii": "dtseries",
    "space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii": "dtseries_smoothed",
    "desc-dcan_qc.hdf5": "dcan_hdf5",
    "motion.tsv": "motion_tsv",
}
FMRIPREP_FILE = "space-fsLR_den-91k_bold.dtseries.nii"


def find_runs(xcpd_dir: str | Path, task: str = "rest", atlases: tuple[str, ...] = (),
              fmriprep_dir: str | Path | None = None) -> pd.DataFrame:
    """One row per sub/ses/task/run with paths (NaN when a file is missing).

    Runs are discovered from any of the XCP-D files listed in ``XCPD_FILES``, so a run
    with e.g. a DCAN file but no dtseries still appears (with the dtseries missing).
    """
    xcpd_dir = Path(xcpd_dir)
    rows: dict[tuple, dict] = {}
    for suffix, col in XCPD_FILES.items():
        for p in sorted(xcpd_dir.glob(f"sub-*/ses-*/func/sub-*_task-{task}_run-*_{suffix}")):
            m = RUN_RE.match(p.name)
            if not m:
                continue
            key = (m["subject"], m["session"], m["task"], m["run"])
            rows.setdefault(key, dict(zip(("subject", "session", "task", "run"), key)))[col] = str(p)
    if not rows:
        raise FileNotFoundError(f"No task-{task} XCP-D run files under {xcpd_dir}")
    df = pd.DataFrame(rows.values())
    for col in XCPD_FILES.values():
        if col not in df:
            df[col] = pd.NA
    for atlas in atlases:
        col = f"ts_{atlas}"
        df[col] = [_existing(xcpd_dir / f"sub-{r.subject}/ses-{r.session}/func/"
                             f"sub-{r.subject}_ses-{r.session}_task-{r.task}_run-{r.run}"
                             f"_space-fsLR_atlas-{atlas}_timeseries.tsv") for r in df.itertuples()]
    if fmriprep_dir is not None:
        fp = Path(fmriprep_dir)
        df["fmriprep_dtseries"] = [_existing(fp / f"sub-{r.subject}/ses-{r.session}/func/"
                                             f"sub-{r.subject}_ses-{r.session}_task-{r.task}_run-{r.run}_{FMRIPREP_FILE}")
                                   for r in df.itertuples()]
    df = df.sort_values(["subject", "session", "run"]).reset_index(drop=True)
    # timepoint = order of sessions within subject (1 = first session)
    order = (df[["subject", "session"]].drop_duplicates()
             .assign(timepoint=lambda d: d.groupby("subject")["session"].rank(method="dense").astype(int)))
    return df.merge(order, on=["subject", "session"], how="left")


def _existing(path: Path):
    return str(path) if path.exists() else pd.NA
