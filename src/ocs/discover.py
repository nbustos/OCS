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


COMBINED_FILES = {
    "space-fsLR_den-91k_desc-denoised_bold.dtseries.nii": "dtseries",
    "space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii": "dtseries_smoothed",
}
_RUN_IN_SOURCE = re.compile(r"_run-([^_]+)_")


def find_combined(xcpd_dir: str | Path, task: str = "rest") -> pd.DataFrame:
    """One row per session that has a run-less (combined) denoised dtseries.

    Columns: subject, session, dtseries, dtseries_smoothed, kind ('combined' = built from >= 2 runs,
    'single' = one unlabeled run, 'unknown'), runs (comma-separated run labels it was built from,
    read from the sidecar: ``Runs`` for files made by scripts/combine_runs.py, else the run entities
    in XCP-D's ``Sources``; empty if unknown).
    """
    import json
    from .combine import runless_kind

    xcpd_dir = Path(xcpd_dir)
    rows = []
    for p in sorted(xcpd_dir.glob(f"sub-*/ses-*/func/sub-*_ses-*_task-{task}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii")):
        sub, ses = p.parts[-4][4:], p.parts[-3][4:]
        row = dict(subject=sub, session=ses, dtseries=str(p), kind=runless_kind(p), runs="")
        sm = p.with_name(p.name.replace("desc-denoised_", "desc-denoisedSmoothed_"))
        row["dtseries_smoothed"] = str(sm) if sm.exists() else pd.NA
        try:
            meta = json.loads(p.with_name(p.name.split(".")[0] + ".json").read_text())
            runs = meta.get("Runs") or sorted({m for s in meta.get("Sources", []) for m in _RUN_IN_SOURCE.findall(str(s))})
            row["runs"] = ",".join(runs)
        except (OSError, ValueError):
            pass
        rows.append(row)
    cols = ["subject", "session", "dtseries", "dtseries_smoothed", "kind", "runs"]
    return pd.DataFrame(rows, columns=cols)
