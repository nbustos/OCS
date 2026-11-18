"""Concatenate XCP-D 0.6 per-run outputs (run-01, run-02, ...) into one file per session.

Mirrors XCP-D 0.6.0's own concatenation (``xcp_d/workflows/concatenation.py``):
- output names drop the run entity (``re.sub("_run-[0-9]+_", "_", name)``)
- CIFTI (dtseries / ptseries) and TSV files are concatenated in run order, no rescaling
- ``*_conmat.tsv`` is recomputed as ``DataFrame.corr()`` of the concatenated parcel time series
- only runs with a denoised dtseries are used (XCP-D drops failed runs)
Files that cannot be concatenated meaningfully (DCAN HDF5, QC CSVs, figures, JSON) are skipped.
Nothing is ever overwritten: every file is written to a temp name and hard-linked into place.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import uuid
from dataclasses import dataclass, field
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

RUN_RE = re.compile(r"^(sub-[^_]+_ses-[^_]+_task-[^_]+)_run-(\d+)_(.+)$")
DENOISED = "space-fsLR_den-91k_desc-denoised_bold.dtseries.nii"


def combined_name(per_run_name: str) -> str:
    """XCP-D's rule for the concatenated filename: drop ``_run-<n>_``."""
    return re.sub("_run-[0-9]+_", "_", per_run_name)


@dataclass
class SessionPlan:
    subject: str
    session: str
    func_dir: Path
    prefix: str                                  # sub-X_ses-Y_task-T
    runs: list[str]                              # run labels used, sorted
    status: str                                  # already_combined | to_combine | single_run | incomplete
    suffixes: dict[str, str] = field(default_factory=dict)   # suffix -> action (concat_cifti/concat_tsv/conmat/skip)
    note: str = ""


def _action(suffix: str) -> str:
    if suffix.endswith((".dtseries.nii", ".ptseries.nii")):
        return "concat_cifti"
    if suffix.endswith("_conmat.tsv"):
        return "conmat"
    if suffix.endswith(("_timeseries.tsv", "motion.tsv", "outliers.tsv")):
        return "concat_tsv"
    return "skip"


def plan_sessions(xcpd_dir: str | Path, task: str = "rest", out_dir: str | Path | None = None,
                  usable: set[tuple[str, str, str]] | None = None) -> list[SessionPlan]:
    """Inspect every session and decide what to do.

    already_combined: a run-less denoised dtseries exists in xcpd_dir or out_dir.
    to_combine:       >= 2 eligible runs (denoised dtseries present, and in ``usable`` if given).
    single_run:       only one eligible run; nothing to combine.
    """
    xcpd_dir = Path(xcpd_dir)
    out_dir = Path(out_dir) if out_dir else None
    plans = []
    for func in sorted(xcpd_dir.glob("sub-*/ses-*/func")):
        files = [p.name for p in func.iterdir() if f"_task-{task}_" in p.name]
        by_run: dict[str, set[str]] = {}
        prefix = None
        for name in files:
            m = RUN_RE.match(name)
            if m:
                prefix = m[1]
                by_run.setdefault(m[2], set()).add(m[3])
        if prefix is None:
            continue
        sub, ses = func.parts[-3][4:], func.parts[-2][4:]
        combined = f"{prefix}_{DENOISED}"
        exists = (func / combined).exists() or (out_dir is not None and (out_dir / func.relative_to(xcpd_dir) / combined).exists())
        eligible = sorted(r for r, s in by_run.items() if DENOISED in s
                          and (usable is None or (sub, ses, r) in usable))
        dropped = sorted(set(by_run) - set(eligible))
        note = f"runs not used: {','.join(dropped)}" if dropped else ""
        if exists:
            status = "already_combined"
        elif len(eligible) >= 2:
            status = "to_combine"
        elif len(eligible) == 1:
            status = "single_run"
        else:
            status = "incomplete"
        plan = SessionPlan(sub, ses, func, prefix, eligible, status, note=note)
        if status == "to_combine":
            common = set.intersection(*(by_run[r] for r in eligible))
            plan.suffixes = {s: _action(s) for s in sorted(common)}
            partial = set.union(*(by_run[r] for r in eligible)) - common
            if partial:
                plan.note += ("; " if plan.note else "") + f"present in only some runs (skipped): {sorted(partial)}"
        plans.append(plan)
    return plans


def _write_new(out: Path, writer) -> None:
    """Write via ``writer(tmp_path)`` then hard-link into place; never replaces an existing file.

    os.link fails if ``out`` exists, so even a file created by someone else between the
    existence check and the write is left untouched.

    Network shares (e.g. SMB) may not support hard links; then the file is created with
    O_CREAT|O_EXCL (also fails if it exists) and filled from the temp file. A partial file
    from a failed copy is removed.
    """
    tmp = out.with_name(f".tmp_{uuid.uuid4().hex}_{out.name}")
    try:
        writer(tmp)
        try:
            os.link(tmp, out)
        except FileExistsError:
            raise
        except OSError:                                   # hard links unsupported here
            fd = os.open(out, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
            try:
                with os.fdopen(fd, "wb") as dst, open(tmp, "rb") as src:
                    shutil.copyfileobj(src, dst, 16 * 1024 * 1024)
            except BaseException:
                out.unlink(missing_ok=True)               # only ever our own just-created file
                raise
    finally:
        tmp.unlink(missing_ok=True)


def _concat_cifti(paths: list[Path], out: Path) -> dict:
    imgs = [nib.load(str(p)) for p in paths]
    ax0 = [im.header.get_axis(0) for im in imgs]
    ax1 = [im.header.get_axis(1) for im in imgs]
    if any(a != ax1[0] for a in ax1[1:]):
        raise ValueError(f"Brain models / parcels differ across runs: {[p.name for p in paths]}")
    steps = {round(float(a.step), 6) for a in ax0}
    if len(steps) != 1:
        raise ValueError(f"TR differs across runs {steps}: {[p.name for p in paths]}")
    data = np.vstack([np.asarray(im.dataobj, dtype=np.float32) for im in imgs])
    series = nib.cifti2.SeriesAxis(start=0.0, step=ax0[0].step, size=data.shape[0], unit=ax0[0].unit)
    _write_new(out, lambda t: nib.Cifti2Image(data, (series, ax1[0])).to_filename(str(t)))
    return {"n_frames_per_run": [int(a.size) for a in ax0], "RepetitionTime": float(ax0[0].step)}


def _concat_tsv(paths: list[Path], out: Path) -> dict:
    dfs = [pd.read_table(p) for p in paths]
    if any(list(d.columns) != list(dfs[0].columns) for d in dfs[1:]):
        raise ValueError(f"Columns differ across runs: {[p.name for p in paths]}")
    _write_new(out, lambda t: pd.concat(dfs, axis=0).to_csv(t, sep="\t", index=False, na_rep="n/a"))
    return {"n_rows_per_run": [len(d) for d in dfs]}


def _order(item):
    """Write order: data files, then conmats (need combined timeseries), then the combined
    denoised dtseries LAST - it marks a session as combined, so an interrupted session is
    picked up again on the next run and only its missing files are written."""
    suffix, action = item
    return (suffix == DENOISED, action == "conmat")


def combine_session(plan: SessionPlan, out_root: Path, xcpd_dir: Path) -> list[dict]:
    """Write the concatenated files for one session. Existing files are never replaced.
    Returns one log row per output."""
    out_func = Path(out_root) / plan.func_dir.relative_to(xcpd_dir)
    out_func.mkdir(parents=True, exist_ok=True)
    log = []
    for suffix, action in sorted(plan.suffixes.items(), key=_order):
        if action == "skip":
            log.append(dict(file=suffix, action="skip"))
            continue
        srcs = [plan.func_dir / f"{plan.prefix}_run-{r}_{suffix}" for r in plan.runs]
        out = out_func / combined_name(srcs[0].name)
        if out.exists():
            log.append(dict(file=out.name, action="exists_kept"))
            continue
        if action == "concat_cifti":
            meta = _concat_cifti(srcs, out)
        elif action == "concat_tsv":
            meta = _concat_tsv(srcs, out)
        else:  # conmat from the combined parcel timeseries, as XCP-D does
            ts_path = out_func / combined_name(srcs[0].name).replace(
                "_measure-pearsoncorrelation_conmat.tsv", "_timeseries.tsv")
            if not ts_path.exists():
                log.append(dict(file=out.name, action="skip_no_timeseries"))
                continue
            corr = pd.read_table(ts_path).corr()
            _write_new(out, lambda t: corr.to_csv(t, sep="\t", na_rep="n/a", index_label="Node"))
            meta = {"ComputedFrom": ts_path.name}
            srcs = [ts_path]
        sidecar_path = out.with_name(out.name.split(".")[0] + ".json")
        if sidecar_path.exists():
            log.append(dict(file=sidecar_path.name, action="sidecar_exists_kept"))
        else:
            sidecar = {"Description": "Run-concatenated XCP-D 0.6 output (OCS scripts/combine_runs.py); "
                                      "concatenation follows XCP-D 0.6.0 (plain concatenation in run order).",
                       "Sources": [str(s) for s in srcs], "Runs": plan.runs, **meta}
            _write_new(sidecar_path, lambda t: t.write_text(json.dumps(sidecar, indent=2)))
        log.append(dict(file=out.name, action=action, **{k: str(v) for k, v in meta.items()}))
    return log
