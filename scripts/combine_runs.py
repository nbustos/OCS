#!/usr/bin/env python
"""Combine XCP-D rest run-01 + run-02 (+ ...) into one set of files per session.

Only sessions WITHOUT an existing combined output are processed. Default is a dry run.

Examples (on cortex, `conda activate ocs`):
    python scripts/combine_runs.py                       # report only, writes nothing
    python scripts/combine_runs.py --run                 # write to derivatives/xcp_d_combined
    python scripts/combine_runs.py --run --run-qc notebooks/dfs/run_qc.csv   # only QC-usable runs
    python scripts/combine_runs.py --run --subject marsonc1001 --subject marsonc1003

Outputs mirror XCP-D 0.6.0 concatenation: run entity dropped from names, CIFTI/TSV files
concatenated in run order, conmats recomputed from the combined parcel time series.
Existing files are never overwritten unless --overwrite is given. A JSON sidecar listing the
source files is written next to each output, and a log CSV is written to the output root.
"""
import argparse
import datetime as dt
import sys
from pathlib import Path

import pandas as pd

from ocs import XCPD_DIR
from ocs.combine import combine_session, plan_sessions


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xcpd-dir", type=Path, default=XCPD_DIR, help=f"XCP-D derivatives (default: {XCPD_DIR})")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="where combined files go (default: <xcpd-dir>/../xcp_d_combined)")
    ap.add_argument("--task", default="rest")
    ap.add_argument("--subject", action="append", help="limit to subject label(s), without 'sub-'")
    ap.add_argument("--run-qc", type=Path, help="run_qc.csv from the main notebook; only runs with usable==True are combined")
    ap.add_argument("--run", action="store_true", help="actually write files (default: dry run)")
    ap.add_argument("--overwrite", action="store_true", help="replace combined files that already exist in --out-dir")
    a = ap.parse_args(argv)

    out_dir = a.out_dir or a.xcpd_dir.parent / "xcp_d_combined"
    usable = None
    if a.run_qc:
        q = pd.read_csv(a.run_qc, dtype={"subject": str, "session": str, "run": str})
        usable = {(r.subject, r.session, r.run) for r in q.itertuples() if bool(r.usable)}
        print(f"Using {len(usable)} usable runs from {a.run_qc}")

    plans = plan_sessions(a.xcpd_dir, a.task, out_dir, usable)
    if a.subject:
        plans = [p for p in plans if p.subject in set(a.subject)]
    if not plans:
        sys.exit(f"No task-{a.task} sessions found under {a.xcpd_dir}")

    report = pd.DataFrame([dict(subject=p.subject, session=p.session, status=p.status,
                                runs=",".join(p.runs), n_files=sum(v != "skip" for v in p.suffixes.values()),
                                note=p.note) for p in plans])
    print(f"\nXCP-D dir : {a.xcpd_dir}\nOutput dir: {out_dir}\n")
    print(report.to_string(index=False))
    print("\n" + report.status.value_counts().to_string())

    todo = [p for p in plans if p.status == "to_combine"]
    if todo:
        actions = pd.Series([v for p in todo for v in p.suffixes.values()]).value_counts()
        print(f"\nPer session to combine ({len(todo)} sessions), file actions:\n" + actions.to_string())
    if not a.run:
        print("\nDry run - nothing written. Re-run with --run to combine.")
        return

    rows = []
    for i, p in enumerate(todo, 1):
        print(f"[{i}/{len(todo)}] sub-{p.subject} ses-{p.session} runs {','.join(p.runs)}", flush=True)
        try:
            for r in combine_session(p, out_dir, a.xcpd_dir, overwrite=a.overwrite):
                rows.append(dict(subject=p.subject, session=p.session, runs=",".join(p.runs), **r))
        except Exception as e:  # keep going; report at the end
            rows.append(dict(subject=p.subject, session=p.session, runs=",".join(p.runs), file="", action="ERROR", error=str(e)))
            print(f"   ERROR: {e}")
    log = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / f"combine_runs_log_{dt.datetime.now():%Y%m%d_%H%M%S}.csv"
    log.to_csv(log_path, index=False)
    print(f"\nDone. {int((log.action == 'ERROR').sum())} errors. Log: {log_path}")
    print(log.action.value_counts().to_string())


if __name__ == "__main__":
    main()
