# OCS

Resting-state analyses of the OCS (marsonc) XCP-D outputs: seed-based
connectivity on fsLR 91k CIFTI (`desc-denoised[Smoothed]_bold.dtseries.nii`),
building on [rsX-Analyzer](https://github.com/nbustos/rsX-Analyzer).

## Setup (cortex or laptop)

```bash
git clone https://github.com/nbustos/OCS.git
cd OCS
./setup_env.sh            # creates/updates the `ocs` conda env + Jupyter kernel
conda activate ocs
```

- **Add a package:** add it to `environment.yml` (conda-forge, or under `pip:`),
  commit, then `git pull && ./setup_env.sh` on each machine.
- **Develop rsx alongside:** `RSX_PATH=/path/to/rsX-Analyzer ./setup_env.sh`
  installs it editable instead of from GitHub.
- **Connectome Workbench** (`wb_command`) is not on conda-forge; install it
  separately if needed. The setup script reports whether it is on `PATH`.

## Main notebook

`notebooks/OCS_REST_main.ipynb` (kernel **Python (ocs)**):

1. Inventory of all XCP-D run files
2. Header checks: TR, run lengths, which FD threshold XCP-D applied
3. Motion QC from the DCAN HDF5 files and run/session inclusion
4. ROI QC for seed ROIs (coverage, signal amplitude, homogeneity, global-signal coupling, tSNR)
5. Final analysis sample (`dfs/final_sample.csv`, `dfs/run_qc.csv`)

QC thresholds are set in CELL 1 (and CELL 10 for ROI). Tables go to `notebooks/dfs/`
and figures to `notebooks/figures/`; both are git-ignored. **Notebook outputs are
stripped on commit** by `nbstripout` (enabled by `setup_env.sh`), so results and
subject IDs never reach this public repo. Re-run `./setup_env.sh` after pulling if
`nbstripout --status` says it is not installed.

Reusable code lives in `src/ocs/` (`discover.py`, `motion.py`, `roi.py`).

## Data paths

The XCP-D derivatives path defaults to the Synology mount on cortex. Override
per machine with:

```bash
export OCS_XCPD_DIR=/path/to/derivatives/xcp_d
```

No participant data is tracked in this repo (see `.gitignore`).
