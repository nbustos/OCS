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

1. Inventory of all XCP-D run files, with a sample overview figure (subjects / sessions / runs found,
   runs missing motion files)
2. Header checks: TR, run lengths, which FD threshold XCP-D applied
3. Motion QC from the DCAN HDF5 files and run/session inclusion (runs and sessions need >= 4 min of
   usable data); the sample overview is redrawn with the counts passing motion QC
4. ROI QC for seed ROIs (coverage, signal amplitude, homogeneity, global-signal coupling); XCP-D outputs only
5. Final analysis sample (`dfs/final_sample.csv`, `dfs/run_qc.csv`)
6. One analysis input per session (`dfs/analysis_inputs.csv`): the combined-runs file by default,
   or the best single run by motion with `SESSION_INPUT = 'best_run'` (CELL 1)

7. Seed-based connectivity (`derivatives/seedfc/`, git-ignored) for every seed, plus seed ↔ parcel
   (ROI-to-ROI) connectivity with the tripartite networks (fronto-parietal control, default mode, salience;
   Schaefer 7-network parcels of XCP-D's 4S356 atlas), then group
   differences, pre -> post change and regressors (age, symptom severity) defined in `MODELS`
   (CELL 1). Phenotypes: one row per subject x timepoint in `notebooks/pheno/pheno.csv`
   (git-ignored; format in `resources/pheno_template.csv`, or set `OCS_PHENO`).

QC thresholds are set in CELL 1 (and CELL 12 for ROI). Tables go to `notebooks/dfs/`
and figures to `notebooks/figures/`; both are git-ignored. **Notebook outputs are
stripped on commit** by `nbstripout` (enabled by `setup_env.sh`), so results and
subject IDs never reach this public repo. Re-run `./setup_env.sh` after pulling if
`nbstripout --status` says it is not installed.

Reusable code lives in `src/ocs/` (`discover.py`, `motion.py`, `roi.py`, `seedfc.py`, `glm.py`, `xfm.py`).

## Combining runs (run-01 + run-02)

`scripts/combine_runs.py` concatenates each session's rest runs for sessions that
XCP-D did not already combine, following XCP-D 0.6.0's own concatenation (run entity
dropped from names; dtseries/ptseries/TSVs stacked in run order; conmats recomputed from
the combined parcel time series). DCAN HDF5 and QC CSVs are not combined; motion QC
stays per run.

```bash
python scripts/combine_runs.py                 # dry run: status of every session
python scripts/combine_runs.py --run           # write combined files into xcp_d/ (next to the runs)
python scripts/combine_runs.py --run --run-qc notebooks/dfs/run_qc.csv   # usable runs only
python scripts/combine_runs.py --test -v       # test mode: first 10 subjects, list every file
python scripts/combine_runs.py --test --run    # combine only the first 10 subjects
```

**Test mode** (`--test`, or `--test N`) processes only the first 10 subjects (sorted by ID),
writing into `xcp_d/` like a full run; a later full run skips them. `-v` lists each file and
its action.

Combined files are added to the XCP-D `func/` folders with XCP-D's own run-less names
(`--out-dir` to write elsewhere). **Nothing is ever overwritten**: each file is written to
a hidden temp name and linked into place, which fails if the target exists (on shares
without hard links, an exclusive create is used instead). Sessions XCP-D already combined
are skipped. A run-less denoised file is only treated as "already combined" if its sidecar lists
≥ 2 runs; one built from a single unlabeled run (or with run-less ReHo next to it) is reported as
`single_run_unlabeled`. The combined denoised dtseries is written last, so an interrupted session is
resumed on the next run (missing files only). Each output gets a JSON sidecar listing its
source files; a log CSV goes to `logs/` (git-ignored).

## ROI masks and template space

The CIFTI 91k subcortical voxels are in MNI152NLin6Asym (2 mm) and fMRIPrep has already registered
each subject there, so masks need no per-subject warping — but every NIfTI ROI must be in
NLin6Asym. Murty 2014 is on that grid; CIT168 is used via a copy warped from 2009cAsym
(`python scripts/warp_cit168_to_nlin6.py`, see `resources/README.md`). The NM-MRI masks (SPM DARTEL
"MNI", Cassidy 2019) fit the NLin6Asym T1 within 1 mm; their SNc/SNr/VTA/PBP/LC files carry an
x-mirrored header relative to `SN_mask`/`CC_mask` (same voxel array), so the notebook reads them with
`SN_mask`'s header (`affine_from` in `ROIS`). CELL 14 Fig 3 overlays the midbrain masks on the template.

## Data paths

The XCP-D derivatives path defaults to the Synology mount on cortex. Override
per machine with:

```bash
export OCS_XCPD_DIR=/path/to/derivatives/xcp_d
```

No participant data is tracked in this repo (see `.gitignore`).
