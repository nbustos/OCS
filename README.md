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

## Data paths

The XCP-D derivatives path defaults to the Synology mount on cortex. Override
per machine with:

```bash
export OCS_XCPD_DIR=/path/to/derivatives/xcp_d
```

No participant data is tracked in this repo (see `.gitignore`).
