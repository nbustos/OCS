#!/usr/bin/env bash
# Create (or update) the `ocs` conda env, install rsx-analyzer + this repo,
# and register a Jupyter kernel.
#
# Usage:
#   ./setup_env.sh                                # create or update env
#   RSX_PATH=~/Software/XCPD ./setup_env.sh       # use a local rsX-Analyzer checkout (editable)
#   ENV_NAME=ocs-test ./setup_env.sh              # different env name
#
# Add a package: put it in environment.yml (conda-forge, or under pip:), then re-run.
set -euo pipefail

ENV_NAME="${ENV_NAME:-ocs}"
RSX_GIT="git+https://github.com/nbustos/rsX-Analyzer.git"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

command -v conda >/dev/null || { echo "conda not found on PATH"; exit 1; }
SOLVER=conda; command -v mamba >/dev/null && SOLVER=mamba

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo ">> Updating env '$ENV_NAME'"
  $SOLVER env update -n "$ENV_NAME" -f "$REPO_DIR/environment.yml" --prune
else
  echo ">> Creating env '$ENV_NAME'"
  $SOLVER env create -n "$ENV_NAME" -f "$REPO_DIR/environment.yml"
fi

RUN="conda run --no-capture-output -n $ENV_NAME"

if [[ -n "${RSX_PATH:-}" ]]; then
  echo ">> Installing rsx-analyzer (editable) from $RSX_PATH"
  $RUN pip install -e "$RSX_PATH"
else
  echo ">> Installing rsx-analyzer from GitHub"
  $RUN pip install --upgrade "$RSX_GIT"
fi

echo ">> Installing ocs package (editable)"
$RUN pip install -e "$REPO_DIR"

echo ">> Enabling nbstripout git filter (notebook outputs are never committed)"
(cd "$REPO_DIR" && $RUN nbstripout --install --attributes .gitattributes)

echo ">> Registering Jupyter kernel"
$RUN python -m ipykernel install --user --name "$ENV_NAME" --display-name "Python ($ENV_NAME)"

echo ">> Sanity check"
$RUN python - <<'PY'
import nibabel, nilearn, rsx_analyzer, ocs
print(f"nibabel {nibabel.__version__} | nilearn {nilearn.__version__} | rsx_analyzer + ocs import OK")
PY

if command -v wb_command >/dev/null; then
  echo ">> Connectome Workbench: $(wb_command -version | grep -i '^version' || echo found)"
else
  echo "!! wb_command not on PATH. Not required for the Python pipeline, but needed for"
  echo "   wb_view / CIFTI smoothing / dscalar resampling. Install from"
  echo "   https://www.humanconnectome.org/software/get-connectome-workbench"
fi

echo "Done. Activate with:  conda activate $ENV_NAME"
