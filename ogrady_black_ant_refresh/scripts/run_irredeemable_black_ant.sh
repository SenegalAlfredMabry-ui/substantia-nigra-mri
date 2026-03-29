#!/bin/bash -l
#SBATCH -J irredeemable_metrics
#SBATCH -p cloud
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o sbatch-%j.out
#SBATCH -e sbatch-%j.err

set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
cd "$repo_root"

export MATLAB_ROOT=${MATLAB_ROOT:-"$repo_root/Software"}
if [[ ! -x "$MATLAB_ROOT/bin/matlab" ]]; then
  echo "[ERROR] MATLAB binary not found at $MATLAB_ROOT/bin/matlab" >&2
  exit 1
fi

"$MATLAB_ROOT/bin/matlab" -nodisplay -nosplash -nodesktop \
  -r "addpath('$repo_root/Scripts'); run('Scripts/irredeemable_black_ant_V1.m'); exit;"
