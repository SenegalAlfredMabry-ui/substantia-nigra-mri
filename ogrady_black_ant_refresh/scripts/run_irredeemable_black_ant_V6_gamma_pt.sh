#!/bin/bash -l
#SBATCH -J irredeemable_gamma_pt
#SBATCH -p cloud
#SBATCH -c 4
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -o sbatch-%j.out
#SBATCH -e sbatch-%j.err

set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
cd "$repo_root"

# Gamma defaults: SN prob 0.15, keep CP/PT masks unthresholded, skip PT brainstem clamp.
export IRREDEEMABLE_SN_PROB_FIXED=${IRREDEEMABLE_SN_PROB_FIXED:-0.15}
unset IRREDEEMABLE_CTRL_PROB_FIXED
unset IRREDEEMABLE_PT_PROB_FIXED
unset IRREDEEMABLE_PT_INTERSECT_BRAINSTEM

# Let MATLAB pick up all other optional knobs from the environment if set.
export MATLAB_ROOT=${MATLAB_ROOT:-"$repo_root/Software"}
if [[ ! -x "$MATLAB_ROOT/bin/matlab" ]]; then
  echo "[ERROR] MATLAB binary not found at $MATLAB_ROOT/bin/matlab" >&2
  exit 1
fi

"$MATLAB_ROOT/bin/matlab" -nodisplay -nosplash -nodesktop \
  -r "addpath('$repo_root/Scripts'); run('Scripts/irredeemable_black_ant_V6_gamma_pt.m'); exit;"
