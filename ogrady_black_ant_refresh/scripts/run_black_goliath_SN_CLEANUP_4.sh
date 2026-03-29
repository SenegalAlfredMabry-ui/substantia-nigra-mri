#!/bin/bash -l
#SBATCH -J black_goliath_cleanup
#SBATCH -p cloud
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -o sbatch-black-goliath-%j.out
#SBATCH -e sbatch-black-goliath-%j.err
#SBATCH --time=12:00:00

set -euo pipefail

MATLAB_BIN="$(pwd)/Software/bin/matlab"
SCRIPT="black_goliath_SN_CLEANUP_4.m"

if [[ ! -x $MATLAB_BIN ]]; then
  echo "MATLAB binary not found at $MATLAB_BIN" >&2
  exit 1
fi

$MATLAB_BIN -nodisplay -nosplash -nodesktop <<'MAT'
addpath('${SNHEART_ROOT:-/path/to/SNHEART}/Scripts');
run('black_goliath_SN_CLEANUP_4.m');
exit;
MAT
