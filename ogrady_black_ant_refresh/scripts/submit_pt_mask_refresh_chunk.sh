#!/bin/bash -l
# submit_pt_mask_refresh_chunk.sh
# Wrapper to rebuild PT-only masks for the newly rewarped subjects.
set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
cd "$repo_root"

# Public GitHub copy: the original first-batch roster was replaced with examples.
default_subjects=(
  100 101 102
)

if [[ -z ${BLACKANT_MASK_SUBJECTS:-} ]]; then
  BLACKANT_MASK_SUBJECTS="${default_subjects[*]}"
fi
export BLACKANT_MASK_SUBJECTS

# Point the builder at the refreshed warp cache unless overridden.
export BLACKANT_SUPERVILLAIN_WARPS="${BLACKANT_SUPERVILLAIN_WARPS:-$repo_root/MRI_data/WARPS_supervillain}"

echo "[INFO] Launching PT mask refresh for subjects: $BLACKANT_MASK_SUBJECTS"
# build_black_ant_masks_v7.sh already requests 2 nodes via SBATCH directives.
bash Scripts/run_build_black_ant_masks_v7_ptonly.sh "$@"
