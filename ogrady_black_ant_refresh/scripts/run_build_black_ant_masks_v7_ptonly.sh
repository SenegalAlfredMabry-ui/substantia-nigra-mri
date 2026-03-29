#!/bin/bash -l
set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
cd "$repo_root"

export BLACKANT_MASK_LABELS="PTCTRL"

exec bash Scripts/run_build_black_ant_masks_v7.sh "$@"
