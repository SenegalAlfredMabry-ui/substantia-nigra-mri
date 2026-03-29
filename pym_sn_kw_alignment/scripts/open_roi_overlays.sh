#!/usr/bin/env bash
# Launch AFNI with the antpunk ROI underlay/overlay pairs for quick inspection.
# Usage: Scripts/open_roi_overlays.sh [SUBJECT ...]
# With no arguments the public-copy example roster is used.
set -euo pipefail

roi_root="${ROI_ROOT:-$(pwd)/MRI_data/WARPS_antpunk/roi_inspection_20260210}"
default_subs=(100 101 102)
if [[ $# -gt 0 ]]; then
  subjects=($@)
else
  subjects=(${default_subs[@]})
fi

command -v afni >/dev/null 2>&1 || { echo "afni not found in PATH" >&2; exit 1; }

for sub in "${subjects[@]}"; do
  anat="${roi_root}/${sub}_anatQQ.nii"
  overlay="${roi_root}/${sub}_subject_roi_in_template.nii.gz"
  if [[ ! -f ${anat} || ! -f ${overlay} ]]; then
    echo "Skipping ${sub}: missing ${anat} or ${overlay}" >&2
    continue
  fi
  echo "Launching AFNI for ${sub} (underlay=${anat}, overlay=${overlay})"
  AFNI_NIFTI_TYPE_WARN=NO afni "${anat}" "${overlay}" >/dev/null 2>&1 &
  sleep 1
done

echo "All requested AFNI sessions launched." 
