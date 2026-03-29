#!/usr/bin/env bash
# View a single subject's ROI debug artifacts in Freeview.
# Usage: Scripts/open_roi_debug_freeview.sh <SUBJECT>
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <SUBJECT>" >&2
  exit 1
fi
sub=$1
inspect_root="${INSPECT_ROOT:-$(pwd)/tmp_roi_debug/tmp_subject_roi_${sub}}"
if [[ ! -d ${inspect_root} ]]; then
  alt_dir=$(find "$(pwd)/tmp" -type d -name "tmp_subject_roi_${sub}" -print -quit 2>/dev/null || true)
  if [[ -n ${alt_dir} && -d ${alt_dir} ]]; then
    inspect_root="${alt_dir}"
  else
    echo "No debug artifacts found for ${sub} (looked under ${inspect_root})" >&2
    exit 1
  fi
fi
anat="${inspect_root}/anatSS.${sub}.nii"
roi_anat="${inspect_root}/subject_roi_in_anat_final.nii.gz"
roi_template="${inspect_root}/subject_roi_in_template_final.nii.gz"
bbox="${inspect_root}/template_bbox_mask.nii.gz"
command -v freeview >/dev/null 2>&1 || { echo "freeview not found" >&2; exit 1; }
motip_root=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data
kw_mask=""
for candidate in \
  "${motip_root}/TSE/SN/SN_ROI_KW_${sub}_ZP.nii" \
  "${motip_root}/TSE/SN/SN_ROI_KW_${sub}_ZP.nii.gz" \
  "${motip_root}/TSE/SN/SN_ROI_KW_${sub}.nii" \
  "${motip_root}/TSE/SN/SN_ROI_KW_${sub}.nii.gz"
do
  if [[ -f ${candidate} ]]; then
    kw_mask="${candidate}"
    break
  fi
done
zeropad_tse=""
for candidate in \
  "${motip_root}/TSE/zeropad_tse_${sub}.nii" \
  "${motip_root}/TSE/zeropad_tse_${sub}.nii.gz" \
  "${motip_root}/TSE/${sub}_TSE_zeropad.nii" \
  "${motip_root}/TSE/${sub}_TSE_zeropad.nii.gz" \
  "${motip_root}/TSE/${sub}_TSE.nii" \
  "${motip_root}/TSE/${sub}_TSE.nii.gz"
do
  if [[ -f ${candidate} ]]; then
    zeropad_tse="${candidate}"
    break
  fi
done
cmd=(freeview -v "${anat}":grayscale=0,1000)
if [[ -n ${zeropad_tse} ]]; then
  cmd+=(-v "${zeropad_tse}":grayscale=0,1000)
fi
cmd+=(
  -v "${roi_anat}":colormap=heat:opacity=0.4
  -v "${roi_template}":colormap=jet:opacity=0.3
  -v "${bbox}":colormap=Red:opacity=0.1
)
if [[ -f ${kw_mask} ]]; then
  cmd+=(-v "${kw_mask}":colormap=Green:opacity=0.35)
fi
"${cmd[@]}" >/dev/null 2>&1 &
