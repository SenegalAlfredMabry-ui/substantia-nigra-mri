#!/bin/bash
# S.H.I.E.L.D_SN_Warp_v1.sh
# -----------------------------------------------------------------------------
# Refresh subject-specific substantia nigra (SN) ROIs by warping the hand-drawn
# MOTIP masks (TSE space) into each subject's MPRAGE and onward to the
# SSwarper template (MNI152_2009), then roll every warped mask into a
# probabilistic map. Mirrors SN_ROIs_warp_to_standard.sh but adds clearer
# logging/guardrails.
# -----------------------------------------------------------------------------

set -euo pipefail
# Fail fast so a missing warp or ROI never produces partially updated atlases.

motip_root=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data
roi_dir="${motip_root}/TSE/SN"

snheart_root=${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data
subjects_file="${snheart_root}/locations/TSE_SN_ASR_SUBLIST.txt"
if [[ ! -f ${subjects_file} ]]; then
  subjects_file="${snheart_root}/locations/TSE_subjects.txt"
fi
template="${snheart_root}/MNI152_2009_template_SSW.nii.gz"
output_dir="${snheart_root}/ROIs/SN_probabilistic_MNI152"
mkdir -p "${output_dir}"
# SNHEART tree stores every derived artifact (MNI-space per-subject masks + atlas).

echo "Subject list: ${subjects_file}"
mapfile -t subjects <"${subjects_file}"
echo "Total subjects: ${#subjects[@]}"
# Preload the roster so we can provide progress updates and reuse the list for
# cumulative statistics when merging.

for entry in "${subjects[@]}"; do
  # Subject IDs can include prefixes; the final three digits match file names.
  sub=${entry: -3}
  echo "=== ${entry} (${sub}) ==="

  sn_roi="${roi_dir}/SN_ROI_KW_${sub}.nii"
  # Hand-traced SN ROI per subject; skip if the contour never existed.
  if [[ ! -f "${sn_roi}" ]]; then
    echo "  !! Missing ROI ${sn_roi}; skipping"
    continue
  fi

  tse_to_mpr_mat="${motip_root}/TSE/${sub}_TSE_toMPRAGE_mat.aff12.1D"
  # Link the TSE-drawn ROI to anatomy so we can reuse the subject's SSwarper
  # nonlinear warp fields.
  if [[ ! -f "${tse_to_mpr_mat}" ]]; then
    anat_to_tse="${motip_root}/TSE/${sub}_anat_toTSE_mat.aff12.1D"
    if [[ ! -f "${anat_to_tse}" ]]; then
      echo "  !! Missing ${anat_to_tse}; skipping"
      continue
    fi
    echo "  Generating TSE→MPRAGE matrix"
    cat_matvec "${anat_to_tse}" -I >"${tse_to_mpr_mat}"
  fi

  sn_in_mprage="${roi_dir}/SN_inMPRAGE_${sub}.nii"
  if [[ ! -f "${sn_in_mprage}" ]]; then
    echo "  Warping ROI into MPRAGE"
    anat_base="${motip_root}/TSE/${sub}_anat.nii"
    if [[ ! -f "${anat_base}" ]]; then
      echo "  !! Missing anatomy ${anat_base}; skipping"
      continue
    fi
    3dAllineate -overwrite -base "${anat_base}" \
      -1Dmatrix_apply "${tse_to_mpr_mat}" \
      -input "${sn_roi}" \
      -prefix "${sn_in_mprage}" \
      -final NN
  fi

  warp_dir="${motip_root}/TSE/${sub}_warptoN27"
  warp_field="${warp_dir}/anat.un.aff.qw_WARP.nii"
  warp_mat="${warp_dir}/anat.un.aff.Xat.1D"
  if [[ ! -f "${warp_field}" || ! -f "${warp_mat}" ]]; then
    echo "  !! Missing SSwarper warps in ${warp_dir}; skipping"
    continue
  fi
  # These SSwarper products move the subject's anatomy into template space and
  # are inverted below to carry the ROI forward.

  sn_in_template="${output_dir}/SN_inMNI152_${sub}.nii"
  if [[ ! -f "${sn_in_template}" ]]; then
    echo "  Applying nonlinear warp to MNI152"
    3dNwarpApply -overwrite \
      -nwarp "${warp_field} ${warp_mat}" \
      -source "${sn_in_mprage}" \
      -master "${template}" \
      -prefix "${sn_in_template}"
  fi

  sn_thresh="${output_dir}/SN_inMNI152_thresh_${sub}.nii"
  echo "  Thresholding warped ROI (>0.1)"
  3dcalc -overwrite -a "${sn_in_template}" -expr 'ispositive(a-0.1)' -prefix "${sn_thresh}"

done

echo "Merging thresholded masks into probabilistic volume"
cd "${output_dir}"
if ls SN_inMNI152_thresh_*.nii >/dev/null 2>&1; then
  # Aggregate all subjects into a gcount map, then divide by the roster size so
  # the atlas stores true 0–1 probabilities.
  count=$(ls SN_inMNI152_thresh_*.nii | wc -l)
  3dmerge -overwrite -gcount -1thresh 1 -prefix SN_probabilistic_counts.nii SN_inMNI152_thresh_*.nii
  3dcalc -overwrite -a SN_probabilistic_counts.nii -expr "a/${count}" -prefix SN_probabilistic.nii
  rm -f SN_probabilistic_counts.nii
  echo "Wrote ${output_dir}/SN_probabilistic.nii (normalized by ${count} subjects)"
else
  echo "No thresholded masks found; skipping merge"
fi
