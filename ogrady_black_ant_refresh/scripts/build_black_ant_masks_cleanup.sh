#!/bin/bash -l
#SBATCH -J black_ant_masks_cleanup
#SBATCH -p cloud
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -o sbatch-cleanup-%j.out
#SBATCH -e sbatch-cleanup-%j.err
#SBATCH --exclusive

# build_black_ant_masks_cleanup.sh
# -----------------------------------------------------------------------------
# Helper wrapper that mirrors build_black_ant_masks.sh but only touches a
# targeted subject list and understands the refreshed MRI_data/WARPS layout
# (anatQQ.* outputs). Use this to backfill a targeted late-addition cohort
# without reprocessing the rest of the roster.
#
# Usage: set BLACKANT_MASK_SUBJECTS to override the default subject list (space
# separated HUMAN_EBR IDs). The script still expects the per-subject
# anat→TSE affine in MRI_data/TSE and will skip any subject missing that file or
# the SSwarper outputs.
# -----------------------------------------------------------------------------

set -euo pipefail

motip_root=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data
snheart_root=${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data

# Public GitHub copy: replace these example IDs or pass BLACKANT_MASK_SUBJECTS.
subjects_default=(
  HUMAN_EBR-100 HUMAN_EBR-101 HUMAN_EBR-102
)

if [[ -n ${BLACKANT_MASK_SUBJECTS:-} ]]; then
  read -r -a subjects <<<"${BLACKANT_MASK_SUBJECTS}"
else
  subjects=(${subjects_default[@]})
fi

echo "[INFO] Cleanup build for ${#subjects[@]} subjects"

snheart_sn_prob=${snheart_root}/ROIs/SN_probabilistic_MNI152/SN_probabilistic.nii
mni_sn_ctrl=${snheart_root}/ROIs/SN_control_OGRADY_inMNI152.nii
mni_pt_ctrl=${snheart_root}/ROIs/LC_control_OGRADY_inMNI152.nii
mni_template=${snheart_root}/MNI152_2009_template_SSW.nii.gz
out_dir=${snheart_root}/TSE/probabilistic_masks
atlas_dir=${snheart_root}/ROIs/black_ant_templates
mkdir -p "${out_dir}" "${atlas_dir}"

declare -A roi_sources=(
  [SN]="${snheart_sn_prob}"
  [SNCTRL]="${mni_sn_ctrl}"
  [PTCTRL]="${mni_pt_ctrl}"
)

require_roi() {
  local label=$1
  local path=$2
  if [[ ! -f ${path} ]]; then
    echo "[ERROR] Missing ${label} template at ${path}"
    exit 1
  fi
}

require_roi SN "${roi_sources[SN]}"
require_roi "SN control" "${roi_sources[SNCTRL]}"
require_roi "PT control" "${roi_sources[PTCTRL]}"

for label in "${!roi_sources[@]}"; do
  tgt="${atlas_dir}/${label}_probabilistic_inMNI152.nii.gz"
  echo "[INFO] Resampling ${label} ROI to ${tgt}"
  3dresample -overwrite -master "${mni_template}" -input "${roi_sources[$label]}"     -prefix "${tgt}"
  roi_sources[$label]="${tgt}"
done

find_tse() {
  local sub=$1
  local patterns=(
    "${motip_root}/TSE/zeropad_tse_${sub}.nii"
    "${motip_root}/TSE/zeropad_tse_${sub}.nii.gz"
    "${motip_root}/TSE/${sub}_TSE_zeropad.nii"
    "${motip_root}/TSE/${sub}_TSE_zeropad.nii.gz"
    "${motip_root}/TSE/${sub}_TSE.nii"
    "${motip_root}/TSE/${sub}_TSE.nii.gz"
  )
  for f in "${patterns[@]}"; do
    [[ -f ${f} ]] && { echo "${f}"; return 0; }
  done
  return 1
}

lookup_warp_stack() {
  local sub=$1
  local legacy_dir=${motip_root}/TSE/${sub}_warptoN27
  local new_dir=${motip_root}/WARPS/${sub}
  if [[ -f ${legacy_dir}/anat.un.aff.qw_WARP.nii ]]; then
    warp_field=${legacy_dir}/anat.un.aff.qw_WARP.nii
    warp_mat=${legacy_dir}/anat.un.aff.Xat.1D
    anat_master=${legacy_dir}/anat.un.aff.nii
    return 0
  fi
  if [[ -f ${new_dir}/anatQQ.${sub}_WARP.nii ]]; then
    warp_field=${new_dir}/anatQQ.${sub}_WARP.nii
    warp_mat=${new_dir}/anatQQ.${sub}.aff12.1D
    anat_master=${new_dir}/anatQQ.${sub}.nii
    return 0
  fi
  return 1
}

for entry in "${subjects[@]}"; do
  [[ -z ${entry} ]] && continue
  if [[ ${#entry} -gt 3 ]]; then
    sub=${entry: -3}
  else
    sub=${entry}
  fi
  echo "=== Cleanup subject ${entry} (${sub}) ==="
  tse_vol=$(find_tse "${sub}" || true)
  if [[ -z ${tse_vol} ]]; then
    echo "  !! No TSE volume for ${sub}; skipping"
    continue
  fi
  if ! lookup_warp_stack "${sub}"; then
    echo "  !! Missing SSwarper outputs in both legacy and WARPS paths; skipping"
    continue
  fi
  anat_to_tse=${motip_root}/TSE/${sub}_anat_toTSE_mat.aff12.1D
  if [[ ! -f ${anat_to_tse} ]]; then
    echo "  !! Missing ${anat_to_tse}; skipping"
    continue
  fi
  for label in "${!roi_sources[@]}"; do
    src_roi=${roi_sources[$label]}
    tmp_anat=$(mktemp --suffix=_${label}_${sub}.nii.gz)
    out_subject="${out_dir}/${label,,}_groupmask_in_${sub}.nii.gz"
    echo "  -> ${label}: template->anat"
    3dNwarpApply -overwrite -nwarp "INV(${warp_field}) INV(${warp_mat})"       -source "${src_roi}" -master "${anat_master}" -ainterp wsinc5 -prefix "${tmp_anat}"
    echo "     ${label}: anat->TSE (${out_subject})"
    3dAllineate -overwrite -base "${tse_vol}" -1Dmatrix_apply "${anat_to_tse}"       -input "${tmp_anat}" -prefix "${out_subject}" -final wsinc5 -float
    rm -f "${tmp_anat}"
  done

done

echo "[INFO] Cleanup build complete"
