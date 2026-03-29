#!/bin/bash -l
#SBATCH -J black_ant_masks_v5
#SBATCH -p cloud
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -c 4
#SBATCH -o sbatch-%j.out
#SBATCH -e sbatch-%j.err
#SBATCH --exclusive
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=FAIL

# build_black_ant_masks_v5.sh
# -----------------------------------------------------------------------------
# PURPOSE
#   Regenerate the MOTIP subject-space neuromelanin probability masks that feed
#   every Black ANT metrics run. The collaborator asked for a reproducible
#   "template (N27) → MNI152 → subject anat → subject TSE" chain so the masks
#   match the SSwarper transforms bundled with ANT MAN.
#
# WORKFLOW OVERVIEW
#   1. Locate the group ROIs (SN, SN-control, PT-control). SN now prefers the
#      normalized MNI152 atlas created by S.H.I.E.L.D., while the control
#      atlases still fall back to the original MOTIP N27 volumes until
#      replacements exist.
#   2. Resample each ROI into the SSwarper reference template
#      (MNI152_2009_template_SSW) once so every subject starts from the same
#      canonical grid.
#   3. For every subject, walk the transforms backwards: use the inverse
#      SSwarper warp (template→anat) followed by the saved anat→TSE affine to
#      land each ROI in the native TSE space. This mirrors how ANT MAN pushes
#      brains into template space but reverses the direction so the masks align
#      with each subject’s raw/zeropad images.
#   4. After each subject’s ROIs land in TSE space, subtract the SN mask from
#      the SN-control and PT-control masks so the backgrounds never include the
#      target voxels, then stage the accompanying brainstem/DC masks the
#      Irredeemable brainstem CNR depends on.
#   5. Write the subject-space masks under 77_SNHEART (MRI_data/TSE/
#      probabilistic_masks_v5) so MATLAB/R jobs can read local copies without
#      mounting 33_MOTIP2018.
#
# REQUIREMENTS
#   - AFNI binaries for 3dAllineate + 3dNwarpApply (same stack used by ANT MAN).
#   - The historical MOTIP directory mounted at ${MOTIP_ROOT}/...
#   - Prior SSwarper outputs (anat.un.aff* + warp field) and the anat→TSE
#     matrices produced by the zeropad preprocessing stage.
# -----------------------------------------------------------------------------

set -euo pipefail
trap 'echo "[ERROR] Command ${BASH_COMMAND} failed at line ${LINENO}" >&2' ERR
set -x

# Ensure AFNI binaries (3dAllineate/3dNwarpApply/3dcalc) are available on the
# compute nodes. Loading the module mirrors how ANT MAN + Visioneer scripts
# prime the AFNI toolchain before launching.
if command -v module >/dev/null 2>&1; then
  module load afni
fi
if ! command -v 3dAllineate >/dev/null 2>&1; then
  echo "[ERROR] AFNI/3dAllineate not found on PATH. Did module load afni succeed?" >&2
  exit 1
fi

# Canonical roots for MOTIP sources (zeropad TSE volumes + warps) and SNHEART
# targets (where we stash the warped masks). Keeping paths explicit avoids
# surprises when the script is invoked from arbitrary working directories.
motip_root=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data
snheart_root=${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data

# Subject roster lives alongside SNHEART so we can curate it without touching
# the legacy MOTIP tree. Use TSE_SN_ASR_SUBLIST when present to restrict to the
# 108 traced subjects feeding the atlas, otherwise fall back to TSE_subjects.
subjects_file=${snheart_root}/locations/TSE_SN_ASR_SUBLIST.txt
if [[ ! -f ${subjects_file} ]]; then
  subjects_file=${snheart_root}/locations/TSE_subjects.txt
fi

# Preferred SN atlas lives in the refreshed SNHEART directory (already in
# MNI152 space and normalized to probabilities). Control atlases also ship
# inside SNHEART so the entire workflow stays on the same template.
snheart_sn_prob=${snheart_root}/ROIs/SN_probabilistic_MNI152/SN_probabilistic.nii
# Prefer the O'Grady SN-control + locus coeruleus (PT-control) atlases that
# already live in MNI152 space to match the refreshed SN probabilistic inputs.
mni_sn_ctrl=${snheart_root}/ROIs/SN_control_OGRADY_inMNI152.nii
mni_pt_ctrl=${snheart_root}/ROIs/LC_control_OGRADY_inMNI152.nii

# SSwarper template (MNI152) so we match the exact reference grid that powered
# the original ANT MAN runs. Resampling into this space makes the per-subject
# inverse warp straightforward (template<->anat pairs already exist).
mni_template=${snheart_root}/MNI152_2009_template_SSW.nii.gz

# Destination directory (subject-space probabilistic masks) lives inside the
# SNHEART TSE tree so downstream jobs can read them next to the raw volumes.
# Allow BLACKANT_OUTPUT_DIR/BLACKANT_ATLAS_DIR overrides so experiments do not
# clobber the production masks and templates.
out_dir=${BLACKANT_OUTPUT_DIR:-${snheart_root}/TSE/probabilistic_masks_v5}
mkdir -p "${out_dir}"

# Keep resampled MNI-space atlases in a dedicated folder under ROIs so we know
# which files describe the shared template versus the per-subject masks.
atlas_dir=${BLACKANT_ATLAS_DIR:-${snheart_root}/ROIs/black_ant_templates_v5}
mkdir -p "${atlas_dir}"

# Manifest captures mask/background availability for Thunderbolts analysis.
manifest_path=${BLACKANT_MASK_MANIFEST:-${out_dir}/black_ant_mask_manifest_v5.tsv}
printf "subject\tsn_voxels\tsnctrl_voxels\tptctrl_voxels\tbrainstem_mask\tbrainstem4v_mask\tdc_mask\tnotes\n" >"${manifest_path}"

# Optional control-mask probability cutoff mirroring the v3 tightening. When
# unset we leave the floating-point probabilities untouched.
ctrl_prob_threshold=${BLACKANT_CTRL_PROB_THRESHOLD:-}
if [[ -n ${ctrl_prob_threshold} ]]; then
  echo "[INFO] Control mask threshold enabled: >= ${ctrl_prob_threshold}"
fi
# Allow experiments that keep SN voxels inside the control masks by setting
# BLACKANT_SUBTRACT_SN=0. Default = 1 (always subtract SN from controls).
subtract_sn=${BLACKANT_SUBTRACT_SN:-1}
if [[ ${subtract_sn} -eq 0 ]]; then
  echo "[INFO] SN-control/PT-control subtraction disabled for this build"
fi

# Stage DC/brainstem masks that Irredeemable uses for the brainstem CNR.
brainstem_dir=${BLACKANT_BRAINSTEM_DIR:-${snheart_root}/TSE/brainstem_masks_v5}
mkdir -p "${brainstem_dir}"
brainstem_roots=(
  "${brainstem_dir}"
  "${motip_root}/TSE/SN"
  "${motip_root}/TSE"
  "${snheart_root}/ants_processing/upsampled/TSE"
)

first_existing_mask() {
  local base=$1
  for ext in .nii .nii.gz; do
    if [[ -f ${base}${ext} ]]; then
      echo "${base}${ext}"
      return 0
    fi
  done
  return 1
}

mask_voxel_count() {
  local mask_path=$1
  if [[ -z ${mask_path} || ! -f ${mask_path} ]]; then
    echo 0
    return 0
  fi
  local count status
  set +e
  count=$(3dBrickStat -count -nonzero "${mask_path}" 2>/dev/null)
  status=$?
  set -e
  if [[ ${status} -ne 0 ]] || [[ -z ${count} ]]; then
    echo 0
  else
    awk 'END{print $NF+0}' <<<"${count}"
  fi
}

# Scratch space for intermediate anat-space volumes. Ensures we do not leave
# clutter in the repo if the script terminates early.
tmp_dir=$(mktemp -d)
cleanup() {
  rm -rf "${tmp_dir}"
}
trap cleanup EXIT

has_dest_mask() {
  local label=$1
  local sub=$2
  local base="${brainstem_dir}/${label}_${sub}"
  for ext in .nii .nii.gz; do
    if [[ -f ${base}${ext} ]]; then
      return 0
    fi
  done
  return 1
}

copy_mask_from_sources() {
  local label=$1
  local sub=$2
  if ! has_dest_mask "${label}" "${sub}"; then
    for root in "${brainstem_roots[@]}"; do
      [[ -z ${root} || ! -d ${root} ]] && continue
      for ext in .nii .nii.gz; do
        local candidate="${root}/${label}_${sub}${ext}"
        if [[ -f ${candidate} ]]; then
          local dest="${brainstem_dir}/${label}_${sub}${ext}"
          cp -f "${candidate}" "${dest}"
          echo "  -> Staged ${label}_${sub}${ext} from ${root}"
          return 0
        fi
      done
    done
    echo "  !! Missing ${label}_${sub} under ${brainstem_roots[*]}"
    return 1
  fi
}

stage_brainstem_masks() {
  local sub=$1
  copy_mask_from_sources "DC_mask_inTSE" "${sub}" || true
  copy_mask_from_sources "brainstem_mask_inTSE" "${sub}" || true
  copy_mask_from_sources "brainstem_4v_mask_inTSE" "${sub}" || true
}

record_manifest_entry() {
  local sub=$1
  local sn_mask snctrl_mask ptctrl_mask brainstem_mask brainstem4v_mask dc_mask
  sn_mask="$(first_existing_mask "${out_dir}/sn_groupmask_in_${sub}")"
  snctrl_mask="$(first_existing_mask "${out_dir}/snctrl_groupmask_in_${sub}")"
  ptctrl_mask="$(first_existing_mask "${out_dir}/ptctrl_groupmask_in_${sub}")"
  brainstem_mask="$(first_existing_mask "${brainstem_dir}/brainstem_mask_inTSE_${sub}")"
  brainstem4v_mask="$(first_existing_mask "${brainstem_dir}/brainstem_4v_mask_inTSE_${sub}")"
  dc_mask="$(first_existing_mask "${brainstem_dir}/DC_mask_inTSE_${sub}")"

  local sn_voxels snctrl_voxels ptctrl_voxels notes=()
  sn_voxels=$(mask_voxel_count "${sn_mask}")
  snctrl_voxels=$(mask_voxel_count "${snctrl_mask}")
  ptctrl_voxels=$(mask_voxel_count "${ptctrl_mask}")
  [[ -z ${brainstem_mask} ]] && notes+=('brainstem_missing')
  [[ -z ${brainstem4v_mask} ]] && notes+=('brainstem4v_missing')
  [[ -z ${dc_mask} ]] && notes+=('dc_missing')
  local note_str
  if [[ ${#notes[@]} -eq 0 ]]; then
    note_str="ok"
  else
    note_str=$(IFS=,; echo "${notes[*]}")
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sub}" "${sn_voxels}" "${snctrl_voxels}" "${ptctrl_voxels}" \
    "${brainstem_mask:-missing}" "${brainstem4v_mask:-missing}" "${dc_mask:-missing}" \
    "${note_str}" >>"${manifest_path}"
}

# Map ROI labels to their source files. Associative array lets us track which
# templates survived the existence checks and rewrites the path once they are
# resampled into MNI space.
declare -A roi_sources=(
  [SN]="${snheart_sn_prob}"
  [SNCTRL]="${mni_sn_ctrl}"
  [PTCTRL]="${mni_pt_ctrl}"
)

require_roi() {
  local label=$1
  local path=$2
  if [[ ! -f ${path} ]]; then
    echo "[ERROR] Missing ${label} template at ${path}. Aborting mask build."
    exit 1
  fi
}

require_roi "SN" "${roi_sources[SN]}"
require_roi "SN control" "${roi_sources[SNCTRL]}"
require_roi "PT control" "${roi_sources[PTCTRL]}"

for label in "${!roi_sources[@]}"; do
  src=${roi_sources[$label]}
  if [[ -f ${src} ]]; then
    # Resample the atlas ROI into the SSwarper template once. We use wsinc5 to
    # preserve soft probability edges and match ANT MAN interpolation choices.
    tgt=${atlas_dir}/${label}_probabilistic_inMNI152.nii.gz
    echo "[INFO] Resampling ${label} ROI from N27 -> MNI152 (${tgt})"
    3dresample -overwrite -master "${mni_template}" -input "${src}" \
      -prefix "${tgt}"
    roi_sources[$label]="${tgt}"
  else
    # Missing template volumes frequently happen for control masks. Flag them
    # loudly so the operator knows the downstream per-subject loop will skip
    # that ROI entirely.
    echo "[WARN] Missing ${label} template ${src}; skipping"
    unset "roi_sources[$label]"
  fi
done

if [[ -n ${BLACKANT_MASK_SUBJECTS:-} ]]; then
  # Optional environment variable to constrain the workload to a curated list.
  # Useful when debugging individual IDs or rebuilding a partial cohort.
  read -r -a subjects <<<"${BLACKANT_MASK_SUBJECTS}"
  echo "[INFO] Restricting to explicit subjects: ${subjects[*]}"
else
  # Default: run the full zeropad roster recorded in 33_MOTIP2018.
  mapfile -t subjects <"${subjects_file}"
fi
echo "[INFO] Processing ${#subjects[@]} MOTIP subjects"

find_tse() {
  local sub=$1
  # Zeropad TSE filenames evolved over time; probe the common permutations so
  # we can land on the right volume without manual symlinks.
  local patterns=(
    "${motip_root}/TSE/zeropad_tse_${sub}.nii"
    "${motip_root}/TSE/zeropad_tse_${sub}.nii.gz"
    "${motip_root}/TSE/${sub}_TSE_zeropad.nii"
    "${motip_root}/TSE/${sub}_TSE_zeropad.nii.gz"
    "${motip_root}/TSE/${sub}_TSE.nii"
    "${motip_root}/TSE/${sub}_TSE.nii.gz"
  )
  for f in "${patterns[@]}"; do
    if [[ -f ${f} ]]; then
      echo "${f}"
      return 0
    fi
  done
  return 1
}

lookup_warp_stack() {
  local sub=$1
  local legacy_dir=${motip_root}/TSE/${sub}_warptoN27
  local new_dir=${motip_root}/WARPS/${sub}
  warp_source=""
  if [[ -f ${legacy_dir}/anat.un.aff.qw_WARP.nii ]]; then
    warp_field=${legacy_dir}/anat.un.aff.qw_WARP.nii
    warp_mat=${legacy_dir}/anat.un.aff.Xat.1D
    anat_master=${legacy_dir}/anat.un.aff.nii
    warp_source=${legacy_dir}
    return 0
  fi
  if [[ -f ${new_dir}/anatQQ.${sub}_WARP.nii ]]; then
    warp_field=${new_dir}/anatQQ.${sub}_WARP.nii
    warp_mat=${new_dir}/anatQQ.${sub}.aff12.1D
    anat_master=${new_dir}/anatQQ.${sub}.nii
    warp_source=${new_dir}
    return 0
  fi
  return 1
}

for entry in "${subjects[@]}"; do
  # Roster sometimes lists 6-digit IDs (e.g., 202128) where the last three
  # digits are the canonical MOTIP participant number. Normalize here so the
  # downstream file lookups succeed.
  if [[ ${#entry} -gt 3 ]]; then
    sub=${entry: -3}
  else
    sub=${entry}
  fi
  echo "=== Subject ${entry} (${sub}) ==="
  tse_vol=$(find_tse "${sub}" || true)
  if [[ -z ${tse_vol} ]]; then
    echo "  !! No TSE volume for ${sub}; skipping"
    continue
  fi
  if ! lookup_warp_stack "${sub}"; then
    echo "  !! Missing SSwarper outputs for ${sub}; skipping"
    continue
  fi
  anat_to_tse=${motip_root}/TSE/${sub}_anat_toTSE_mat.aff12.1D
  if [[ ! -f ${anat_to_tse} ]]; then
    echo "  !! Missing ${anat_to_tse}; skipping"
    continue
  fi
  [[ -n ${warp_source:-} ]] && echo "  -> Using SSwarper stack from ${warp_source}"

  sn_mask_path="${out_dir}/sn_groupmask_in_${sub}.nii.gz"
  for label in "${!roi_sources[@]}"; do
    src_roi=${roi_sources[$label]}
    [[ -f ${src_roi} ]] || continue
    tmp_anat="${tmp_dir}/${label}_${sub}_anat.nii.gz"
    out_subject="${out_dir}/${label,,}_groupmask_in_${sub}.nii.gz"
    # Step 1: move the MNI-space ROI back into subject anatomical space by
    # applying the inverse of the SSwarper transforms. Using the anat.un.aff
    # volume as master ensures we preserve the SSwarper resampling grid.
    echo "  -> ${label}: template->anat"
    3dNwarpApply -overwrite -nwarp "INV(${warp_field}) INV(${warp_mat})" \
      -source "${src_roi}" -master "${anat_master}" -ainterp wsinc5 -prefix "${tmp_anat}"
    # Step 2: apply the precomputed anat→TSE affine so the ROI lands in the
    # native zeropad TSE space. wsinc5 matches the upstream interpolation and
    # keeps the probabilistic intensities stable.
    echo "     ${label}: anat->TSE (${out_subject})"
    3dAllineate -overwrite -base "${tse_vol}" -1Dmatrix_apply "${anat_to_tse}" \
      -cmass -input "${tmp_anat}" -prefix "${out_subject}" -final wsinc5 -float
    if [[ ${label} == "SNCTRL" || ${label} == "PTCTRL" ]]; then
      if [[ ${subtract_sn} -eq 0 ]]; then
        echo "     ${label}: skipping SN overlap subtraction (BLACKANT_SUBTRACT_SN=0)"
      else
        if [[ -f ${sn_mask_path} ]]; then
          tmp_ctrl="${tmp_dir}/${label,,}_${sub}_tmp.nii.gz"
          mv "${out_subject}" "${tmp_ctrl}"
          echo "     ${label}: removing SN overlap"
          3dcalc -overwrite -a "${tmp_ctrl}" -b "${sn_mask_path}" \
            -expr 'a*(1-step(b))' -prefix "${out_subject}"
        else
          echo "     [WARN] SN mask missing for ${sub}; cannot subtract overlap"
        fi
      fi
    fi
    if [[ ${label} == "SNCTRL" ]]; then
      if [[ -n ${ctrl_prob_threshold} ]]; then
        thr_out="${out_dir}/snctrl_thresholdmask_in_${sub}.nii.gz"
        echo "     ${label}: threshold >= ${ctrl_prob_threshold} (${thr_out})"
        3dcalc -overwrite -a "${out_subject}" \
          -expr "step(a-${ctrl_prob_threshold})" -prefix "${thr_out}"
      fi
    fi
  done

  stage_brainstem_masks "${sub}"
  record_manifest_entry "${sub}"

done

echo "[INFO] Outputs written to ${out_dir}"
echo "[INFO] Brainstem/DC masks staged under ${brainstem_dir}"
