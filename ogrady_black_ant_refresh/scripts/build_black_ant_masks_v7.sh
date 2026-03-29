#!/bin/bash -l
#SBATCH -J black_ant_masks_v7
#SBATCH -p cloud
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -c 4
#SBATCH -o sbatch-%j.out
#SBATCH -e sbatch-%j.err
#SBATCH --exclusive
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=FAIL

# build_black_ant_masks_v7.sh
# -----------------------------------------------------------------------------
# PURPOSE
#   Refresh the MOTIP subject-space neuromelanin probability masks that feed
#   the Black ANT metrics pipeline. V7 keeps the validated V6 workflow (template
#   chain, AFNI guards, manifest logging, staged DC/brainstem masks) but stops
#   subtracting the SN mask from the PT-control background so that ROI remains
#   identical to the source atlas.
#
# WORKFLOW OVERVIEW
#   1. Locate the group ROIs (SN, SN-control, PT-control) plus the refreshed
#      SNHEART probabilistic atlas.
#   2. Resample each ROI into the SSwarper template once (MNI152_2009_SSW) so
#      every subject uses the same canonical grid.
#   3. For each subject: invert the SSwarper transforms to land the templates
#      in anatomical space, then apply the anat→TSE affine to reach native TSE.
#   4. Subtract the SN mask from the control masks (unless explicitly disabled)
#      and emit optional thresholded masks when requested.
#   5. Stage subject-specific brainstem/DC masks so the Irredeemable CNR loops
#      have local inputs, then append an entry to the manifest describing the
#      files we wrote. Voxel counts are intentionally omitted in V6.
#
# REQUIREMENTS
#   - AFNI binaries for 3dAllineate + 3dNwarpApply (module load afni).
#   - Historical MOTIP directory mounted at ${MOTIP_ROOT}/...
#   - Prior SSwarper outputs (anat.un.aff*/anatQQ.*) and anat→TSE matrices.
# -----------------------------------------------------------------------------

set -euo pipefail
trap 'echo "[ERROR] Command ${BASH_COMMAND} failed at line ${LINENO}" >&2' ERR

if command -v module >/dev/null 2>&1; then
  module load afni
fi
if ! command -v 3dAllineate >/dev/null 2>&1; then
  echo "[ERROR] AFNI/3dAllineate not on PATH; did module load afni run?" >&2
  exit 1
fi

motip_root=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data
snheart_root=${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data
snheart_repo=${SNHEART_ROOT:-/path/to/SNHEART}
wallet_warp_root=${BLACKANT_SUPERVILLAIN_WARPS:-${snheart_repo}/MRI_data/WARPS_supervillain}

subjects_file=${snheart_root}/locations/TSE_SN_ASR_SUBLIST.txt
if [[ ! -f ${subjects_file} ]]; then
  subjects_file=${snheart_root}/locations/TSE_subjects.txt
fi

snheart_sn_prob=${snheart_root}/ROIs/SN_probabilistic_MNI152/SN_probabilistic.nii
mni_sn_ctrl=${snheart_root}/ROIs/SN_control_OGRADY_inMNI152.nii
mni_pt_ctrl=${snheart_root}/ROIs/LC_control_OGRADY_inMNI152.nii
mni_template=${snheart_root}/MNI152_2009_template_SSW.nii.gz

ctrl_prob_threshold=${BLACKANT_CTRL_PROB_THRESHOLD:-}
if [[ -n ${ctrl_prob_threshold} ]]; then
  echo "[INFO] Control threshold enabled: >= ${ctrl_prob_threshold}"
fi

out_dir=${BLACKANT_OUTPUT_DIR:-${snheart_root}/TSE/probabilistic_masks_v7}
mkdir -p "${out_dir}"
atlas_dir=${BLACKANT_ATLAS_DIR:-${snheart_root}/ROIs/black_ant_templates_v7}
mkdir -p "${atlas_dir}"
manifest_path=${BLACKANT_MASK_MANIFEST:-${out_dir}/black_ant_mask_manifest_v7.tsv}
printf "subject\tsn_voxels\tsnctrl_voxels\tptctrl_voxels\tbrainstem_mask\tbrainstem4v_mask\tdc_mask\tnotes\n" >"${manifest_path}"

subtract_sn=${BLACKANT_SUBTRACT_SN:-1}
if [[ ${subtract_sn} -eq 0 ]]; then
  echo "[INFO] SN subtraction disabled (BLACKANT_SUBTRACT_SN=0)"
fi

sn_sub_threshold=${BLACKANT_SN_SUB_THRESHOLD:-0.15}
echo "[INFO] SN subtraction threshold >= ${sn_sub_threshold} (BLACKANT_SN_SUB_THRESHOLD)"

brainstem_dir=${BLACKANT_BRAINSTEM_DIR:-${snheart_root}/TSE/brainstem_masks_v7}
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
  local output
  if ! output=$(3dBrickStat -count -nonzero "${mask_path}" 2>/dev/null); then
    echo 0
    return 0
  fi
  awk 'END{print $NF+0}' <<<"${output}"
}

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
  if has_dest_mask "${label}" "${sub}"; then
    return 0
  fi
  for root in "${brainstem_roots[@]}"; do
    [[ -z ${root} || ! -d ${root} ]] && continue
    for ext in .nii .nii.gz; do
      local candidate="${root}/${label}_${sub}${ext}"
      if [[ -f ${candidate} ]]; then
        cp -f "${candidate}" "${brainstem_dir}/${label}_${sub}${ext}"
        echo "  -> Staged ${label}_${sub}${ext} from ${root}"
        return 0
      fi
    done
  done
  echo "  !! Missing ${label}_${sub} under ${brainstem_roots[*]}"
  return 1
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

  local notes=()
  [[ -z ${sn_mask} ]] && notes+=(sn_missing)
  [[ -z ${snctrl_mask} ]] && notes+=(snctrl_missing)
  [[ -z ${ptctrl_mask} ]] && notes+=(ptctrl_missing)
  [[ -z ${brainstem_mask} ]] && notes+=(brainstem_missing)
  [[ -z ${brainstem4v_mask} ]] && notes+=(brainstem4v_missing)
  [[ -z ${dc_mask} ]] && notes+=(dc_missing)

  local sn_voxels snctrl_voxels ptctrl_voxels
  sn_voxels=$(mask_voxel_count "${sn_mask}")
  snctrl_voxels=$(mask_voxel_count "${snctrl_mask}")
  ptctrl_voxels=$(mask_voxel_count "${ptctrl_mask}")

  local note_str
  if [[ ${#notes[@]} -eq 0 ]]; then
    note_str="ok"
  else
    note_str=$(IFS=,; echo "${notes[*]}")
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sub}" "${sn_voxels}" "${snctrl_voxels}" "${ptctrl_voxels}" \
    "${brainstem_mask:-missing}" "${brainstem4v_mask:-missing}" \
    "${dc_mask:-missing}" "${note_str}" >>"${manifest_path}"
}

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
  local wallet_dir=${wallet_warp_root}/${sub}
  warp_source=""
  # Prefer the refreshed supervillain cache before older warp stacks.
  if [[ -f ${wallet_dir}/anatQQ.${sub}_WARP.nii ]]; then
    warp_field=${wallet_dir}/anatQQ.${sub}_WARP.nii
    warp_mat=${wallet_dir}/anatQQ.${sub}.aff12.1D
    anat_master=${wallet_dir}/anatQQ.${sub}.nii
    warp_source=${wallet_dir}
    return 0
  fi
  if [[ -f ${new_dir}/anatQQ.${sub}_WARP.nii ]]; then
    warp_field=${new_dir}/anatQQ.${sub}_WARP.nii
    warp_mat=${new_dir}/anatQQ.${sub}.aff12.1D
    anat_master=${new_dir}/anatQQ.${sub}.nii
    warp_source=${new_dir}
    return 0
  fi
  if [[ -f ${legacy_dir}/anat.un.aff.qw_WARP.nii ]]; then
    warp_field=${legacy_dir}/anat.un.aff.qw_WARP.nii
    warp_mat=${legacy_dir}/anat.un.aff.Xat.1D
    anat_master=${legacy_dir}/anat.un.aff.nii
    warp_source=${legacy_dir}
    return 0
  fi
  return 1
}

tmp_dir=$(mktemp -d)
cleanup_tmp() {
  rm -rf "${tmp_dir}"
}
trap cleanup_tmp EXIT

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
    tgt=${atlas_dir}/${label}_probabilistic_inMNI152.nii.gz
    echo "[INFO] Resampling ${label} ROI from N27 -> MNI152 (${tgt})"
    3dresample -overwrite -master "${mni_template}" -input "${src}" -prefix "${tgt}"
    roi_sources[$label]="${tgt}"
  else
    echo "[WARN] Missing ${label} template ${src}; skipping"
    unset "roi_sources[$label]"
  fi
done
if [[ -n ${BLACKANT_MASK_LABELS:-} ]]; then
  read -r -a requested_labels <<<"${BLACKANT_MASK_LABELS}"
  ordered_labels=()
  for raw_label in "${requested_labels[@]}"; do
    upper_label=${raw_label^^}
    case ${upper_label} in
      SN|SNCTRL|PTCTRL)
        ordered_labels+=("${upper_label}")
        ;;
      *)
        echo "[WARN] Ignoring unknown label '${raw_label}' in BLACKANT_MASK_LABELS" >&2
        ;;
    esac
  done
  if [[ ${#ordered_labels[@]} -eq 0 ]]; then
    echo "[WARN] BLACKANT_MASK_LABELS produced no valid labels; falling back to all" >&2
    ordered_labels=(SN SNCTRL PTCTRL)
  else
    echo "[INFO] Restricting ROI build to labels: ${ordered_labels[*]}"
  fi
else
  ordered_labels=(SN SNCTRL PTCTRL)
fi

# SNCTRL subtraction needs an SN mask; add it automatically when missing.
need_sn_for_ctrl=0
for label in "${ordered_labels[@]}"; do
  if [[ ${label^^} == SNCTRL ]]; then
    need_sn_for_ctrl=1
    break
  fi
done
if [[ ${need_sn_for_ctrl} -eq 1 ]]; then
  has_sn=0
  for label in "${ordered_labels[@]}"; do
    if [[ ${label^^} == SN ]]; then
      has_sn=1
      break
    fi
  done
  if [[ ${has_sn} -eq 0 ]]; then
    ordered_labels=(SN "${ordered_labels[@]}")
    echo "[INFO] Added SN label so SNCTRL subtraction can run"
  fi
fi

if [[ -n ${BLACKANT_MASK_SUBJECTS:-} ]]; then
  read -r -a subjects <<<"${BLACKANT_MASK_SUBJECTS}"
  echo "[INFO] Restricting to explicit subjects: ${subjects[*]}"
else
  mapfile -t subjects <"${subjects_file}"
fi
echo "[INFO] Processing ${#subjects[@]} MOTIP subjects"

for entry in "${subjects[@]}"; do
  [[ -z ${entry} ]] && continue
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
  sn_threshold_mask="${tmp_dir}/sn_thresholdmask_${sub}.nii.gz"
  rm -f "${sn_threshold_mask}"
  for label in "${ordered_labels[@]}"; do
    [[ -n ${roi_sources[$label]:-} ]] || continue
    src_roi=${roi_sources[$label]}
    [[ -f ${src_roi} ]] || continue
    tmp_anat="${tmp_dir}/${label}_${sub}_anat.nii.gz"
    out_subject="${out_dir}/${label,,}_groupmask_in_${sub}.nii.gz"
    echo "  -> ${label}: template->anat"
    3dNwarpApply -overwrite -nwarp "INV(${warp_field}) INV(${warp_mat})" \
      -source "${src_roi}" -master "${anat_master}" -ainterp wsinc5 -prefix "${tmp_anat}"
    echo "     ${label}: anat->TSE (${out_subject})"
    3dAllineate -overwrite -base "${tse_vol}" -1Dmatrix_apply "${anat_to_tse}" \
      -cmass -input "${tmp_anat}" -prefix "${out_subject}" -final wsinc5 -float
    if [[ ${label} == "SN" ]]; then
      3dcalc -overwrite -a "${out_subject}" \
        -expr "step(a-${sn_sub_threshold})" -prefix "${sn_threshold_mask}"
    fi
    if [[ ${label} == "SNCTRL" ]]; then
      if [[ ${subtract_sn} -eq 0 ]]; then
        echo "     ${label}: skipping SN overlap subtraction"
      else
        if [[ -s ${sn_threshold_mask} ]]; then
          tmp_ctrl="${tmp_dir}/${label,,}_${sub}_tmp.nii.gz"
          mv "${out_subject}" "${tmp_ctrl}"
          echo "     ${label}: removing SN overlap"
          3dcalc -overwrite -a "${tmp_ctrl}" -b "${sn_threshold_mask}" \
            -expr 'a*(1-step(b))' -prefix "${out_subject}"
        else
          echo "     [WARN] SN mask missing for ${sub}; cannot subtract overlap"
        fi
      fi
    fi
    if [[ ${label} == "SNCTRL" && -n ${ctrl_prob_threshold} ]]; then
      thr_out="${out_dir}/snctrl_thresholdmask_in_${sub}.nii.gz"
      echo "     ${label}: threshold >= ${ctrl_prob_threshold} (${thr_out})"
      3dcalc -overwrite -a "${out_subject}" \
        -expr "step(a-${ctrl_prob_threshold})" -prefix "${thr_out}"
    fi
  done

  stage_brainstem_masks "${sub}"
  record_manifest_entry "${sub}"

done

echo "[INFO] Outputs written to ${out_dir}"
echo "[INFO] Brainstem/DC masks staged under ${brainstem_dir}"
echo "[INFO] Manifest saved to ${manifest_path}"
