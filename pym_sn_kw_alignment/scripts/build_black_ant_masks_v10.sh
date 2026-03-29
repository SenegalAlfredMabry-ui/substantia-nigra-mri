#!/bin/bash -l
#SBATCH -J black_ant_masks_v10
#SBATCH -p cloud
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -c 4
#SBATCH -o sbatch-%j.out
#SBATCH -e sbatch-%j.err
#SBATCH --exclusive
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=FAIL

# build_black_ant_masks_v10.sh
# -----------------------------------------------------------------------------
# PURPOSE
#   Rebuild the MOTIP subject-space neuromelanin probability masks (SN/SNCTRL)
#   and warp every PT-control atlas variant using the refreshed
#   supervillain_wallet warps. V10 keeps the v9 guards, stages the traced
#   nigrosome masks (N1–N5, bilaterally), and prefers the localized antpunk
#   template→TSE warps for the Virgo pilot subjects (falling back to
#   supervillain stacks when ANTPUNK is unavailable). Outputs land under
#   probabilistic_masks_v10 so Irredeemable V10 can ingest the
#   nigrosome-ready directories directly.
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
#   - Historical MOTIP directory mounted at `${MOTIP_ROOT}`.
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
require_wallet_warps=${BLACKANT_REQUIRE_WALLET_WARPS:-1}
antpunk_root=${BLACKANT_ANTPUNK_ROOT:-${snheart_repo}/MRI_data/WARPS_antpunk}
antpunk_manifest=${BLACKANT_ANTPUNK_MANIFEST:-${antpunk_root}/warp_manifest_antpunk.tsv}
prefer_antpunk=${BLACKANT_PREFER_ANTPUNK:-1}
antpunk_required_status=${BLACKANT_ANTPUNK_REQUIRED_STATUS:-ANTPUNK}
roi_offset_json=${BLACKANT_ROI_OFFSET_JSON:-}
apply_roi_offsets=${BLACKANT_APPLY_ROI_OFFSETS:-1}
roi_offset_interp=${BLACKANT_ROI_OFFSET_INTERP:-wsinc5}
roi_offset_sign=${BLACKANT_ROI_OFFSET_SIGN:--1}
# Public GitHub copy: the original pilot-priority roster was replaced with
# example IDs. Override with BLACKANT_PRIORITY_SUBJECTS for a real run.
priority_subjects_default="100 101 102"
priority_subjects_env=${BLACKANT_PRIORITY_SUBJECTS:-}
priority_subjects=(${priority_subjects_env:-$priority_subjects_default})
template_to_anat_chain=""
anat_to_tse_override=""
warp_stack_label=""
declare -A offset_dx=()
declare -A offset_dy=()
declare -A offset_dz=()
declare -A offset_warned=()
declare -A subject_offset_applied=()
declare -A subject_offset_failed=()

if [[ ${prefer_antpunk} -eq 1 ]]; then
  if [[ -n ${antpunk_required_status} ]]; then
    echo "[INFO] Preferring antpunk warps when qc_status=${antpunk_required_status}"
  else
    echo "[INFO] Preferring antpunk warps when available (no qc_status requirement)"
  fi
fi
if [[ -n ${roi_offset_json} && ${apply_roi_offsets} -eq 1 ]]; then
  if [[ ! -f ${roi_offset_json} ]]; then
    echo "[ERROR] BLACKANT_ROI_OFFSET_JSON is set but file is missing: ${roi_offset_json}" >&2
    exit 1
  fi
  if [[ ! ${roi_offset_sign} =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
    echo "[ERROR] BLACKANT_ROI_OFFSET_SIGN must be numeric (got ${roi_offset_sign})" >&2
    exit 1
  fi
  echo "[INFO] SSwarper ROI offsets enabled from ${roi_offset_json} (interp=${roi_offset_interp}, sign=${roi_offset_sign})"
fi

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

out_root=${BLACKANT_OUTPUT_DIR:-${snheart_root}/TSE/probabilistic_masks_v10}
base_mask_dir=${out_root}/base
mkdir -p "${base_mask_dir}"
atlas_dir=${BLACKANT_ATLAS_DIR:-${snheart_root}/ROIs/black_ant_templates_v10}
mkdir -p "${atlas_dir}"
manifest_path=${BLACKANT_MASK_MANIFEST:-${out_root}/black_ant_mask_manifest_v10.tsv}
printf "subject\tsn_voxels\tsnctrl_voxels\tptctrl_voxels\tbrainstem_mask\tbrainstem4v_mask\tdc_mask\tnotes\n" >"${manifest_path}"
pt_variant_manifest=${BLACKANT_PT_VARIANT_MANIFEST:-${out_root}/pt_variant_manifest_v10.tsv}
pt_summary=${BLACKANT_PT_VARIANT_SUMMARY:-${out_root}/pt_variant_summary_v10.tsv}
enable_nigrosomes=${BLACKANT_ENABLE_NIGROSOMES:-1}

subtract_sn=${BLACKANT_SUBTRACT_SN:-1}
if [[ ${subtract_sn} -eq 0 ]]; then
  echo "[INFO] SN subtraction disabled (BLACKANT_SUBTRACT_SN=0)"
fi

sn_sub_threshold=${BLACKANT_SN_SUB_THRESHOLD:-0.15}
echo "[INFO] SN subtraction threshold >= ${sn_sub_threshold} (BLACKANT_SN_SUB_THRESHOLD)"

# By default, keep manual/legacy N1 and N2 masks intact instead of clipping
# them to the SN mask; set BLACKANT_CLIP_N1N2=1 to restore the previous
# behaviour.
clip_n1n2=${BLACKANT_CLIP_N1N2:-0}
if [[ ${clip_n1n2} -eq 1 ]]; then
  echo "[INFO] SN clipping enabled for N1/N2 masks (BLACKANT_CLIP_N1N2=1)"
else
  echo "[INFO] SN clipping disabled for N1/N2 masks"
fi

brainstem_dir=${BLACKANT_BRAINSTEM_DIR:-${snheart_root}/TSE/brainstem_masks_v10}
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

find_template_file() {
  local stem=$1
  local candidates=("${stem}" "${stem}.nii" "${stem}.nii.gz")
  for candidate in "${candidates[@]}"; do
    if [[ -f ${candidate} ]]; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

extract_subject_id() {
  local entry=$1
  if [[ -z ${entry} ]]; then
    echo ""
    return 0
  fi
  if [[ ${entry} == *-* ]]; then
    echo "${entry##*-}"
    return 0
  fi
  echo "${entry}"
}

get_antpunk_qc_status() {
  local sub=$1
  local qc=""
  local status_file=${antpunk_root}/${sub}/qc_status.txt
  if [[ -f ${status_file} ]]; then
    qc=$(<"${status_file}")
  fi
  if [[ -z ${qc} && -f ${antpunk_manifest} ]]; then
    qc=$(awk -F $'\t' -v target="${sub}" 'NR==1{next} $1==target {print $6; exit}' "${antpunk_manifest}" 2>/dev/null)
  fi
  echo "${qc}"
}

mask_voxel_count() {
  local mask_path=$1
  if [[ -z ${mask_path} || ! -f ${mask_path} ]]; then
    echo 0
    return 0
  fi
  local output
  if ! output=$(3dBrickStat -count -non-zero "${mask_path}" 2>/dev/null); then
    echo 0
    return 0
  fi
  awk 'END{print $NF+0}' <<<"${output}"
}

strip_nii_ext() {
  local name
  name=$(basename "$1")
  name=${name%.nii.gz}
  name=${name%.nii}
  echo "${name}"
}

load_subject_offsets() {
  [[ -n ${roi_offset_json} && ${apply_roi_offsets} -eq 1 ]] || return 0
  local line_count=0
  while IFS=$'\t' read -r sub dx dy dz; do
    [[ -z ${sub} ]] && continue
    offset_dx["${sub}"]=${dx}
    offset_dy["${sub}"]=${dy}
    offset_dz["${sub}"]=${dz}
    line_count=$((line_count + 1))
  done < <(python3 - "${roi_offset_json}" <<'PY'
import json
import math
import sys

path = sys.argv[1]
with open(path, encoding="utf-8") as f:
    data = json.load(f)
for sub in sorted(data.keys(), key=lambda x: (len(x), x)):
    entry = data.get(sub) or {}
    delta = entry.get("delta_mm")
    if not isinstance(delta, list) or len(delta) != 3:
        continue
    try:
        vals = [float(v) for v in delta]
    except Exception:
        continue
    if any(math.isnan(v) for v in vals):
        continue
    print(f"{sub}\t{vals[0]:.6f}\t{vals[1]:.6f}\t{vals[2]:.6f}")
PY
)
  echo "[INFO] Loaded ${line_count} per-subject offsets"
}

apply_subject_offset_tse() {
  local sub=$1
  local mask_path=$2
  local label=${3:-mask}
  local interp=${4:-${roi_offset_interp}}
  [[ -n ${roi_offset_json} && ${apply_roi_offsets} -eq 1 ]] || return 0
  if [[ ${warp_stack_label:-} != "sswarper" ]]; then
    return 0
  fi
  if [[ ! -f ${mask_path} ]]; then
    return 0
  fi
  if [[ -z ${offset_dx[${sub}]:-} || -z ${offset_dy[${sub}]:-} || -z ${offset_dz[${sub}]:-} ]]; then
    if [[ -z ${offset_warned[${sub}]:-} ]]; then
      echo "     [WARN] No ROI offset entry for ${sub}; leaving ${label} unshifted"
      offset_warned["${sub}"]=1
    fi
    return 0
  fi
  local apply_dx apply_dy apply_dz
  apply_dx=$(awk -v v="${offset_dx[${sub}]}" -v s="${roi_offset_sign}" 'BEGIN{printf "%.6f", v*s}')
  apply_dy=$(awk -v v="${offset_dy[${sub}]}" -v s="${roi_offset_sign}" 'BEGIN{printf "%.6f", v*s}')
  apply_dz=$(awk -v v="${offset_dz[${sub}]}" -v s="${roi_offset_sign}" 'BEGIN{printf "%.6f", v*s}')
  local sign_tag
  sign_tag=$(tr -c '0-9A-Za-z._-' '_' <<<"${roi_offset_sign}")
  local mat_file="${tmp_dir}/offset_${sub}_${sign_tag}.1D"
  if [[ ! -f ${mat_file} ]]; then
    printf "1 0 0 %s 0 1 0 %s 0 0 1 %s\n" \
      "${apply_dx}" "${apply_dy}" "${apply_dz}" > "${mat_file}"
  fi
  local stem
  stem=$(strip_nii_ext "${mask_path}")
  local shifted="${tmp_dir}/${stem}_offset_${sub}.nii.gz"
  if 3dAllineate -overwrite -base "${mask_path}" -input "${mask_path}" \
      -1Dmatrix_apply "${mat_file}" -prefix "${shifted}" -final "${interp}" -float >/dev/null 2>&1; then
    mv -f "${shifted}" "${mask_path}"
    subject_offset_applied["${sub}"]=1
    echo "     ${label}: applied sswarper offset dx=${apply_dx},dy=${apply_dy},dz=${apply_dz} (raw=${offset_dx[${sub}]},${offset_dy[${sub}]},${offset_dz[${sub}]}; sign=${roi_offset_sign})"
  else
    subject_offset_failed["${sub}"]=1
    rm -f "${shifted}"
    echo "     [WARN] ${label}: failed applying sswarper offset for ${sub}"
  fi
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

warp_nigrosome_masks() {
  local sub=$1
  local sn_clip_mask=${2:-}
  local assigned_mask=${3:-}
  local non_n2_overlap_mask=${4:-}
  local clip_available=0
  if [[ -n ${sn_clip_mask} && -s ${sn_clip_mask} ]]; then
    clip_available=1
  else
    echo "  !! SN threshold mask missing for ${sub}; nigrosomes will not be SN-clipped"
  fi
  declare -A local_counts=()
  declare -A trimmed_counts=()
  for label in "${nigrosome_labels[@]}"; do
    local src=${nigrosome_sources[$label]:-}
    if [[ -z ${src} ]]; then
      local_counts[$label]=0
      trimmed_counts[$label]=0
      continue
    fi
    local tmp_anat="${tmp_dir}/${label}_${sub}_anat.nii.gz"
    local out_subject="${base_mask_dir}/${label}_groupmask_in_${sub}.nii.gz"
    local is_n1_variant=0
    local is_n2_variant=0
    if [[ ${label} == nigrosome1_* ]]; then
      is_n1_variant=1
    elif [[ ${label} == nigrosome2_* ]]; then
      is_n2_variant=1
    fi
    echo "  -> ${label}: template->anat"
    3dNwarpApply -overwrite -nwarp "${template_to_anat_chain}" \
      -source "${src}" -master "${anat_master}" -ainterp wsinc5 -prefix "${tmp_anat}"
    echo "     ${label}: anat->TSE (${out_subject})"
    3dAllineate -overwrite -base "${tse_vol}" -1Dmatrix_apply "${anat_to_tse}" \
      -cmass -input "${tmp_anat}" -prefix "${out_subject}" -final wsinc5 -float
    apply_subject_offset_tse "${sub}" "${out_subject}" "${label}" "wsinc5"
    local skip_clip=0
    if [[ ${clip_available} -eq 1 ]]; then
      if [[ ${clip_n1n2} -eq 0 && ( ${is_n1_variant} -eq 1 || ${is_n2_variant} -eq 1 ) ]]; then
        skip_clip=1
      fi
    fi
    if [[ ${clip_available} -eq 1 && ${skip_clip} -eq 0 ]]; then
      local clipped="${tmp_dir}/${label}_${sub}_snclipped.nii.gz"
      echo "     ${label}: clipping to SN mask"
      3dcalc -overwrite -a "${out_subject}" -b "${sn_clip_mask}" \
        -expr 'a*step(b)' -prefix "${clipped}"
      mv -f "${clipped}" "${out_subject}"
    fi
    local before_exclusive
    before_exclusive=$(mask_voxel_count "${out_subject}")
    local after_exclusive=${before_exclusive}
    local trimmed=0
    local overlap_source=""
    if [[ ${is_n2_variant} -eq 1 ]]; then
      overlap_source=${non_n2_overlap_mask}
    else
      overlap_source=${assigned_mask}
    fi
    if [[ -n ${overlap_source} && -f ${overlap_source} ]]; then
      local tmp_exclusive="${tmp_dir}/${label}_${sub}_exclusive.nii.gz"
      echo "     ${label}: removing overlap with earlier nigrosomes"
      3dcalc -overwrite -a "${out_subject}" -b "${overlap_source}" \
        -expr 'a*(1-step(b))' -prefix "${tmp_exclusive}"
      mv -f "${tmp_exclusive}" "${out_subject}"
      after_exclusive=$(mask_voxel_count "${out_subject}")
      trimmed=$((before_exclusive - after_exclusive))
      if [[ ${trimmed} -lt 0 ]]; then
        trimmed=0
      fi
      if [[ ${trimmed} -gt 0 ]]; then
        echo "     ${label}: trimmed ${trimmed} overlapping voxels"
      fi
    fi
    if [[ -n ${assigned_mask} ]]; then
      if [[ -f ${assigned_mask} ]]; then
        local tmp_assigned="${tmp_dir}/assigned_${sub}.nii.gz"
        3dcalc -overwrite -a "${assigned_mask}" -b "${out_subject}" \
          -expr 'step(a)+step(b)' -prefix "${tmp_assigned}" >/dev/null 2>&1 || true
        3dcalc -overwrite -a "${tmp_assigned}" -expr 'step(a)' -prefix "${assigned_mask}" >/dev/null 2>&1 || true
      else
        3dcalc -overwrite -a "${out_subject}" -expr 'step(a)' -prefix "${assigned_mask}" >/dev/null 2>&1 || true
      fi
    fi
    if [[ ${is_n2_variant} -eq 0 && -n ${non_n2_overlap_mask} ]]; then
      if [[ -f ${non_n2_overlap_mask} ]]; then
        local tmp_non_n2="${tmp_dir}/assigned_non_n2_${sub}.nii.gz"
        3dcalc -overwrite -a "${non_n2_overlap_mask}" -b "${out_subject}" \
          -expr 'step(a)+step(b)' -prefix "${tmp_non_n2}" >/dev/null 2>&1 || true
        3dcalc -overwrite -a "${tmp_non_n2}" -expr 'step(a)' -prefix "${non_n2_overlap_mask}" >/dev/null 2>&1 || true
      else
        3dcalc -overwrite -a "${out_subject}" -expr 'step(a)' -prefix "${non_n2_overlap_mask}" >/dev/null 2>&1 || true
      fi
    fi
    local_counts[$label]=${after_exclusive}
    trimmed_counts[$label]=${trimmed}
  done
  {
    printf "%s" "${sub}"
    for id in "${nigro_ids[@]}"; do
      local idx=${id#N}
      for hemi in "${nigro_hemis[@]}"; do
        local label="nigrosome${idx}_${hemi,,}"
        printf "\t%s\t%s" "${local_counts[$label]:-0}" "${trimmed_counts[$label]:-0}"
      done
      local left_label="nigrosome${idx}_l"
      local right_label="nigrosome${idx}_r"
      local left_count=${local_counts[$left_label]:-0}
      local right_count=${local_counts[$right_label]:-0}
      local left_trim=${trimmed_counts[$left_label]:-0}
      local right_trim=${trimmed_counts[$right_label]:-0}
      printf "\t%s\t%s" $((left_count + right_count)) $((left_trim + right_trim))
    done
    for variant in "${enabled_n2_variants[@]}"; do
      for hemi in "${nigro_hemis[@]}"; do
        local label="nigrosome2_${variant}_${hemi,,}"
        printf "\t%s\t%s" "${local_counts[$label]:-0}" "${trimmed_counts[$label]:-0}"
      done
      local left_label="nigrosome2_${variant}_l"
      local right_label="nigrosome2_${variant}_r"
      local left_count=${local_counts[$left_label]:-0}
      local right_count=${local_counts[$right_label]:-0}
      local left_trim=${trimmed_counts[$left_label]:-0}
      local right_trim=${trimmed_counts[$right_label]:-0}
      printf "\t%s\t%s" $((left_count + right_count)) $((left_trim + right_trim))
    done
    printf "\n"
  } >>"${nigrosome_manifest}"
}

link_base_masks_to_variants() {
  local sub=$1
  local sn_path="${base_mask_dir}/sn_groupmask_in_${sub}.nii.gz"
  local snctrl_path="${base_mask_dir}/snctrl_groupmask_in_${sub}.nii.gz"
  for variant in "${!pt_variant_dirs[@]}"; do
    local variant_dir=${pt_variant_dirs[$variant]}
    if [[ -f ${sn_path} ]]; then
      ln -sf "${sn_path}" "${variant_dir}/sn_groupmask_in_${sub}.nii.gz"
    fi
    if [[ -f ${snctrl_path} ]]; then
      ln -sf "${snctrl_path}" "${variant_dir}/snctrl_groupmask_in_${sub}.nii.gz"
    fi
  done
}

record_pt_variant_manifest() {
  local sub=$1
  local variant=$2
  local mask_path=$3
  local voxels=$4
  local notes=$5
  printf "%s\t%s\t%s\t%s\t%s\n" \
    "${sub}" "${variant}" "${mask_path}" "${voxels}" "${notes}" >>"${pt_variant_manifest}"
}

record_manifest_entry() {
  local sub=$1
  local sn_mask snctrl_mask ptctrl_mask brainstem_mask brainstem4v_mask dc_mask
  sn_mask="$(first_existing_mask "${base_mask_dir}/sn_groupmask_in_${sub}" || true)"
  snctrl_mask="$(first_existing_mask "${base_mask_dir}/snctrl_groupmask_in_${sub}" || true)"
  ptctrl_mask="$(first_existing_mask "${base_mask_dir}/ptctrl_groupmask_in_${sub}" || true)"
  brainstem_mask="$(first_existing_mask "${brainstem_dir}/brainstem_mask_inTSE_${sub}" || true)"
  brainstem4v_mask="$(first_existing_mask "${brainstem_dir}/brainstem_4v_mask_inTSE_${sub}" || true)"
  dc_mask="$(first_existing_mask "${brainstem_dir}/DC_mask_inTSE_${sub}" || true)"

  local notes=()
  [[ -z ${sn_mask} ]] && notes+=(sn_missing)
  [[ -z ${snctrl_mask} ]] && notes+=(snctrl_missing)
  [[ -z ${ptctrl_mask} ]] && notes+=(ptctrl_missing)
  [[ -z ${brainstem_mask} ]] && notes+=(brainstem_missing)
  [[ -z ${brainstem4v_mask} ]] && notes+=(brainstem4v_missing)
  [[ -z ${dc_mask} ]] && notes+=(dc_missing)
  if [[ -n ${warp_stack_label:-} ]]; then
    notes+=(warp_${warp_stack_label})
  fi
  if [[ -n ${subject_offset_applied[${sub}]:-} ]]; then
    notes+=(offset_applied)
  fi
  if [[ -n ${subject_offset_failed[${sub}]:-} ]]; then
    notes+=(offset_failed)
  fi

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

build_pt_variant_mask() {
  local variant=$1
  local sub=$2
  [[ -n ${pt_variant_sources[$variant]:-} ]] || return 0
  local src_roi=${pt_variant_sources[$variant]}
  local tmp_anat="${tmp_dir}/ptctrl_${variant}_${sub}_anat.nii.gz"
  local variant_dir=${pt_variant_dirs[$variant]}
  local out_subject="${variant_dir}/ptctrl_groupmask_in_${sub}.nii.gz"
  echo "  -> PT ${variant}: template->anat"
  3dNwarpApply -overwrite -nwarp "${template_to_anat_chain}" \
    -source "${src_roi}" -master "${anat_master}" -ainterp wsinc5 -prefix "${tmp_anat}"
  echo "     PT ${variant}: anat->TSE (${out_subject})"
  3dAllineate -overwrite -base "${tse_vol}" -1Dmatrix_apply "${anat_to_tse}" \
    -cmass -input "${tmp_anat}" -prefix "${out_subject}" -final wsinc5 -float
  local voxels
  voxels=$(mask_voxel_count "${out_subject}")
  local note="ok"
  if [[ ${voxels} -eq 0 ]]; then
    note="empty"
  fi
  record_pt_variant_manifest "${sub}" "${variant}" "${out_subject}" "${voxels}" "${note}"
  if [[ -n ${primary_pt_variant} && ${variant} == ${primary_pt_variant} ]]; then
    ln -sf "${out_subject}" "${base_mask_dir}/ptctrl_groupmask_in_${sub}.nii.gz"
  fi
}

build_pt_superior_mask() {
  local sub=$1
  if [[ ${enable_pt_superior} -ne 1 ]]; then
    return 0
  fi
  local variant_dir=${pt_variant_dirs[SUPERIOR]}
  local sn_mask="${base_mask_dir}/sn_groupmask_in_${sub}.nii.gz"
  local brainstem_mask
  brainstem_mask=$(first_existing_mask "${brainstem_dir}/brainstem_mask_inTSE_${sub}")
  local out_subject="${variant_dir}/ptctrl_groupmask_in_${sub}.nii.gz"
  local note="missing"
  local voxels=0
  if [[ -z ${brainstem_mask} ]]; then
    note="brainstem_missing"
  elif [[ ! -f ${sn_mask} ]]; then
    note="sn_missing"
  else
    local tmp_superior="${tmp_dir}/pt_superior_${sub}.nii.gz"
    3dcalc -overwrite -a "${brainstem_mask}" -b "${sn_mask}" \
      -expr 'ispositive(a)*step(b)' -prefix "${tmp_superior}" >/dev/null 2>&1 || true
    if [[ -s ${tmp_superior} ]]; then
      cp -f "${tmp_superior}" "${out_subject}"
      voxels=$(mask_voxel_count "${out_subject}")
      if [[ ${voxels} -gt 0 ]]; then
        note="ok"
      else
        note="empty"
      fi
    else
      note="empty"
      rm -f "${tmp_superior}"
      rm -f "${out_subject}"
    fi
  fi
  record_pt_variant_manifest "${sub}" "SUPERIOR" "${out_subject}" "${voxels}" "${note}"
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
    "${snheart_root}/TSE/zeropad_tse_${sub}.nii"
    "${snheart_root}/TSE/zeropad_tse_${sub}.nii.gz"
    "${snheart_root}/TSE/${sub}_TSE_zeropad.nii"
    "${snheart_root}/TSE/${sub}_TSE_zeropad.nii.gz"
    "${snheart_root}/TSE/${sub}_TSE.nii"
    "${snheart_root}/TSE/${sub}_TSE.nii.gz"
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
  template_to_anat_chain=""
  anat_to_tse_override=""
  warp_stack_label=""
  local legacy_dir=${motip_root}/TSE/${sub}_warptoN27
  local new_dir=${motip_root}/WARPS/${sub}
  local wallet_dir=${wallet_warp_root}/${sub}
  warp_source=""
  local sswarper_found=0
  # Prefer the refreshed supervillain cache before older warp stacks.
  if [[ -f ${wallet_dir}/anatQQ.${sub}_WARP.nii ]]; then
    warp_field=${wallet_dir}/anatQQ.${sub}_WARP.nii
    warp_mat=${wallet_dir}/anatQQ.${sub}.aff12.1D
    anat_master=${wallet_dir}/anatQQ.${sub}.nii
    warp_source=${wallet_dir}
    sswarper_found=1
  elif [[ ${require_wallet_warps} -eq 0 && -f ${new_dir}/anatQQ.${sub}_WARP.nii ]]; then
    warp_field=${new_dir}/anatQQ.${sub}_WARP.nii
    warp_mat=${new_dir}/anatQQ.${sub}.aff12.1D
    anat_master=${new_dir}/anatQQ.${sub}.nii
    warp_source=${new_dir}
    sswarper_found=1
  elif [[ ${require_wallet_warps} -eq 0 && -f ${legacy_dir}/anat.un.aff.qw_WARP.nii ]]; then
    warp_field=${legacy_dir}/anat.un.aff.qw_WARP.nii
    warp_mat=${legacy_dir}/anat.un.aff.Xat.1D
    anat_master=${legacy_dir}/anat.un.aff.nii
    warp_source=${legacy_dir}
    sswarper_found=1
  fi
  if [[ ${sswarper_found} -eq 0 ]]; then
    return 1
  fi

  if [[ ${prefer_antpunk} -eq 1 ]]; then
    local antpunk_dir=${antpunk_root}/${sub}
    local components_file=${antpunk_dir}/template_to_tse_components.txt
    if [[ -f ${components_file} ]]; then
      local qc_status
      qc_status=$(get_antpunk_qc_status "${sub}")
      local allow_antpunk=1
      if [[ -n ${antpunk_required_status} && ${qc_status} != ${antpunk_required_status} ]]; then
        allow_antpunk=0
      fi
      if [[ ${allow_antpunk} -eq 1 ]]; then
        read -r template_to_anat_override anat_to_tse_file <"${components_file}"
        if [[ -n ${template_to_anat_override} ]]; then
          template_to_anat_chain="${template_to_anat_override}"
          anat_to_tse_override=${anat_to_tse_file:-}
          warp_source=${antpunk_dir}
          warp_stack_label="antpunk"
          if [[ -n ${qc_status} ]]; then
            echo "  -> Using antpunk warp (qc_status=${qc_status})"
          else
            echo "  -> Using antpunk warp (qc_status unknown)"
          fi
          return 0
        fi
      else
        echo "  -> Antpunk warp present but qc_status=${qc_status}; requiring ${antpunk_required_status}" >&2
      fi
    fi
  fi

  template_to_anat_chain="INV(${warp_field}) INV(${warp_mat})"
  warp_stack_label="sswarper"
  return 0
}

tmp_dir=$(mktemp -d)
cleanup_tmp() {
  rm -rf "${tmp_dir}"
}
trap cleanup_tmp EXIT
load_subject_offsets

declare -A roi_sources=(
  [SN]="${snheart_sn_prob}"
  [SNCTRL]="${mni_sn_ctrl}"
)

nigro_ids=(N1 N2)
nigro_hemis=(L R)
n2_extra_variants=(test approach)
enabled_n2_variants=()
nigrosome_labels=()
declare -A nigrosome_sources=()
nigrosome_manifest=${BLACKANT_NIGROSOME_MANIFEST:-${out_root}/nigrosome_manifest_v10.tsv}
has_nigrosome_templates=0
if [[ ${enable_nigrosomes} -eq 1 ]]; then
  declare -A n2_variant_prefix=(
    [test]="Test_SN_N2"
    [approach]="Approach_SN_N2"
  )
  declare -A n2_variant_enabled=()
  nigrosome_generator=${BLACKANT_NIGROSOME_GENERATOR:-${snheart_repo}/Scripts/Agent_Anesthesia.py}
  default_nigro_python=${snheart_repo}/venv/bin/python3
  if [[ -n ${BLACKANT_NIGROSOME_PYTHON:-} ]]; then
    nigrosome_python=${BLACKANT_NIGROSOME_PYTHON}
  elif [[ -x ${default_nigro_python} ]]; then
    nigrosome_python=${default_nigro_python}
  else
    nigrosome_python=$(command -v python3)
  fi
  if [[ -z ${nigrosome_python} || ! -x ${nigrosome_python} ]]; then
    echo "[ERROR] Unable to locate a Python interpreter for nigrosome generation; set BLACKANT_NIGROSOME_PYTHON" >&2
    exit 1
  fi
  auto_generate_nigrosomes=${BLACKANT_GENERATE_NIGROSOMES:-0}
  nigrosome_template_dir_default=${snheart_root}/ROIs/nigrosome_templates_auto
  nigrosome_template_dir=${BLACKANT_NIGROSOME_TEMPLATE_DIR:-${nigrosome_template_dir_default}}
  nigrosome_metadata_path=${BLACKANT_NIGROSOME_METADATA:-${nigrosome_template_dir}/summary.json}

  if [[ -n ${BLACKANT_NIGROSOME_ROOT:-} ]]; then
    nigrosome_root=${BLACKANT_NIGROSOME_ROOT}
  else
    if [[ ${auto_generate_nigrosomes} -eq 1 ]]; then
      if [[ ! -x ${nigrosome_generator} ]]; then
        echo "[ERROR] Nigrosome generator ${nigrosome_generator} not executable; set BLACKANT_GENERATE_NIGROSOMES=0 or provide BLACKANT_NIGROSOME_ROOT" >&2
        exit 1
      fi
      echo "[INFO] Generating canonical nigrosome templates via ${nigrosome_generator}"
      mkdir -p "${nigrosome_template_dir}"
      "${nigrosome_python}" "${nigrosome_generator}" \
        --sn-prob "${snheart_sn_prob}" \
        --output-dir "${nigrosome_template_dir}" \
        --metadata "${nigrosome_metadata_path}"
      nigrosome_root=${nigrosome_template_dir}
      echo "[INFO] Nigrosome templates written to ${nigrosome_root} (metadata: ${nigrosome_metadata_path})"
    else
      nigrosome_root=${snheart_root}/ROIs
    fi
  fi

  echo "[INFO] Using nigrosome templates under ${nigrosome_root}"
  for id in "${nigro_ids[@]}"; do
    for hemi in "${nigro_hemis[@]}"; do
      label="nigrosome${id#N}_${hemi,,}"
      base_path="${nigrosome_root}/SN_${id}_${hemi}"
      if src_path=$(find_template_file "${base_path}"); then
        nigrosome_sources[${label}]="${src_path}"
        nigrosome_labels+=("${label}")
      else
        echo "[WARN] Missing nigrosome template ${base_path}{.nii/.nii.gz}; ${label} will be skipped" >&2
      fi
      if [[ ${id} == N2 ]]; then
        for variant in "${n2_extra_variants[@]}"; do
          prefix=${n2_variant_prefix[${variant}]}
          [[ -z ${prefix} ]] && continue
          extra_label="nigrosome${id#N}_${variant}_${hemi,,}"
          extra_base="${nigrosome_root}/${prefix}_${hemi}"
          if src_variant=$(find_template_file "${extra_base}"); then
            nigrosome_sources[${extra_label}]="${src_variant}"
            nigrosome_labels+=("${extra_label}")
            n2_variant_enabled[${variant}]=1
          else
            echo "[WARN] Missing ${variant^^} N2 template ${extra_base}{.nii/.nii.gz}; ${extra_label} will be skipped" >&2
          fi
        done
      fi
    done
  done
  for variant in "${n2_extra_variants[@]}"; do
    if [[ ${n2_variant_enabled[${variant}]:-0} -eq 1 ]]; then
      enabled_n2_variants+=("${variant}")
    fi
  done
  has_nigrosome_templates=${#nigrosome_labels[@]}
  if [[ ${has_nigrosome_templates} -eq 0 ]]; then
    echo "[WARN] No nigrosome templates found under ${nigrosome_root}; manifest will remain header-only" >&2
  else
    echo "[INFO] Enabling nigrosome warps for ${has_nigrosome_templates} templates"
  fi
else
  echo "[INFO] Nigrosome warps disabled (BLACKANT_ENABLE_NIGROSOMES=0)"
fi
{
  printf "subject"
  for id in "${nigro_ids[@]}"; do
    lower=${id,,}
    for hemi in "${nigro_hemis[@]}"; do
      printf "\t%s_%s_voxels\t%s_%s_trimmed" "${lower}" "${hemi,,}" "${lower}" "${hemi,,}"
    done
    printf "\t%s_total_voxels\t%s_total_trimmed" "${lower}" "${lower}"
  done
  for variant in "${enabled_n2_variants[@]}"; do
    for hemi in "${nigro_hemis[@]}"; do
      printf "\t%s_%s_%s_voxels\t%s_%s_%s_trimmed" \
        "n2" "${variant}" "${hemi,,}" "n2" "${variant}" "${hemi,,}"
    done
    printf "\t%s_%s_total_voxels\t%s_%s_total_trimmed" "n2" "${variant}" "n2" "${variant}"
  done
  printf "\n"
} >"${nigrosome_manifest}"

declare -A pt_variant_candidates=(
  [OGRADY]="${snheart_root}/ROIs/LC_control_OGRADY_inMNI152.nii"
  [MOTIP]="${snheart_root}/ROIs/PT_control_MOTIP_inMNI152.nii.gz"
  [MOTIP_ATLASWARP]="${snheart_root}/ROIs/PT_control_MOTIP_inMNI152_atlaswarp.nii.gz"
  [KW]="${snheart_root}/ROIs/PT_control_KW_inMNI152.nii.gz"
  [SAM]="${snheart_root}/ROIs/MNI152_CONTROL_PT_SAM.nii"
)
if [[ -n ${BLACKANT_PT_VARIANTS:-} ]]; then
  read -r -a requested_pt_variants <<<"${BLACKANT_PT_VARIANTS}"
else
  requested_pt_variants=(KW)
fi
declare -A pt_variant_sources=()
pt_variant_labels=()
for raw_variant in "${requested_pt_variants[@]}"; do
  upper_variant=${raw_variant^^}
  if [[ -n ${pt_variant_candidates[${upper_variant}]:-} ]]; then
    pt_variant_sources[${upper_variant}]=${pt_variant_candidates[${upper_variant}]}
    pt_variant_labels+=("${upper_variant}")
  else
    echo "[WARN] Unknown PT variant ${raw_variant}; skipping" >&2
  fi
done
if [[ ${#pt_variant_labels[@]} -eq 0 ]]; then
  echo "[WARN] No valid PT variants requested; defaulting to KW" >&2
  pt_variant_labels=(KW)
  pt_variant_sources[KW]=${pt_variant_candidates[KW]}
fi
IFS=$'\n' pt_variant_labels=($(sort -u <<<"${pt_variant_labels[*]}") )
unset IFS
declare -A pt_variant_dirs=()
for variant in "${pt_variant_labels[@]}"; do
  dir="${out_root}/PT_${variant}"
  mkdir -p "${dir}"
  pt_variant_dirs[${variant}]="${dir}"
done
printf "subject\tvariant\tmask_path\tvoxels\tnotes\n" >"${pt_variant_manifest}"
primary_pt_variant=${BLACKANT_PRIMARY_PT_VARIANT:-KW}
if [[ -n ${primary_pt_variant} && -z ${pt_variant_sources[$primary_pt_variant]:-} ]]; then
  echo "[WARN] Requested primary PT variant ${primary_pt_variant} missing; base directory will not include PT masks" >&2
  primary_pt_variant=""
fi
enable_pt_superior=${BLACKANT_ENABLE_PT_SUPERIOR:-1}
if [[ ${enable_pt_superior} -eq 1 ]]; then
  pt_variant_dirs[SUPERIOR]="${out_root}/PT_SUPERIOR"
  mkdir -p "${pt_variant_dirs[SUPERIOR]}"
fi

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
for label in "${nigrosome_labels[@]}"; do
  src=${nigrosome_sources[$label]:-}
  [[ -z ${src} ]] && continue
  tgt=${atlas_dir}/${label}_probabilistic_inMNI152.nii.gz
  echo "[INFO] Resampling ${label} -> ${tgt}"
  3dresample -overwrite -master "${mni_template}" -input "${src}" -prefix "${tgt}"
  nigrosome_sources[$label]="${tgt}"
done
 for variant in "${pt_variant_labels[@]}"; do
  src=${pt_variant_sources[$variant]}
  if [[ -f ${src} ]]; then
    tgt=${atlas_dir}/PTCTRL_${variant}_probabilistic_inMNI152.nii.gz
    echo "[INFO] Resampling PT variant ${variant} (${tgt})"
    3dresample -overwrite -master "${mni_template}" -input "${src}" -prefix "${tgt}"
    pt_variant_sources[$variant]="${tgt}"
  else
    echo "[WARN] Missing PT variant ${variant} source ${src}; dropping"
    unset "pt_variant_sources[$variant]"
  fi
done
filtered_pt_variants=()
for variant in "${pt_variant_labels[@]}"; do
  if [[ -n ${pt_variant_sources[$variant]:-} ]]; then
    filtered_pt_variants+=("${variant}")
  fi
done
pt_variant_labels=("${filtered_pt_variants[@]}")
if [[ ${#pt_variant_labels[@]} -eq 0 ]]; then
  echo "[WARN] No PT variants survived template checks; PT mask builds will be skipped" >&2
fi
if [[ -n ${BLACKANT_MASK_LABELS:-} ]]; then
  read -r -a requested_labels <<<"${BLACKANT_MASK_LABELS}"
  ordered_labels=()
  for raw_label in "${requested_labels[@]}"; do
    upper_label=${raw_label^^}
    case ${upper_label} in
      SN|SNCTRL)
        ordered_labels+=("${upper_label}")
        ;;
      *)
        echo "[WARN] Ignoring unknown label '${raw_label}' in BLACKANT_MASK_LABELS" >&2
        ;;
    esac
  done
  if [[ ${#ordered_labels[@]} -eq 0 ]]; then
    echo "[WARN] BLACKANT_MASK_LABELS produced no valid labels; falling back to all" >&2
    ordered_labels=(SN SNCTRL)
  else
    echo "[INFO] Restricting ROI build to labels: ${ordered_labels[*]}"
  fi
else
  ordered_labels=(SN SNCTRL)
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

if [[ ${#priority_subjects[@]} -gt 0 ]]; then
  declare -a reordered_subjects=()
  declare -A subject_used=()
  for p in "${priority_subjects[@]}"; do
    [[ -z ${p} ]] && continue
    for idx in "${!subjects[@]}"; do
      [[ -n ${subject_used[$idx]:-} ]] && continue
      entry=${subjects[$idx]}
      [[ -z ${entry} ]] && continue
      sub=$(extract_subject_id "${entry}")
      if [[ ${sub} == ${p} ]]; then
        reordered_subjects+=("${entry}")
        subject_used[$idx]=1
        break
      fi
    done
  done
  for idx in "${!subjects[@]}"; do
    if [[ -z ${subject_used[$idx]:-} ]]; then
      reordered_subjects+=("${subjects[$idx]}")
    fi
  done
  subjects=("${reordered_subjects[@]}")
fi
echo "[INFO] Processing ${#subjects[@]} subjects"

for entry in "${subjects[@]}"; do
  [[ -z ${entry} ]] && continue
  sub=$(extract_subject_id "${entry}")
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
  if [[ -z ${template_to_anat_chain:-} ]]; then
    echo "  !! No template→anat warp chain for ${sub}; skipping"
    continue
  fi
  if [[ -n ${anat_to_tse_override:-} ]]; then
    anat_to_tse=${anat_to_tse_override}
  else
    anat_to_tse=${motip_root}/TSE/${sub}_anat_toTSE_mat.aff12.1D
    if [[ ! -f ${anat_to_tse} ]]; then
      anat_to_tse=${motip_root}/TSE/HUMAN_EBR-${sub}_anat_toTSE_mat.aff12.1D
    fi
    if [[ ! -f ${anat_to_tse} ]]; then
      anat_to_tse=${snheart_root}/align/${sub}_MPRAGEANATtoTSE_${sub}_mat.aff12.1D
    fi
  fi
  if [[ ! -f ${anat_to_tse} ]]; then
    echo "  !! Missing ${anat_to_tse}; skipping"
    continue
  fi
  if [[ -n ${warp_source:-} ]]; then
    if [[ ${warp_stack_label} == "antpunk" ]]; then
      echo "  -> Using antpunk stack from ${warp_source}"
    else
      echo "  -> Using SSwarper stack from ${warp_source}"
    fi
  fi

  sn_mask_path="${base_mask_dir}/sn_groupmask_in_${sub}.nii.gz"
  sn_threshold_mask="${tmp_dir}/sn_thresholdmask_${sub}.nii.gz"
  rm -f "${sn_threshold_mask}"
  for label in "${ordered_labels[@]}"; do
    [[ -n ${roi_sources[$label]:-} ]] || continue
    src_roi=${roi_sources[$label]}
    [[ -f ${src_roi} ]] || continue
    tmp_anat="${tmp_dir}/${label}_${sub}_anat.nii.gz"
    out_subject="${base_mask_dir}/${label,,}_groupmask_in_${sub}.nii.gz"
    echo "  -> ${label}: template->anat"
    3dNwarpApply -overwrite -nwarp "${template_to_anat_chain}" \
      -source "${src_roi}" -master "${anat_master}" -ainterp wsinc5 -prefix "${tmp_anat}"
    echo "     ${label}: anat->TSE (${out_subject})"
    3dAllineate -overwrite -base "${tse_vol}" -1Dmatrix_apply "${anat_to_tse}" \
      -cmass -input "${tmp_anat}" -prefix "${out_subject}" -final wsinc5 -float
    if [[ ${label} == "SN" ]]; then
      apply_subject_offset_tse "${sub}" "${out_subject}" "${label}" "wsinc5"
    fi
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
      thr_out="${base_mask_dir}/snctrl_thresholdmask_in_${sub}.nii.gz"
      echo "     ${label}: threshold >= ${ctrl_prob_threshold} (${thr_out})"
      3dcalc -overwrite -a "${out_subject}" \
        -expr "step(a-${ctrl_prob_threshold})" -prefix "${thr_out}"
    fi
  done

  stage_brainstem_masks "${sub}"
  link_base_masks_to_variants "${sub}"
  for variant in "${!pt_variant_sources[@]}"; do
    build_pt_variant_mask "${variant}" "${sub}"
  done
  build_pt_superior_mask "${sub}"
  if [[ ${has_nigrosome_templates} -gt 0 ]]; then
    subject_nigro_assigned="${tmp_dir}/nigrosome_${sub}_assigned.nii.gz"
    subject_nigro_assigned_non_n2="${tmp_dir}/nigrosome_${sub}_assigned_non_n2.nii.gz"
    rm -f "${subject_nigro_assigned}" "${subject_nigro_assigned_non_n2}"
    warp_nigrosome_masks "${sub}" "${sn_threshold_mask}" "${subject_nigro_assigned}" "${subject_nigro_assigned_non_n2}"
    rm -f "${subject_nigro_assigned}" "${subject_nigro_assigned_non_n2}"
  fi
  record_manifest_entry "${sub}"

done

echo "[INFO] Outputs written to ${base_mask_dir}"
echo "[INFO] Brainstem/DC masks staged under ${brainstem_dir}"
echo "[INFO] Manifest saved to ${manifest_path}"
echo "[INFO] PT variant manifest saved to ${pt_variant_manifest}"
if [[ ${has_nigrosome_templates} -gt 0 ]]; then
  echo "[INFO] Nigrosome manifest (with trimmed voxel counts) saved to ${nigrosome_manifest}"
else
  echo "[INFO] Nigrosome manifest header saved to ${nigrosome_manifest} (no templates warped)"
fi

python3 <<PY
import csv
import statistics as stats
from collections import defaultdict

manifest_path = "${pt_variant_manifest}"
summary_path = "${pt_summary}"
volumes = defaultdict(list)
with open(manifest_path, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        try:
            vol = float(row['voxels'])
        except ValueError:
            continue
        volumes[row['variant']].append(vol)

with open(summary_path, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['variant', 'n', 'min', 'median', 'max', 'mean'])
    for variant in sorted(volumes):
        vals = sorted(volumes[variant])
        if not vals:
            writer.writerow([variant, 0, '', '', '', ''])
            continue
        writer.writerow([
            variant,
            len(vals),
            vals[0],
            stats.median(vals),
            vals[-1],
            sum(vals)/len(vals),
        ])
PY
echo "[INFO] PT variant summary written to ${pt_summary}"
