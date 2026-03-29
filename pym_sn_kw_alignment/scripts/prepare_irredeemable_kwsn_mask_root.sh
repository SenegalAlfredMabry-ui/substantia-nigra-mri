#!/usr/bin/env bash
set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
motip_root=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data

source_root=${KWBACKUP_SOURCE_ROOT:-${repo_root}/MRI_data/TSE/probabilistic_masks_v10}
out_root=${KWBACKUP_OUTPUT_ROOT:-${repo_root}/MRI_data/TSE/probabilistic_masks_kwsn_v9backup}
pt_variants_raw=${KWBACKUP_PT_VARIANTS:-KW}
subjects_override=${KWBACKUP_SUBJECTS:-}

usage() {
  cat <<'EOF'
Usage: prepare_irredeemable_kwsn_mask_root.sh

Environment:
  KWBACKUP_SOURCE_ROOT   Existing black-ant mask root providing SNCTRL/PT masks.
                         Default: MRI_data/TSE/probabilistic_masks_v10
  KWBACKUP_OUTPUT_ROOT   Output root for the KW-backed Irredeemable masks.
                         Default: MRI_data/TSE/probabilistic_masks_kwsn_v9backup
  KWBACKUP_PT_VARIANTS   PT variants to stage (space-delimited). Default: KW
  KWBACKUP_SUBJECTS      Optional explicit subject list. Default: discover from
                         source_root/base/sn_groupmask_in_*.nii*
EOF
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

[[ -d ${source_root} ]] || { echo "[ERROR] Source root not found: ${source_root}" >&2; exit 1; }
[[ -d ${motip_root}/TSE/SN ]] || { echo "[ERROR] MOTIP SN dir missing: ${motip_root}/TSE/SN" >&2; exit 1; }

read -r -a pt_variants <<<"${pt_variants_raw}"
[[ ${#pt_variants[@]} -gt 0 ]] || pt_variants=(KW)

resolve_kw_mask() {
  local sub=$1
  local candidates=(
    "${motip_root}/TSE/SN/SN_ROI_KW_${sub}_ZP.nii"
    "${motip_root}/TSE/SN/SN_ROI_KW_${sub}_ZP.nii.gz"
    "${motip_root}/TSE/SN/SN_ROI_KW_${sub}.nii"
    "${motip_root}/TSE/SN/SN_ROI_KW_${sub}.nii.gz"
  )
  local p
  for p in "${candidates[@]}"; do
    if [[ -f ${p} ]]; then
      printf '%s\n' "${p}"
      return 0
    fi
  done
  return 1
}

locate_source_mask() {
  local dir=$1
  local stem=$2
  local candidates=("${dir}/${stem}.nii.gz" "${dir}/${stem}.nii")
  local p
  for p in "${candidates[@]}"; do
    if [[ -f ${p} ]]; then
      printf '%s\n' "${p}"
      return 0
    fi
  done
  return 1
}

subjects=()
if [[ -n ${subjects_override} ]]; then
  read -r -a subjects <<<"${subjects_override}"
else
  while IFS= read -r path; do
    name=$(basename "${path}")
    sub=${name#sn_groupmask_in_}
    sub=${sub%.nii.gz}
    sub=${sub%.nii}
    subjects+=("${sub}")
  done < <(find "${source_root}/base" -maxdepth 1 -type f \( -name 'sn_groupmask_in_*.nii' -o -name 'sn_groupmask_in_*.nii.gz' \) | sort)
fi

[[ ${#subjects[@]} -gt 0 ]] || { echo "[ERROR] No subjects resolved" >&2; exit 1; }

mkdir -p "${out_root}/base"

manifest=${out_root}/kwsn_backup_manifest.tsv
printf "subject\tkw_sn\tbase_snctrl\tbase_ptctrl\tstatus\tnotes\n" > "${manifest}"

for sub in "${subjects[@]}"; do
  kw_path=$(resolve_kw_mask "${sub}" || true)
  base_dir="${out_root}/base"
  snctrl_src=$(locate_source_mask "${source_root}/base" "snctrl_groupmask_in_${sub}" || true)
  ptctrl_src=$(locate_source_mask "${source_root}/base" "ptctrl_groupmask_in_${sub}" || true)
  notes=()

  rm -f "${base_dir}/sn_groupmask_in_${sub}.nii" "${base_dir}/sn_groupmask_in_${sub}.nii.gz"
  rm -f "${base_dir}/snctrl_groupmask_in_${sub}.nii" "${base_dir}/snctrl_groupmask_in_${sub}.nii.gz"
  rm -f "${base_dir}/ptctrl_groupmask_in_${sub}.nii" "${base_dir}/ptctrl_groupmask_in_${sub}.nii.gz"

  if [[ -n ${kw_path} ]]; then
    ln -sf "${kw_path}" "${base_dir}/sn_groupmask_in_${sub}.$([[ ${kw_path} == *.nii.gz ]] && echo 'nii.gz' || echo 'nii')"
  else
    notes+=(missing_kw)
  fi
  if [[ -n ${snctrl_src} ]]; then
    ln -sf "${snctrl_src}" "${base_dir}/snctrl_groupmask_in_${sub}.$([[ ${snctrl_src} == *.nii.gz ]] && echo 'nii.gz' || echo 'nii')"
  else
    notes+=(missing_snctrl)
  fi
  if [[ -n ${ptctrl_src} ]]; then
    ln -sf "${ptctrl_src}" "${base_dir}/ptctrl_groupmask_in_${sub}.$([[ ${ptctrl_src} == *.nii.gz ]] && echo 'nii.gz' || echo 'nii')"
  else
    notes+=(missing_ptctrl_base)
  fi

  for variant in "${pt_variants[@]}"; do
    label=${variant^^}
    src_variant_dir="${source_root}/PT_${label}"
    out_variant_dir="${out_root}/PT_${label}"
    mkdir -p "${out_variant_dir}"
    snctrl_variant=$(locate_source_mask "${src_variant_dir}" "snctrl_groupmask_in_${sub}" || true)
    ptctrl_variant=$(locate_source_mask "${src_variant_dir}" "ptctrl_groupmask_in_${sub}" || true)

    rm -f "${out_variant_dir}/sn_groupmask_in_${sub}.nii" "${out_variant_dir}/sn_groupmask_in_${sub}.nii.gz"
    rm -f "${out_variant_dir}/snctrl_groupmask_in_${sub}.nii" "${out_variant_dir}/snctrl_groupmask_in_${sub}.nii.gz"
    rm -f "${out_variant_dir}/ptctrl_groupmask_in_${sub}.nii" "${out_variant_dir}/ptctrl_groupmask_in_${sub}.nii.gz"

    if [[ -n ${kw_path} ]]; then
      ln -sf "${kw_path}" "${out_variant_dir}/sn_groupmask_in_${sub}.$([[ ${kw_path} == *.nii.gz ]] && echo 'nii.gz' || echo 'nii')"
    fi
    if [[ -n ${snctrl_variant} ]]; then
      ln -sf "${snctrl_variant}" "${out_variant_dir}/snctrl_groupmask_in_${sub}.$([[ ${snctrl_variant} == *.nii.gz ]] && echo 'nii.gz' || echo 'nii')"
    else
      notes+=("missing_snctrl_${label}")
    fi
    if [[ -n ${ptctrl_variant} ]]; then
      ln -sf "${ptctrl_variant}" "${out_variant_dir}/ptctrl_groupmask_in_${sub}.$([[ ${ptctrl_variant} == *.nii.gz ]] && echo 'nii.gz' || echo 'nii')"
    else
      notes+=("missing_ptctrl_${label}")
    fi
  done

  status=ok
  if [[ ${#notes[@]} -gt 0 ]]; then
    status=warn
  fi
  note_str=${notes[*]:-}
  note_str=${note_str// /,}
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sub}" "${kw_path:-missing}" "${snctrl_src:-missing}" "${ptctrl_src:-missing}" "${status}" "${note_str:-ok}" >> "${manifest}"
done

echo "[INFO] KW-backed mask root prepared at ${out_root}"
echo "[INFO] Manifest written to ${manifest}"
