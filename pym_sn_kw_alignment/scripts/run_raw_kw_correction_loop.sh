#!/usr/bin/env bash
set -euo pipefail

ROOT=${SNHEART_ROOT:-/path/to/SNHEART}
TMP_DIR=${ROOT}/tmp
MOTIP_SN=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/TSE/SN
ALIGN_SCRIPT=${ROOT}/Scripts/check_sn_alignment_metrics.py
BUILD_SCRIPT=${ROOT}/Scripts/build_black_ant_masks_v10.sh
V10_BASE=${ROOT}/MRI_data/TSE/probabilistic_masks_v10/base
BUILD_CHUNK_SIZE=${RAWKW_BUILD_CHUNK_SIZE:-18}
ROI_OFFSET_SIGN=${RAWKW_ROI_OFFSET_SIGN:--1}
BUILD_NODES=${RAWKW_BUILD_NODES:-1}
BUILD_NTASKS=${RAWKW_BUILD_NTASKS:-1}
RUN_TS=$(date -u +%Y%m%dT%H%M%SZ)
RUN_DIR=${TMP_DIR}/raw_kw_loop_${RUN_TS}

HIGHHAT_WITH_KW=${TMP_DIR}/highhat_with_kw_subjects.txt
HIGHHAT_NO_KW=${TMP_DIR}/highhat_no_kw.txt
CRASH_OFFSETS=${TMP_DIR}/crash_offsets_subjects.txt
CRASH_RIDE=${TMP_DIR}/crash_ride_subjects.txt

mkdir -p "${RUN_DIR}"

require_file() {
  local path=$1
  [[ -f ${path} ]] || { echo "[ERROR] Missing required file: ${path}" >&2; exit 1; }
}

require_file "${ALIGN_SCRIPT}"
require_file "${BUILD_SCRIPT}"
require_file "${HIGHHAT_WITH_KW}"
require_file "${CRASH_OFFSETS}"
[[ -d ${MOTIP_SN} ]] || { echo "[ERROR] Missing MOTIP KW directory: ${MOTIP_SN}" >&2; exit 1; }
[[ -d ${V10_BASE} ]] || { echo "[ERROR] Missing v10 base directory: ${V10_BASE}" >&2; exit 1; }
command -v sbatch >/dev/null 2>&1 || { echo "[ERROR] sbatch is not available" >&2; exit 1; }

build_union_roster() {
  local out_file=$1
  shift
  : > "${out_file}"
  for src in "$@"; do
    [[ -f ${src} ]] || continue
    awk 'NF {print $1}' "${src}" >> "${out_file}"
  done
  sort -u -n "${out_file}" -o "${out_file}"
}

find_kw_mask() {
  local sub=$1
  local candidates=(
    "${MOTIP_SN}/SN_ROI_KW_${sub}_ZP.nii"
    "${MOTIP_SN}/SN_ROI_KW_${sub}_ZP.nii.gz"
    "${MOTIP_SN}/SN_ROI_KW_${sub}.nii"
    "${MOTIP_SN}/SN_ROI_KW_${sub}.nii.gz"
  )
  for c in "${candidates[@]}"; do
    if [[ -f ${c} ]]; then
      printf '%s\n' "${c}"
      return 0
    fi
  done
  return 1
}

prepare_group_inputs() {
  local group=$1
  local roster=$2
  local group_dir=${RUN_DIR}/${group}
  local manual_dir=${group_dir}/manual_raw_kw
  local valid_file=${group_dir}/${group}_rawkw_subjects.txt
  local missing_file=${group_dir}/${group}_missing_or_skipped.tsv
  local pre_tsv=${group_dir}/${group}_rawkw_vs_v10_pre.tsv
  local offset_json=${group_dir}/${group}_offsets_rawkw_v10.json
  local offset_summary=${group_dir}/${group}_offsets_summary.tsv

  mkdir -p "${group_dir}" "${manual_dir}"
  : > "${valid_file}"
  printf "subject\treason\n" > "${missing_file}"

  while read -r sub; do
    [[ -z ${sub} ]] && continue
    local kw_path=""
    kw_path=$(find_kw_mask "${sub}" || true)
    local baseline=${V10_BASE}/sn_groupmask_in_${sub}.nii.gz
    [[ -f ${baseline} ]] || baseline=${V10_BASE}/sn_groupmask_in_${sub}.nii

    if [[ -z ${kw_path} ]]; then
      printf "%s\tmissing_raw_kw\n" "${sub}" >> "${missing_file}"
      continue
    fi
    if [[ ! -f ${baseline} ]]; then
      printf "%s\tmissing_v10_mask\n" "${sub}" >> "${missing_file}"
      continue
    fi
    ln -sf "${kw_path}" "${manual_dir}/sn_groupmask_in_${sub}.nii"
    printf "%s\n" "${sub}" >> "${valid_file}"
  done < "${roster}"

  if [[ ! -s ${valid_file} ]]; then
    echo "[ERROR] ${group}: no valid raw-KW subjects available after filtering" >&2
    exit 1
  fi

  python "${ALIGN_SCRIPT}" \
    --subjects-file "${valid_file}" \
    --manual-dir "${manual_dir}" \
    --baseline-dir "${V10_BASE}" \
    --output "${pre_tsv}"

  python - "${pre_tsv}" "${offset_json}" "${offset_summary}" <<'PY'
import csv
import json
import math
import sys

in_tsv, out_json, out_summary = sys.argv[1:4]
offsets = {}
rows = []

with open(in_tsv, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        sub = r["subject"]
        notes = r.get("notes", "")
        if "manual_missing" in notes or "baseline_missing" in notes:
            continue
        try:
            dx = float(r["delta_x_mm"])
            dy = float(r["delta_y_mm"])
            dz = float(r["delta_z_mm"])
            dm = float(r["delta_mag_mm"])
        except ValueError:
            continue
        if any(math.isnan(v) for v in (dx, dy, dz, dm)):
            continue
        offsets[sub] = {
            "delta_mm": [dx, dy, dz],
            "delta_mag_mm": dm,
        }
        rows.append((sub, dx, dy, dz, dm))

with open(out_json, "w", encoding="utf-8") as f:
    json.dump(offsets, f, indent=2, sort_keys=True)

with open(out_summary, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["subject", "delta_x_mm", "delta_y_mm", "delta_z_mm", "delta_mag_mm"])
    for row in rows:
        w.writerow(row)
PY

  echo "${manual_dir}|${valid_file}|${pre_tsv}|${offset_json}|${offset_summary}|${missing_file}"
}

join_subjects() {
  local file=$1
  paste -sd' ' "${file}"
}

submit_build_job() {
  local job_name=$1
  local subjects_file=$2
  local offset_json=$3
  local dep=${4:-}
  local log_file=${RUN_DIR}/${job_name}.log
  local subject_line
  subject_line=$(join_subjects "${subjects_file}")
  local -a dep_args=()
  if [[ -n ${dep} ]]; then
    dep_args+=(--dependency="${dep}")
  fi
  BLACKANT_MASK_SUBJECTS="${subject_line}" \
  BLACKANT_PREFER_ANTPUNK=0 \
  BLACKANT_ROI_OFFSET_JSON="${offset_json}" \
  BLACKANT_APPLY_ROI_OFFSETS=1 \
  BLACKANT_ROI_OFFSET_SIGN="${ROI_OFFSET_SIGN}" \
  sbatch --parsable "${dep_args[@]}" --nodes="${BUILD_NODES}" --ntasks="${BUILD_NTASKS}" \
    --job-name="${job_name}" --output="${log_file}" "${BUILD_SCRIPT}"
}

submit_build_chunks() {
  local group=$1
  local subjects_file=$2
  local offset_json=$3
  local chunk_size=$4
  mapfile -t subjects < <(awk 'NF {print $1}' "${subjects_file}")
  local total=${#subjects[@]}
  if (( total == 0 )); then
    echo "[ERROR] ${group}: no subjects to submit for build" >&2
    exit 1
  fi
  local jobids=()
  local i=0
  local chunk_idx=1
  while (( i < total )); do
    local chunk_file=${RUN_DIR}/${group}_build_chunk_${chunk_idx}_subjects.txt
    printf "%s\n" "${subjects[@]:i:chunk_size}" > "${chunk_file}"
    local job_name="build_v10_${group}_rawkw_sswarper_chunk${chunk_idx}"
    local jid
    jid=$(submit_build_job "${job_name}" "${chunk_file}" "${offset_json}")
    jobids+=("${jid}")
    printf "[INFO] submitted %s job=%s subjects_file=%s\n" "${job_name}" "${jid}" "${chunk_file}" >&2
    i=$(( i + chunk_size ))
    chunk_idx=$(( chunk_idx + 1 ))
  done
  local jobs_file=${RUN_DIR}/${group}_build_jobids.txt
  printf "%s\n" "${jobids[@]}" > "${jobs_file}"
  local dep="afterok:$(IFS=:; echo "${jobids[*]}")"
  local jid_csv
  jid_csv=$(IFS=,; echo "${jobids[*]}")
  echo "${dep}|${jid_csv}"
}

submit_postcheck_job() {
  local job_name=$1
  local dep=$2
  local subjects_file=$3
  local manual_dir=$4
  local output_tsv=$5
  local log_file=${RUN_DIR}/${job_name}.log
  local cmd
  cmd="cd ${ROOT}; python ${ALIGN_SCRIPT} --subjects-file ${subjects_file} --manual-dir ${manual_dir} --baseline-dir ${V10_BASE} --output ${output_tsv}"
  sbatch --parsable --dependency="${dep}" --job-name="${job_name}" --output="${log_file}" --wrap="${cmd}"
}

highhat_union=${RUN_DIR}/highhat_all_subjects.txt
crash_union=${RUN_DIR}/crash_all_subjects.txt
build_union_roster "${highhat_union}" "${HIGHHAT_WITH_KW}" "${HIGHHAT_NO_KW}"
build_union_roster "${crash_union}" "${CRASH_OFFSETS}" "${CRASH_RIDE}"

IFS='|' read -r highhat_manual highhat_valid highhat_pre highhat_json highhat_offsets_summary highhat_missing \
  <<< "$(prepare_group_inputs highhat "${highhat_union}")"
IFS='|' read -r crash_manual crash_valid crash_pre crash_json crash_offsets_summary crash_missing \
  <<< "$(prepare_group_inputs crash "${crash_union}")"

IFS='|' read -r highhat_build_dep highhat_build_jobs \
  <<< "$(submit_build_chunks highhat "${highhat_valid}" "${highhat_json}" "${BUILD_CHUNK_SIZE}")"
IFS='|' read -r crash_build_dep crash_build_jobs \
  <<< "$(submit_build_chunks crash "${crash_valid}" "${crash_json}" "${BUILD_CHUNK_SIZE}")"

highhat_post_tsv=${RUN_DIR}/highhat_rawkw_vs_v10_post.tsv
highhat_post_jid=$(submit_postcheck_job qc_highhat_rawkw_post "${highhat_build_dep}" "${highhat_valid}" "${highhat_manual}" "${highhat_post_tsv}")

crash_post_tsv=${RUN_DIR}/crash_rawkw_vs_v10_post.tsv
crash_post_jid=$(submit_postcheck_job qc_crash_rawkw_post "${crash_build_dep}" "${crash_valid}" "${crash_manual}" "${crash_post_tsv}")

cat > "${RUN_DIR}/run_summary.tsv" <<EOF
group	pre_metrics	offset_json	missing_report	valid_subjects	warp_mode	build_jobs	postcheck_job	post_metrics
highhat	${highhat_pre}	${highhat_json}	${highhat_missing}	${highhat_valid}	sswarper_only	${highhat_build_jobs}	${highhat_post_jid}	${highhat_post_tsv}
crash	${crash_pre}	${crash_json}	${crash_missing}	${crash_valid}	sswarper_only	${crash_build_jobs}	${crash_post_jid}	${crash_post_tsv}
EOF

echo "[DONE] Raw-KW correction loop started."
echo "[DONE] Run directory: ${RUN_DIR}"
echo "[DONE] Summary: ${RUN_DIR}/run_summary.tsv"
