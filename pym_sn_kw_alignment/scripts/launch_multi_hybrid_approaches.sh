#!/usr/bin/env bash
set -euo pipefail

ROOT=${SNHEART_ROOT:-/path/to/SNHEART}
TMP_DIR=${ROOT}/tmp
BUILD_SCRIPT=${ROOT}/Scripts/build_black_ant_masks_v10.sh
ALIGN_SCRIPT=${ROOT}/Scripts/check_sn_alignment_metrics.py

usage() {
  cat <<'EOF'
Usage: launch_multi_hybrid_approaches.sh [RUN_DIR]

Environment variables:
  MULTIOPT_RUN_DIR            Explicit run directory (default: positional arg, else newest raw_kw_loop_*).
  MULTIOPT_SUBJECT_MODE       frontier | all | file (default: frontier).
  MULTIOPT_SUBJECTS_FILE      Subject list when MULTIOPT_SUBJECT_MODE=file.
  MULTIOPT_CHUNK_SIZE         Subjects per chunk (default: 12).
  MULTIOPT_NODES              sbatch --nodes (default: 1).
  MULTIOPT_NTASKS             sbatch --ntasks (default: 1).
  MULTIOPT_APPROACH_SPECS     Comma list name:key:scale triplets.
                              Default:
                              strict_full:strict_final:1.0,strict_gain70:strict_final:0.7,strict_gain40:strict_final:0.4
  MULTIOPT_MASK_LABELS        BLACKANT_MASK_LABELS value (default: SN).
  MULTIOPT_ROI_OFFSET_SIGN    BLACKANT_ROI_OFFSET_SIGN (default: 1).
  MULTIOPT_ISOLATE_ATLAS      1=use per-chunk BLACKANT_ATLAS_DIR (default: 1).

Notes:
  - Each approach writes to a separate BLACKANT_OUTPUT_DIR, so approaches can run concurrently.
  - With MULTIOPT_ISOLATE_ATLAS=1, each chunk also uses a private atlas cache directory.
  - Approach key values:
      strict_final | strict | conservative_final | conservative
EOF
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

run_dir_arg=${1:-}
if [[ -n ${MULTIOPT_RUN_DIR:-} ]]; then
  RUN_DIR=${MULTIOPT_RUN_DIR}
elif [[ -n ${run_dir_arg} ]]; then
  RUN_DIR=${run_dir_arg}
else
  RUN_DIR=$(ls -dt "${TMP_DIR}"/raw_kw_loop_* 2>/dev/null | head -n 1 || true)
fi

[[ -n ${RUN_DIR} ]] || { echo "[ERROR] Could not resolve RUN_DIR" >&2; exit 1; }
[[ -d ${RUN_DIR} ]] || { echo "[ERROR] RUN_DIR does not exist: ${RUN_DIR}" >&2; exit 1; }
[[ -x ${BUILD_SCRIPT} ]] || { echo "[ERROR] Missing/executable build script: ${BUILD_SCRIPT}" >&2; exit 1; }
[[ -f ${ALIGN_SCRIPT} ]] || { echo "[ERROR] Missing align checker: ${ALIGN_SCRIPT}" >&2; exit 1; }
command -v sbatch >/dev/null 2>&1 || { echo "[ERROR] sbatch not found" >&2; exit 1; }
command -v python >/dev/null 2>&1 || command -v python3 >/dev/null 2>&1 || { echo "[ERROR] python not found" >&2; exit 1; }

PYTHON_BIN=${PYTHON_BIN:-python3}
if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
  PYTHON_BIN=python
fi

SUBJECT_MODE=${MULTIOPT_SUBJECT_MODE:-frontier}
SUBJECTS_FILE=${MULTIOPT_SUBJECTS_FILE:-}
CHUNK_SIZE=${MULTIOPT_CHUNK_SIZE:-12}
NODES=${MULTIOPT_NODES:-1}
NTASKS=${MULTIOPT_NTASKS:-1}
MASK_LABELS=${MULTIOPT_MASK_LABELS:-SN}
ROI_OFFSET_SIGN=${MULTIOPT_ROI_OFFSET_SIGN:-1}
ISOLATE_ATLAS=${MULTIOPT_ISOLATE_ATLAS:-1}
APPROACH_SPECS=${MULTIOPT_APPROACH_SPECS:-strict_full:strict_final:1.0,strict_gain70:strict_final:0.7,strict_gain40:strict_final:0.4}

highhat_all=${RUN_DIR}/highhat/highhat_rawkw_subjects.txt
crash_all=${RUN_DIR}/crash/crash_rawkw_subjects.txt
highhat_pre=${RUN_DIR}/highhat/highhat_rawkw_vs_v10_pre.tsv
crash_pre=${RUN_DIR}/crash/crash_rawkw_vs_v10_pre.tsv
highhat_manual=${RUN_DIR}/highhat/manual_raw_kw
crash_manual=${RUN_DIR}/crash/manual_raw_kw

for f in "${highhat_all}" "${crash_all}" "${highhat_pre}" "${crash_pre}"; do
  [[ -f ${f} ]] || { echo "[ERROR] Missing required run artifact: ${f}" >&2; exit 1; }
done
for d in "${highhat_manual}" "${crash_manual}"; do
  [[ -d ${d} ]] || { echo "[ERROR] Missing required manual dir: ${d}" >&2; exit 1; }
done

resolve_json_for_key() {
  local key=$1
  if [[ ${key} == file=* ]]; then
    echo "${key#file=}"
    return 0
  fi
  case "${key}" in
    strict_final)
      if [[ -f ${RUN_DIR}/hybrid_offsets_strict_final.json ]]; then
        echo "${RUN_DIR}/hybrid_offsets_strict_final.json"
      else
        echo "${RUN_DIR}/hybrid_offsets_strict.json"
      fi
      ;;
    strict)
      echo "${RUN_DIR}/hybrid_offsets_strict.json"
      ;;
    conservative_final)
      if [[ -f ${RUN_DIR}/hybrid_offsets_conservative_final.json ]]; then
        echo "${RUN_DIR}/hybrid_offsets_conservative_final.json"
      else
        echo "${RUN_DIR}/hybrid_offsets_conservative.json"
      fi
      ;;
    conservative)
      echo "${RUN_DIR}/hybrid_offsets_conservative.json"
      ;;
    *)
      echo ""
      ;;
  esac
}

selected_all=${RUN_DIR}/multiopt_selected_subjects.txt
case "${SUBJECT_MODE}" in
  frontier)
    snapshot_tsv=${RUN_DIR}/hybrid_live_direction_snapshot.tsv
    [[ -f ${snapshot_tsv} ]] || { echo "[ERROR] frontier mode needs ${snapshot_tsv}" >&2; exit 1; }
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} {d=$(h["direction"]); s=$(h["subject"]); if (d!="closer" && s!="") print s}' \
      "${snapshot_tsv}" | sort -u -n > "${selected_all}"
    ;;
  all)
    cat "${highhat_all}" "${crash_all}" | awk 'NF{print $1}' | sort -u -n > "${selected_all}"
    ;;
  file)
    [[ -n ${SUBJECTS_FILE} ]] || { echo "[ERROR] MULTIOPT_SUBJECTS_FILE is required for mode=file" >&2; exit 1; }
    [[ -f ${SUBJECTS_FILE} ]] || { echo "[ERROR] Subject file not found: ${SUBJECTS_FILE}" >&2; exit 1; }
    awk 'NF{print $1}' "${SUBJECTS_FILE}" | sort -u -n > "${selected_all}"
    ;;
  *)
    echo "[ERROR] Unsupported MULTIOPT_SUBJECT_MODE=${SUBJECT_MODE} (use frontier|all|file)" >&2
    exit 1
    ;;
esac

if [[ ! -s ${selected_all} ]]; then
  echo "[ERROR] No subjects selected (mode=${SUBJECT_MODE})" >&2
  exit 1
fi

multiopt_dir=${RUN_DIR}/multiopt_outputs
offset_dir=${RUN_DIR}/multiopt_offsets
mkdir -p "${multiopt_dir}" "${offset_dir}"

summary_tsv=${RUN_DIR}/multiopt_submission_summary.tsv
printf "approach\toffset_json\tgroup\tsubjects_file\tbuild_jobs\tpostcheck_job\tbaseline_dir\tpost_tsv\n" > "${summary_tsv}"

make_scaled_json() {
  local src_json=$1
  local scale=$2
  local out_json=$3
  local approach_name=$4
  "${PYTHON_BIN}" - "${src_json}" "${scale}" "${out_json}" "${approach_name}" <<'PY'
import json
import sys

src, scale_raw, out, approach = sys.argv[1:5]
scale = float(scale_raw)
with open(src, encoding="utf-8") as f:
    data = json.load(f)
for sub, entry in data.items():
    delta = entry.get("delta_mm")
    if not isinstance(delta, list) or len(delta) != 3:
        continue
    try:
        vals = [float(v) * scale for v in delta]
    except Exception:
        continue
    entry["delta_mm"] = vals
    entry["multiopt_scale"] = scale
    entry["multiopt_approach"] = approach
with open(out, "w", encoding="utf-8") as f:
    json.dump(data, f, indent=2, sort_keys=True)
PY
}

submit_group_chunks() {
  local approach=$1
  local group=$2
  local subjects_file=$3
  local offset_json=$4
  local output_dir=$5
  local approach_dir=$6
  local jobs_file=${approach_dir}/${group}_build_jobids.txt
  local atlas_root=${approach_dir}/atlas_cache
  : > "${jobs_file}"
  mapfile -t subjects < <(awk 'NF{print $1}' "${subjects_file}")
  local total=${#subjects[@]}
  if (( total == 0 )); then
    echo ""
    return 0
  fi
  local i=0
  local chunk=1
  while (( i < total )); do
    local chunk_file=${approach_dir}/${group}_chunk_${chunk}.txt
    printf "%s\n" "${subjects[@]:i:CHUNK_SIZE}" > "${chunk_file}"
    local chunk_subjects
    chunk_subjects=$(paste -sd' ' "${chunk_file}")
    local job_name="mopt_${approach}_${group}_c${chunk}"
    local log_file=${approach_dir}/${job_name}.log
    local chunk_atlas_dir="${atlas_root}/${group}_chunk_${chunk}"
    mkdir -p "${chunk_atlas_dir}"
    local jid
    if [[ ${ISOLATE_ATLAS} -eq 1 ]]; then
      jid=$(BLACKANT_MASK_SUBJECTS="${chunk_subjects}" \
        BLACKANT_PREFER_ANTPUNK=0 \
        BLACKANT_ROI_OFFSET_JSON="${offset_json}" \
        BLACKANT_APPLY_ROI_OFFSETS=1 \
        BLACKANT_ROI_OFFSET_SIGN="${ROI_OFFSET_SIGN}" \
        BLACKANT_MASK_LABELS="${MASK_LABELS}" \
        BLACKANT_OUTPUT_DIR="${output_dir}" \
        BLACKANT_ATLAS_DIR="${chunk_atlas_dir}" \
        sbatch --parsable --nodes="${NODES}" --ntasks="${NTASKS}" \
        --job-name="${job_name}" --output="${log_file}" "${BUILD_SCRIPT}")
    else
      jid=$(BLACKANT_MASK_SUBJECTS="${chunk_subjects}" \
        BLACKANT_PREFER_ANTPUNK=0 \
        BLACKANT_ROI_OFFSET_JSON="${offset_json}" \
        BLACKANT_APPLY_ROI_OFFSETS=1 \
        BLACKANT_ROI_OFFSET_SIGN="${ROI_OFFSET_SIGN}" \
        BLACKANT_MASK_LABELS="${MASK_LABELS}" \
        BLACKANT_OUTPUT_DIR="${output_dir}" \
        sbatch --parsable --nodes="${NODES}" --ntasks="${NTASKS}" \
        --job-name="${job_name}" --output="${log_file}" "${BUILD_SCRIPT}")
    fi
    echo "${jid}" >> "${jobs_file}"
    printf "[INFO] submitted %s -> %s (%s) atlas=%s\n" "${job_name}" "${jid}" "${chunk_file}" "${chunk_atlas_dir}" >&2
    i=$(( i + CHUNK_SIZE ))
    chunk=$(( chunk + 1 ))
  done
  paste -sd, "${jobs_file}"
}

submit_postcheck() {
  local approach=$1
  local group=$2
  local dep_csv=$3
  local subjects_file=$4
  local manual_dir=$5
  local baseline_dir=$6
  local out_tsv=$7
  local approach_dir=$8
  if [[ -z ${dep_csv} ]]; then
    echo ""
    return 0
  fi
  local dep="afterok:${dep_csv//,/:}"
  local job_name="mopt_qc_${approach}_${group}"
  local log_file=${approach_dir}/${job_name}.log
  local cmd
  cmd="cd ${ROOT}; ${PYTHON_BIN} ${ALIGN_SCRIPT} --subjects-file ${subjects_file} --manual-dir ${manual_dir} --baseline-dir ${baseline_dir} --output ${out_tsv}"
  sbatch --parsable --dependency="${dep}" --job-name="${job_name}" --output="${log_file}" --wrap="${cmd}"
}

IFS=',' read -r -a specs <<< "${APPROACH_SPECS}"
if (( ${#specs[@]} == 0 )); then
  echo "[ERROR] No approaches parsed from MULTIOPT_APPROACH_SPECS" >&2
  exit 1
fi

for spec in "${specs[@]}"; do
  IFS=':' read -r approach_name source_key scale <<< "${spec}"
  [[ -n ${approach_name} && -n ${source_key} && -n ${scale} ]] || {
    echo "[ERROR] Invalid approach spec: ${spec} (expected name:key:scale)" >&2
    exit 1
  }
  src_json=$(resolve_json_for_key "${source_key}")
  [[ -n ${src_json} ]] || { echo "[ERROR] Unknown source key in spec: ${spec}" >&2; exit 1; }
  [[ -f ${src_json} ]] || { echo "[ERROR] Missing source JSON for ${spec}: ${src_json}" >&2; exit 1; }

  offset_json=${offset_dir}/${approach_name}.json
  make_scaled_json "${src_json}" "${scale}" "${offset_json}" "${approach_name}"

  approach_dir=${multiopt_dir}/${approach_name}
  output_dir=${approach_dir}/probabilistic_masks_v10
  baseline_dir=${output_dir}/base
  mkdir -p "${approach_dir}" "${output_dir}" "${baseline_dir}"

  highhat_selected=${approach_dir}/highhat_subjects.txt
  crash_selected=${approach_dir}/crash_subjects.txt
  awk 'NR==FNR{s[$1]=1;next} ($1 in s){print $1}' "${selected_all}" "${highhat_all}" | sort -u -n > "${highhat_selected}"
  awk 'NR==FNR{s[$1]=1;next} ($1 in s){print $1}' "${selected_all}" "${crash_all}" | sort -u -n > "${crash_selected}"

  highhat_jobs=$(submit_group_chunks "${approach_name}" highhat "${highhat_selected}" "${offset_json}" "${output_dir}" "${approach_dir}")
  crash_jobs=$(submit_group_chunks "${approach_name}" crash "${crash_selected}" "${offset_json}" "${output_dir}" "${approach_dir}")

  highhat_post_tsv=${approach_dir}/highhat_post.tsv
  crash_post_tsv=${approach_dir}/crash_post.tsv
  highhat_post_job=$(submit_postcheck "${approach_name}" highhat "${highhat_jobs}" "${highhat_selected}" "${highhat_manual}" "${baseline_dir}" "${highhat_post_tsv}" "${approach_dir}")
  crash_post_job=$(submit_postcheck "${approach_name}" crash "${crash_jobs}" "${crash_selected}" "${crash_manual}" "${baseline_dir}" "${crash_post_tsv}" "${approach_dir}")

  printf "%s\t%s\thighhat\t%s\t%s\t%s\t%s\t%s\n" \
    "${approach_name}" "${offset_json}" "${highhat_selected}" "${highhat_jobs}" "${highhat_post_job}" "${baseline_dir}" "${highhat_post_tsv}" >> "${summary_tsv}"
  printf "%s\t%s\tcrash\t%s\t%s\t%s\t%s\t%s\n" \
    "${approach_name}" "${offset_json}" "${crash_selected}" "${crash_jobs}" "${crash_post_job}" "${baseline_dir}" "${crash_post_tsv}" >> "${summary_tsv}"

  echo "[DONE] approach=${approach_name} source=${source_key} scale=${scale} offset_json=${offset_json}"
  echo "       highhat_subjects=$(wc -l < "${highhat_selected}") crash_subjects=$(wc -l < "${crash_selected}")"
done

echo "[DONE] Multi-approach submissions complete"
echo "[DONE] summary=${summary_tsv}"
echo "[HINT] live snapshot: bash ${ROOT}/Scripts/snapshot_multi_hybrid_direction.sh ${RUN_DIR}"
