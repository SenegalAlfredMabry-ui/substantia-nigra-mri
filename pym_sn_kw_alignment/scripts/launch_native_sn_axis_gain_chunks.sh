#!/usr/bin/env bash
set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
shared_root=${SNHEART_ROOT:-/path/to/SNHEART}
run_dir=${AXIS_GAIN_RUN_DIR:-$(ls -dt "${shared_root}"/tmp/raw_kw_loop_* 2>/dev/null | head -n 1 || true)}
[[ -n ${run_dir} && -d ${run_dir} ]] || { echo "[ERROR] Could not resolve raw KW run dir" >&2; exit 1; }

prev_root=${AXIS_GAIN_PREV_ROOT:-$(ls -dt "${shared_root}"/tmp/native_sn_translation_gain_* 2>/dev/null | head -n 1 || true)}
summary_file=${AXIS_GAIN_SUMMARY_FILE:-${prev_root}/best_gain_per_subject_merged.tsv}
[[ -f ${summary_file} ]] || { echo "[ERROR] Could not resolve previous translation summary" >&2; exit 1; }

out_root=${AXIS_GAIN_ROOT:-${repo_root}/tmp/native_sn_axis_gain_$(date -u +%Y%m%dT%H%M%SZ)}
chunk_count=${AXIS_GAIN_CHUNK_COUNT:-4}
gain_list=${AXIS_GAIN_GAINS:-"-0.5 -0.25 0 0.25 0.5 0.75 1.0"}
runner=${repo_root}/Scripts/run_native_sn_axis_gain_chunk.sh
prep=${repo_root}/Scripts/prepare_native_sn_rigid_rescue_inputs.py
subjects_file=${AXIS_GAIN_SUBJECTS_FILE:-}

mkdir -p "${out_root}"
python3 "${prep}" --run-dir "${run_dir}" --output-root "${out_root}/inputs"

input_dir="${out_root}/inputs/staged_masks"
manual_dir="${out_root}/inputs/manual_masks"
manifest="${out_root}/inputs/starting_mask_manifest.tsv"

if [[ -z ${subjects_file} ]]; then
  subjects_file="${out_root}/subjects_from_pre_nonimprovers.txt"
  python3 - "${summary_file}" "${subjects_file}" <<'PY'
import csv
import sys
from pathlib import Path

summary_path = Path(sys.argv[1])
subjects_path = Path(sys.argv[2])
rows = list(csv.DictReader(summary_path.open(), delimiter="\t"))
subjects = [row["subject"] for row in rows if row.get("best_stage") == "pre"]
subjects_path.write_text("\n".join(subjects) + "\n", encoding="utf-8")
print(subjects_path)
print(f"count={len(subjects)}")
PY
fi

resolved_subjects="${out_root}/resolved_subjects.txt"
awk 'NR==FNR{keep[$1]=1; next} NR>1 && ($1 in keep) && $3!="missing"{print $1}' \
  "${subjects_file}" "${manifest}" | sort -n > "${resolved_subjects}"

mapfile -t subjects < "${resolved_subjects}"
[[ ${#subjects[@]} -gt 0 ]] || { echo "[ERROR] No staged subjects resolved for axis-gain rescue" >&2; exit 1; }

summary="${out_root}/submission_summary.tsv"
printf "chunk\tjobid\tcount\tgains\tsubjects_file\toutput_root\n" > "${summary}"
chunk_size=$(( (${#subjects[@]} + chunk_count - 1) / chunk_count ))

i=0
chunk=1
while (( i < ${#subjects[@]} )); do
  tag=$(printf "c%02d" "${chunk}")
  chunk_dir="${out_root}/chunks/${tag}"
  mkdir -p "${chunk_dir}"
  chunk_subjects_file="${chunk_dir}/subjects.txt"
  printf "%s\n" "${subjects[@]:i:chunk_size}" > "${chunk_subjects_file}"
  job_name="axis_gain_${tag}"
  jid=$(sbatch --parsable --job-name="${job_name}" \
    --output="${chunk_dir}/${job_name}.out" --error="${chunk_dir}/${job_name}.err" \
    --export=ALL,AXIS_GAIN_SUBJECTS_FILE="${chunk_subjects_file}",AXIS_GAIN_MANUAL_DIR="${manual_dir}",AXIS_GAIN_INPUT_DIR="${input_dir}",AXIS_GAIN_OUTPUT_ROOT="${chunk_dir}",AXIS_GAIN_GAINS="${gain_list}" \
    "${runner}")
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${tag}" "${jid}" "$(wc -l < "${chunk_subjects_file}")" "${gain_list}" "${chunk_subjects_file}" "${chunk_dir}" >> "${summary}"
  i=$(( i + chunk_size ))
  chunk=$(( chunk + 1 ))
done

echo "[INFO] Submission summary written to ${summary}"
