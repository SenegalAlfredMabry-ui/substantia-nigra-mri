#!/usr/bin/env bash
set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
shared_root=${SNHEART_ROOT:-/path/to/SNHEART}
run_dir=${RIGID_RESCUE_RUN_DIR:-$(ls -dt "${shared_root}"/tmp/raw_kw_loop_* 2>/dev/null | head -n 1 || true)}
[[ -n ${run_dir} && -d ${run_dir} ]] || { echo "[ERROR] Could not resolve run dir" >&2; exit 1; }

out_root=${RIGID_RESCUE_ROOT:-${repo_root}/tmp/native_sn_translation_gain_$(date -u +%Y%m%dT%H%M%SZ)}
chunk_count=${RIGID_RESCUE_CHUNK_COUNT:-4}
gain_list=${RIGID_RESCUE_GAINS:-"0.5 1.0 1.5 2.0"}
runner=${repo_root}/Scripts/run_native_sn_translation_gain_chunk.sh
prep=${repo_root}/Scripts/prepare_native_sn_rigid_rescue_inputs.py

mkdir -p "${out_root}"
python3 "${prep}" --run-dir "${run_dir}" --output-root "${out_root}/inputs"

input_dir="${out_root}/inputs/staged_masks"
manual_dir="${out_root}/inputs/manual_masks"
manifest="${out_root}/inputs/starting_mask_manifest.tsv"
summary="${out_root}/submission_summary.tsv"
printf "chunk\tjobid\tcount\tgains\tsubjects_file\toutput_root\n" > "${summary}"

mapfile -t subjects < <(awk -F'\t' 'NR>1 && $3!="missing"{print $1}' "${manifest}" | sort -n)
[[ ${#subjects[@]} -gt 0 ]] || { echo "[ERROR] No staged subjects" >&2; exit 1; }
chunk_size=$(( (${#subjects[@]} + chunk_count - 1) / chunk_count ))

i=0
chunk=1
while (( i < ${#subjects[@]} )); do
  tag=$(printf "c%02d" "${chunk}")
  chunk_dir="${out_root}/chunks/${tag}"
  mkdir -p "${chunk_dir}"
  subjects_file="${chunk_dir}/subjects.txt"
  printf "%s\n" "${subjects[@]:i:chunk_size}" > "${subjects_file}"
  job_name="gain_sn_${tag}"
  jid=$(sbatch --parsable --job-name="${job_name}" \
    --output="${chunk_dir}/${job_name}.out" --error="${chunk_dir}/${job_name}.err" \
    --export=ALL,RIGID_RESCUE_SUBJECTS_FILE="${subjects_file}",RIGID_RESCUE_MANUAL_DIR="${manual_dir}",RIGID_RESCUE_INPUT_DIR="${input_dir}",RIGID_RESCUE_OUTPUT_ROOT="${chunk_dir}",RIGID_RESCUE_GAINS="${gain_list}" \
    "${runner}")
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${tag}" "${jid}" "$(wc -l < "${subjects_file}")" "${gain_list}" "${subjects_file}" "${chunk_dir}" >> "${summary}"
  i=$(( i + chunk_size ))
  chunk=$(( chunk + 1 ))
done

echo "[INFO] Submission summary written to ${summary}"
