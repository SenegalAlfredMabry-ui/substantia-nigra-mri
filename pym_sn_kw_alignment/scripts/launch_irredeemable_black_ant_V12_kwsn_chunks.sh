#!/usr/bin/env bash
set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
shared_root=${SNHEART_ROOT:-/path/to/SNHEART}
runner=${repo_root}/Scripts/run_irredeemable_black_ant_V12_kwsn.sh

usage() {
  cat <<'EOF'
Usage: launch_irredeemable_black_ant_V12_kwsn_chunks.sh

Environment:
  V12_KWSN_RUN_DIR       raw_kw_loop run dir (default: newest under ${MOTIP_ROOT}/.../tmp)
  V12_KWSN_SUBJECTS      Optional explicit subject list. Default: raw-KW cohort from run dir
  V12_KWSN_CHUNK_COUNT   Number of chunks (default: 4)
  V12_KWSN_CHUNK_SIZE    Explicit chunk size override (takes precedence over chunk count)
  V12_KWSN_PT_VARIANTS   PT variants for the run (default: KW)
  V12_KWSN_CORR_MODES    Corr modes for the run (default: legacy)
  V12_KWSN_ROOT          Root output dir (default: MRI_data/analysis/tests/v12_kwsn_chunks)
  V12_KWSN_DRY_RUN       1=print commands only, 0=submit (default: 0)
EOF
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

run_dir=${V12_KWSN_RUN_DIR:-$(ls -dt "${shared_root}"/tmp/raw_kw_loop_* 2>/dev/null | head -n 1 || true)}
[[ -n ${run_dir} && -d ${run_dir} ]] || { echo "[ERROR] Could not resolve V12_KWSN_RUN_DIR" >&2; exit 1; }
[[ -x ${runner} ]] || { echo "[ERROR] Missing runner: ${runner}" >&2; exit 1; }
command -v sbatch >/dev/null 2>&1 || { echo "[ERROR] sbatch not found" >&2; exit 1; }

resolve_subjects() {
  local hh=${run_dir}/highhat/highhat_rawkw_subjects.txt
  local cr=${run_dir}/crash/crash_rawkw_subjects.txt
  awk 'NF{print $1}' "${hh}" "${cr}" | sort -u -n
}

mapfile -t subjects < <(
  if [[ -n ${V12_KWSN_SUBJECTS:-} ]]; then
    tr ' ' '\n' <<<"${V12_KWSN_SUBJECTS}" | awk 'NF{print $1}' | sort -u -n
  else
    resolve_subjects
  fi
)
[[ ${#subjects[@]} -gt 0 ]] || { echo "[ERROR] No subjects resolved" >&2; exit 1; }

chunk_count=${V12_KWSN_CHUNK_COUNT:-4}
chunk_size=${V12_KWSN_CHUNK_SIZE:-}
if [[ -z ${chunk_size} ]]; then
  total=${#subjects[@]}
  chunk_size=$(( (total + chunk_count - 1) / chunk_count ))
fi

pt_variants=${V12_KWSN_PT_VARIANTS:-KW}
corr_modes=${V12_KWSN_CORR_MODES:-legacy}
root_out=${V12_KWSN_ROOT:-${repo_root}/MRI_data/analysis/tests/v12_kwsn_chunks}
chunk_root=${root_out}/chunks
mkdir -p "${chunk_root}"

summary=${root_out}/submission_summary.tsv
printf "chunk\tcount\tsubjects_file\tanalysis_root\tmask_root\tjobid\n" > "${summary}"

dry_run=${V12_KWSN_DRY_RUN:-0}
i=0
chunk=1
while (( i < ${#subjects[@]} )); do
  chunk_tag=$(printf "c%02d" "${chunk}")
  chunk_dir=${chunk_root}/${chunk_tag}
  mkdir -p "${chunk_dir}"
  subjects_file=${chunk_dir}/subjects.txt
  printf "%s\n" "${subjects[@]:i:chunk_size}" > "${subjects_file}"
  subject_list=$(paste -sd' ' "${subjects_file}")
  count=$(wc -w < <(echo "${subject_list}"))
  analysis_root=${chunk_dir}/analysis
  mask_root=${chunk_dir}/masks
  job_name="v12kwsn_${chunk_tag}"
  cmd=(sbatch --parsable --job-name="${job_name}" \
    --output="${chunk_dir}/${job_name}.out" --error="${chunk_dir}/${job_name}.err" \
    --export=ALL,IRREDEEMABLE_SUBJECTS="${subject_list}",IRREDEEMABLE_EXPECT_COUNT="${count}",IRREDEEMABLE_CHUNK_TAG="${chunk_tag}",IRREDEEMABLE_ANALYSIS_ROOT="${analysis_root}",KWBACKUP_OUTPUT_ROOT="${mask_root}",KWBACKUP_RUN_DIR="${run_dir}",IRREDEEMABLE_PT_VARIANTS="${pt_variants}",IRREDEEMABLE_CORR_MODES="${corr_modes}" \
    "${runner}")
  if [[ ${dry_run} -eq 1 ]]; then
    printf '[DRYRUN] %s\n' "${cmd[*]}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${chunk_tag}" "${count}" "${subjects_file}" "${analysis_root}" "${mask_root}" "DRYRUN" >> "${summary}"
  else
    jid=$("${cmd[@]}")
    printf '[INFO] submitted %s -> %s (%s subjects)\n' "${job_name}" "${jid}" "${count}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${chunk_tag}" "${count}" "${subjects_file}" "${analysis_root}" "${mask_root}" "${jid}" >> "${summary}"
  fi
  i=$(( i + chunk_size ))
  chunk=$(( chunk + 1 ))
done

echo "[INFO] Summary written to ${summary}"
