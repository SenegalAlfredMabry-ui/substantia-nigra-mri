#!/usr/bin/env bash
set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
if [[ -n ${AXIS_REFINE_PREV_ROOT:-} ]]; then
  prev_root=${AXIS_REFINE_PREV_ROOT}
else
  prev_root=$(
    find "${repo_root}/tmp" -maxdepth 1 -mindepth 1 -type d -name 'native_sn_axis_gain_20*' ! -name 'native_sn_axis_gain_refine_*' \
      | sort -r \
      | head -n 1
  )
fi
[[ -n ${prev_root} && -d ${prev_root} ]] || { echo "[ERROR] Could not resolve previous axis-gain root" >&2; exit 1; }

out_root=${AXIS_REFINE_ROOT:-${repo_root}/tmp/native_sn_axis_gain_refine_$(date -u +%Y%m%dT%H%M%SZ)}
chunk_count=${AXIS_REFINE_CHUNK_COUNT:-4}
min_delta=${AXIS_REFINE_MIN_DELTA_MM:-4.0}
runner=${repo_root}/Scripts/run_native_sn_axis_gain_refine_chunk.sh
subjects_file=${AXIS_REFINE_SUBJECTS_FILE:-}
previous_summary_override=${AXIS_REFINE_PREVIOUS_SUMMARY:-}
current_dir=${AXIS_REFINE_CURRENT_DIR:-}

x_offsets=${AXIS_REFINE_X_OFFSETS:-"-0.25 -0.125 0 0.125 0.25 0.375 0.5"}
y_offsets=${AXIS_REFINE_Y_OFFSETS:-"-0.25 -0.125 0 0.125 0.25"}
z_offsets=${AXIS_REFINE_Z_OFFSETS:-"-0.5 -0.375 -0.25 -0.125 0 0.125 0.25"}
x_min=${AXIS_REFINE_X_MIN:--0.5}
x_max=${AXIS_REFINE_X_MAX:-1.5}
y_min=${AXIS_REFINE_Y_MIN:--1.0}
y_max=${AXIS_REFINE_Y_MAX:-1.25}
z_min=${AXIS_REFINE_Z_MIN:--1.25}
z_max=${AXIS_REFINE_Z_MAX:-0.5}
enable_rigid=${AXIS_REFINE_ENABLE_RIGID:-0}
rigid_suffix=${AXIS_REFINE_RIGID_SUFFIX:-__rig}
rigid_maxrot_deg=${AXIS_REFINE_RIGID_MAXROT_DEG:-1.5}
rigid_maxshf_mm=${AXIS_REFINE_RIGID_MAXSHF_MM:-0.75}
rigid_cost=${AXIS_REFINE_RIGID_COST:-mi}

mkdir -p "${out_root}"

if [[ -n ${previous_summary_override} ]]; then
  merged_summary=${previous_summary_override}
else
  merged_summary="${out_root}/previous_axis_gain_summary_merged.tsv"
  python3 - "${prev_root}" "${merged_summary}" <<'PY'
import csv
import sys
from pathlib import Path

prev_root = Path(sys.argv[1])
out_path = Path(sys.argv[2])
rows = []
for path in sorted((prev_root / "chunks").glob("c*/best_axis_gain_summary.tsv")):
    with path.open(newline="", encoding="utf-8") as handle:
        rows.extend(csv.DictReader(handle, delimiter="\t"))
rows.sort(key=lambda row: int(row["subject"]))
with out_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "subject",
            "accepted_stage",
            "accepted_delta_mm",
            "pre_delta_mm",
            "improvement_mm",
            "best_candidate_stage",
            "best_candidate_delta_mm",
            "best_gx",
            "best_gy",
            "best_gz",
            "accepted_improvement",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
print(out_path)
print(f"count={len(rows)}")
PY
fi
[[ -f ${merged_summary} ]] || { echo "[ERROR] Previous summary not found: ${merged_summary}" >&2; exit 1; }

if [[ -z ${subjects_file} ]]; then
  subjects_file="${out_root}/subjects_ge_${min_delta}mm.txt"
  python3 - "${merged_summary}" "${subjects_file}" "${min_delta}" <<'PY'
import csv
import sys
from pathlib import Path

summary_path = Path(sys.argv[1])
subjects_path = Path(sys.argv[2])
min_delta = float(sys.argv[3])
subjects = []
with summary_path.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        try:
            delta = float(row["accepted_delta_mm"])
        except Exception:
            continue
        if delta >= min_delta:
            subjects.append(row["subject"])
subjects_path.write_text("\n".join(subjects) + ("\n" if subjects else ""), encoding="utf-8")
print(subjects_path)
print(f"count={len(subjects)}")
PY
fi

mapfile -t subjects < "${subjects_file}"
[[ ${#subjects[@]} -gt 0 ]] || { echo "[ERROR] No subjects selected for axis refinement" >&2; exit 1; }

manual_dir="${prev_root}/inputs/manual_masks"
pre_dir="${prev_root}/inputs/staged_masks"
[[ -d ${manual_dir} && -d ${pre_dir} ]] || { echo "[ERROR] Missing manual/pre dirs under previous axis root" >&2; exit 1; }

summary="${out_root}/submission_summary.tsv"
printf "chunk\tjobid\tcount\tsubjects_file\toutput_root\tmin_delta_mm\tx_offsets\ty_offsets\tz_offsets\tenable_rigid\trigid_maxrot_deg\trigid_maxshf_mm\trigid_cost\n" > "${summary}"

if (( chunk_count > ${#subjects[@]} )); then
  chunk_count=${#subjects[@]}
fi
chunk_size=$(( (${#subjects[@]} + chunk_count - 1) / chunk_count ))

i=0
chunk=1
while (( i < ${#subjects[@]} )); do
  tag=$(printf "c%02d" "${chunk}")
  chunk_dir="${out_root}/chunks/${tag}"
  mkdir -p "${chunk_dir}"
  chunk_subjects_file="${chunk_dir}/subjects.txt"
  printf "%s\n" "${subjects[@]:i:chunk_size}" > "${chunk_subjects_file}"
  job_name="axis_refine_${tag}"
  jid=$(sbatch --parsable --job-name="${job_name}" \
    --output="${chunk_dir}/${job_name}.out" --error="${chunk_dir}/${job_name}.err" \
    --export=ALL,AXIS_REFINE_SUBJECTS_FILE="${chunk_subjects_file}",AXIS_REFINE_MANUAL_DIR="${manual_dir}",AXIS_REFINE_PRE_DIR="${pre_dir}",AXIS_REFINE_PREV_ROOT="${prev_root}",AXIS_REFINE_PREVIOUS_SUMMARY="${merged_summary}",AXIS_REFINE_CURRENT_DIR="${current_dir}",AXIS_REFINE_OUTPUT_ROOT="${chunk_dir}",AXIS_REFINE_X_OFFSETS="${x_offsets}",AXIS_REFINE_Y_OFFSETS="${y_offsets}",AXIS_REFINE_Z_OFFSETS="${z_offsets}",AXIS_REFINE_X_MIN="${x_min}",AXIS_REFINE_X_MAX="${x_max}",AXIS_REFINE_Y_MIN="${y_min}",AXIS_REFINE_Y_MAX="${y_max}",AXIS_REFINE_Z_MIN="${z_min}",AXIS_REFINE_Z_MAX="${z_max}",AXIS_REFINE_ENABLE_RIGID="${enable_rigid}",AXIS_REFINE_RIGID_SUFFIX="${rigid_suffix}",AXIS_REFINE_RIGID_MAXROT_DEG="${rigid_maxrot_deg}",AXIS_REFINE_RIGID_MAXSHF_MM="${rigid_maxshf_mm}",AXIS_REFINE_RIGID_COST="${rigid_cost}" \
    "${runner}")
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${tag}" "${jid}" "$(wc -l < "${chunk_subjects_file}")" "${chunk_subjects_file}" "${chunk_dir}" "${min_delta}" "${x_offsets}" "${y_offsets}" "${z_offsets}" "${enable_rigid}" "${rigid_maxrot_deg}" "${rigid_maxshf_mm}" "${rigid_cost}" >> "${summary}"
  i=$(( i + chunk_size ))
  chunk=$(( chunk + 1 ))
done

echo "[INFO] Submission summary written to ${summary}"
