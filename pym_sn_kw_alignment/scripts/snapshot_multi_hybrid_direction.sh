#!/usr/bin/env bash
set -euo pipefail

ROOT=${SNHEART_ROOT:-/path/to/SNHEART}
TMP_DIR=${ROOT}/tmp
ALIGN_SCRIPT=${ROOT}/Scripts/check_sn_alignment_metrics.py

RUN_DIR=${1:-}
if [[ -z ${RUN_DIR} ]]; then
  RUN_DIR=$(ls -dt "${TMP_DIR}"/raw_kw_loop_* 2>/dev/null | head -n 1 || true)
fi
[[ -n ${RUN_DIR} ]] || { echo "[ERROR] Could not resolve RUN_DIR" >&2; exit 1; }
[[ -d ${RUN_DIR} ]] || { echo "[ERROR] RUN_DIR not found: ${RUN_DIR}" >&2; exit 1; }
[[ -f ${ALIGN_SCRIPT} ]] || { echo "[ERROR] Missing align checker: ${ALIGN_SCRIPT}" >&2; exit 1; }

PYTHON_BIN=${PYTHON_BIN:-python3}
if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
  PYTHON_BIN=python
fi
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || { echo "[ERROR] python missing" >&2; exit 1; }

multiopt_dir=${RUN_DIR}/multiopt_outputs
[[ -d ${multiopt_dir} ]] || { echo "[ERROR] No multiopt outputs found at ${multiopt_dir}" >&2; exit 1; }

ts=$(date -u +%Y%m%dT%H%M%SZ)
out_tsv=${RUN_DIR}/multiopt_direction_snapshot_${ts}.tsv
latest_tsv=${RUN_DIR}/multiopt_direction_snapshot_latest.tsv

printf "timestamp_utc\tapproach\tgroup\twritten\ttotal\tcloser\tfarther\tneutral\tavg_improvement_mm\n" > "${out_tsv}"

compare_counts() {
  local pre_tsv=$1
  local cur_tsv=$2
  "${PYTHON_BIN}" - "${pre_tsv}" "${cur_tsv}" <<'PY'
import csv
import math
import sys

pre_path, cur_path = sys.argv[1:3]
pre = {}
with open(pre_path, newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        pre[r["subject"]] = r

closer = farther = neutral = 0
improvements = []
with open(cur_path, newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        s = r.get("subject", "")
        p = pre.get(s)
        if not p:
            continue
        try:
            pd = float(p["delta_mag_mm"])
            cd = float(r["delta_mag_mm"])
        except Exception:
            continue
        if math.isnan(pd) or math.isnan(cd):
            continue
        imp = pd - cd
        improvements.append(imp)
        if imp > 1e-6:
            closer += 1
        elif imp < -1e-6:
            farther += 1
        else:
            neutral += 1

avg = sum(improvements) / len(improvements) if improvements else 0.0
print(f"{closer}\t{farther}\t{neutral}\t{avg:.6f}")
PY
}

for approach_dir in "${multiopt_dir}"/*; do
  [[ -d ${approach_dir} ]] || continue
  approach=$(basename "${approach_dir}")
  baseline_dir=${approach_dir}/probabilistic_masks_v10/base
  [[ -d ${baseline_dir} ]] || continue

  for group in highhat crash; do
    subjects_file=${approach_dir}/${group}_subjects.txt
    manual_dir=${RUN_DIR}/${group}/manual_raw_kw
    pre_tsv=${RUN_DIR}/${group}/${group}_rawkw_vs_v10_pre.tsv
    [[ -f ${subjects_file} ]] || continue
    [[ -d ${manual_dir} ]] || continue
    [[ -f ${pre_tsv} ]] || continue

    written_file=${approach_dir}/${group}_written_now.txt
    cur_tsv=${approach_dir}/${group}_current_now.tsv
    : > "${written_file}"
    while read -r sub; do
      [[ -z ${sub} ]] && continue
      if [[ -f ${baseline_dir}/sn_groupmask_in_${sub}.nii.gz || -f ${baseline_dir}/sn_groupmask_in_${sub}.nii ]]; then
        echo "${sub}" >> "${written_file}"
      fi
    done < "${subjects_file}"

    total=$(wc -l < "${subjects_file}")
    written=$(wc -l < "${written_file}")
    closer=0
    farther=0
    neutral=0
    avg="0.000000"
    if (( written > 0 )); then
      "${PYTHON_BIN}" "${ALIGN_SCRIPT}" \
        --subjects-file "${written_file}" \
        --manual-dir "${manual_dir}" \
        --baseline-dir "${baseline_dir}" \
        --output "${cur_tsv}" >/dev/null 2>&1 || true
      read -r closer farther neutral avg < <(compare_counts "${pre_tsv}" "${cur_tsv}")
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${approach}" "${group}" "${written}" "${total}" \
      "${closer}" "${farther}" "${neutral}" "${avg}" >> "${out_tsv}"
  done
done

cp -f "${out_tsv}" "${latest_tsv}"
echo "[DONE] snapshot=${out_tsv}"
echo "[DONE] latest=${latest_tsv}"
