#!/bin/bash -l
#SBATCH -p cloud
#SBATCH -c 4
#SBATCH -o sbatch-%j.out
#SBATCH -e sbatch-%j.err

set -euo pipefail

if command -v module >/dev/null 2>&1; then
  module load afni >/dev/null 2>&1 || true
fi
if ! command -v 3dAllineate >/dev/null 2>&1; then
  echo "[ERROR] AFNI/3dAllineate not available on PATH" >&2
  exit 1
fi

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
cd "$repo_root"

subjects_file=${AXIS_GAIN_SUBJECTS_FILE:?AXIS_GAIN_SUBJECTS_FILE required}
manual_dir=${AXIS_GAIN_MANUAL_DIR:?AXIS_GAIN_MANUAL_DIR required}
input_dir=${AXIS_GAIN_INPUT_DIR:?AXIS_GAIN_INPUT_DIR required}
out_root=${AXIS_GAIN_OUTPUT_ROOT:?AXIS_GAIN_OUTPUT_ROOT required}
gain_list=${AXIS_GAIN_GAINS:-"-0.5 -0.25 0 0.25 0.5 0.75 1.0"}

mkdir -p "${out_root}/pre" "${out_root}/best" "${out_root}/tmp" "${out_root}/metrics"

calc_delta() {
  python3 - "$1" "$2" <<'PY'
import sys
from pathlib import Path
import nibabel as nib
import numpy as np


def load_mask(path):
    img = nib.load(str(path))
    data = img.get_fdata() > 0
    coords = np.argwhere(data)
    if coords.size == 0:
        return None
    com = coords.mean(axis=0)
    return nib.affines.apply_affine(img.affine, com)


manual = load_mask(Path(sys.argv[1]))
base = load_mask(Path(sys.argv[2]))
if manual is None or base is None:
    print("nan nan nan")
    sys.exit(0)
delta = manual - base
print(f"{delta[0]:.6f} {delta[1]:.6f} {delta[2]:.6f}")
PY
}

stage_manifest="${out_root}/stage_manifest.tsv"
python3 - "${stage_manifest}" ${gain_list} <<'PY'
import itertools
import sys
from pathlib import Path


def stage_token(value: float) -> str:
    scaled = int(round(abs(value) * 100))
    if abs(value) < 1e-9:
        prefix = "z"
    elif value < 0:
        prefix = "m"
    else:
        prefix = "p"
    return f"{prefix}{scaled:03d}"


manifest_path = Path(sys.argv[1])
gains = [float(x) for x in sys.argv[2:]]
with manifest_path.open("w", encoding="utf-8") as handle:
    handle.write("stage\tgx\tgy\tgz\n")
    for gx, gy, gz in itertools.product(gains, repeat=3):
        if gx == 0.0 and gy == 0.0 and gz == 0.0:
            continue
        stage = f"x{stage_token(gx)}_y{stage_token(gy)}_z{stage_token(gz)}"
        handle.write(f"{stage}\t{gx:.2f}\t{gy:.2f}\t{gz:.2f}\n")
PY

while read -r sid; do
  [[ -z ${sid} ]] && continue
  sid=${sid##* }
  manual_mask="${manual_dir}/sn_groupmask_in_${sid}.nii.gz"
  [[ -f ${manual_mask} ]] || manual_mask="${manual_dir}/sn_groupmask_in_${sid}.nii"
  pre_mask="${input_dir}/sn_groupmask_in_${sid}.nii.gz"
  [[ -f ${pre_mask} ]] || pre_mask="${input_dir}/sn_groupmask_in_${sid}.nii"
  [[ -f ${manual_mask} && -f ${pre_mask} ]] || { echo "[WARN] ${sid}: missing manual or input mask"; continue; }

  pre_out="${out_root}/pre/$(basename "${pre_mask}")"
  ln -sf "${pre_mask}" "${pre_out}"

  read -r dx dy dz <<<"$(calc_delta "${manual_mask}" "${pre_mask}")"
  if [[ ${dx} == "nan" ]]; then
    echo "[WARN] ${sid}: failed COM delta"
    continue
  fi

  while IFS=$'\t' read -r stage gx gy gz; do
    [[ ${stage} == "stage" ]] && continue
    mkdir -p "${out_root}/${stage}"
    read -r sdx sdy sdz <<<"$(python3 - "${dx}" "${dy}" "${dz}" "${gx}" "${gy}" "${gz}" <<'PY'
import sys
dx, dy, dz, gx, gy, gz = map(float, sys.argv[1:7])
print(f"{dx*gx:.6f} {dy*gy:.6f} {dz*gz:.6f}")
PY
)"
    trans_mat="${out_root}/tmp/translate_${stage}_${sid}.1D"
    printf "1 0 0 %s 0 1 0 %s 0 0 1 %s\n" "${sdx}" "${sdy}" "${sdz}" > "${trans_mat}"
    trans_out="${out_root}/${stage}/sn_groupmask_in_${sid}.nii.gz"
    3dAllineate -overwrite -base "${pre_mask}" -input "${pre_mask}" \
      -1Dmatrix_apply "${trans_mat}" -prefix "${trans_out}" -final NN -float >/dev/null 2>&1
  done < "${stage_manifest}"
done < "${subjects_file}"

python3 Scripts/check_sn_alignment_metrics.py --subjects-file "${subjects_file}" \
  --manual-dir "${manual_dir}" --baseline-dir "${out_root}/pre" \
  --output "${out_root}/metrics/pre.tsv"

while IFS=$'\t' read -r stage gx gy gz; do
  [[ ${stage} == "stage" ]] && continue
  python3 Scripts/check_sn_alignment_metrics.py --subjects-file "${subjects_file}" \
    --manual-dir "${manual_dir}" --baseline-dir "${out_root}/${stage}" \
    --output "${out_root}/metrics/${stage}.tsv"
done < "${stage_manifest}"

python3 - "${out_root}" <<'PY'
import csv
import math
import sys
from pathlib import Path


root = Path(sys.argv[1])
stage_manifest = root / "stage_manifest.tsv"
metrics_dir = root / "metrics"

stages = []
with stage_manifest.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        stages.append(row)

tables = {}
for stage in ["pre", *[row["stage"] for row in stages]]:
    path = metrics_dir / f"{stage}.tsv"
    if not path.exists():
        continue
    with path.open(newline="", encoding="utf-8") as handle:
        tables[stage] = {row["subject"]: row for row in csv.DictReader(handle, delimiter="\t")}

subjects = sorted({subject for table in tables.values() for subject in table}, key=int)
best_dir = root / "best"
best_dir.mkdir(exist_ok=True)
summary_path = root / "best_axis_gain_summary.tsv"
counts_path = root / "best_axis_gain_stage_counts.tsv"
summary_rows = []
counts = {}

for subject in subjects:
    pre_row = tables.get("pre", {}).get(subject)
    if not pre_row:
        continue
    try:
        pre_delta = float(pre_row["delta_mag_mm"])
    except Exception:
        continue

    best_stage = "pre"
    best_delta = pre_delta
    best_gx = best_gy = best_gz = 0.0
    best_source_stage = "pre"
    best_source_delta = pre_delta

    for row in stages:
        stage = row["stage"]
        metric_row = tables.get(stage, {}).get(subject)
        if not metric_row:
            continue
        try:
            delta = float(metric_row["delta_mag_mm"])
        except Exception:
            continue
        if delta < best_source_delta:
            best_source_stage = stage
            best_source_delta = delta
            best_gx = float(row["gx"])
            best_gy = float(row["gy"])
            best_gz = float(row["gz"])

    if best_source_delta < pre_delta:
        best_stage = best_source_stage
        best_delta = best_source_delta

    improvement = pre_delta - best_delta
    counts[best_stage] = counts.get(best_stage, 0) + 1

    src = root / best_stage / f"sn_groupmask_in_{subject}.nii.gz"
    if not src.exists():
        src = root / best_stage / f"sn_groupmask_in_{subject}.nii"
    if src.exists():
        dst = best_dir / src.name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        dst.symlink_to(src)

    summary_rows.append(
        [
            subject,
            best_stage,
            f"{best_delta:.6f}",
            f"{pre_delta:.6f}",
            f"{improvement:.6f}",
            best_source_stage,
            f"{best_source_delta:.6f}",
            f"{best_gx:.2f}",
            f"{best_gy:.2f}",
            f"{best_gz:.2f}",
            "1" if best_stage != "pre" else "0",
        ]
    )

with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(
        [
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
        ]
    )
    writer.writerows(summary_rows)

with counts_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["stage", "count"])
    for stage, count in sorted(counts.items()):
        writer.writerow([stage, count])

print(f"WROTE {summary_path}")
print(f"WROTE {counts_path}")
print(
    "SUMMARY",
    f"subjects={len(summary_rows)}",
    f"improved={sum(1 for row in summary_rows if row[-1] == '1')}",
    f"stayed_pre={sum(1 for row in summary_rows if row[1] == 'pre')}",
)
PY
