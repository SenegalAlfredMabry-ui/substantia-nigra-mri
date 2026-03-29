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
cd "${repo_root}"

subjects_file=${AXIS_REFINE_SUBJECTS_FILE:?AXIS_REFINE_SUBJECTS_FILE required}
manual_dir=${AXIS_REFINE_MANUAL_DIR:?AXIS_REFINE_MANUAL_DIR required}
pre_dir=${AXIS_REFINE_PRE_DIR:?AXIS_REFINE_PRE_DIR required}
prev_root=${AXIS_REFINE_PREV_ROOT:?AXIS_REFINE_PREV_ROOT required}
prev_summary=${AXIS_REFINE_PREVIOUS_SUMMARY:?AXIS_REFINE_PREVIOUS_SUMMARY required}
current_dir=${AXIS_REFINE_CURRENT_DIR:-}
out_root=${AXIS_REFINE_OUTPUT_ROOT:?AXIS_REFINE_OUTPUT_ROOT required}

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

mkdir -p "${out_root}/pre" "${out_root}/current" "${out_root}/best" "${out_root}/tmp" "${out_root}/metrics"

calc_delta() {
  python3 - "$1" "$2" <<'PY'
import sys
from pathlib import Path
import nibabel as nib
import numpy as np


def load_mask(path: Path):
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
all_stage_manifest="${out_root}/all_stage_manifest.tsv"
subject_stage_manifest="${out_root}/subject_stage_manifest.tsv"
python3 - "${subjects_file}" "${prev_summary}" "${stage_manifest}" "${subject_stage_manifest}" \
  "${x_offsets}" "${y_offsets}" "${z_offsets}" \
  "${x_min}" "${x_max}" "${y_min}" "${y_max}" "${z_min}" "${z_max}" <<'PY'
import csv
import itertools
import math
import sys
from pathlib import Path


def read_subjects(path: Path) -> list[str]:
    subjects = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            token = line.strip()
            if token:
                subjects.append(token.split()[0])
    return subjects


def parse_offsets(text: str) -> list[float]:
    return [float(tok) for tok in text.split()]


def clamp(value: float, lower: float, upper: float) -> float:
    return min(max(value, lower), upper)


def stage_token(value: float) -> str:
    scaled = int(round(abs(value) * 1000))
    if abs(value) < 1e-9:
        prefix = "z"
    elif value < 0:
        prefix = "m"
    else:
        prefix = "p"
    return f"{prefix}{scaled:04d}"


def stage_name(gx: float, gy: float, gz: float) -> str:
    return f"x{stage_token(gx)}_y{stage_token(gy)}_z{stage_token(gz)}"


subjects_path = Path(sys.argv[1])
summary_path = Path(sys.argv[2])
stage_manifest_path = Path(sys.argv[3])
subject_stage_manifest_path = Path(sys.argv[4])
x_offsets = parse_offsets(sys.argv[5])
y_offsets = parse_offsets(sys.argv[6])
z_offsets = parse_offsets(sys.argv[7])
x_min = float(sys.argv[8])
x_max = float(sys.argv[9])
y_min = float(sys.argv[10])
y_max = float(sys.argv[11])
z_min = float(sys.argv[12])
z_max = float(sys.argv[13])

summary_rows = {}
with summary_path.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        summary_rows[row["subject"]] = row

stage_rows: dict[str, tuple[float, float, float]] = {}
subjects = read_subjects(subjects_path)

with subject_stage_manifest_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(
        [
            "subject",
            "stage",
            "gx",
            "gy",
            "gz",
            "current_stage",
            "current_delta_mm",
            "current_gx",
            "current_gy",
            "current_gz",
        ]
    )

    for subject in subjects:
        row = summary_rows.get(subject)
        if row is None:
            continue
        current_stage = row["accepted_stage"]
        current_delta = float(row["accepted_delta_mm"])
        center_x = float(row.get("best_gx") or 0.0)
        center_y = float(row.get("best_gy") or 0.0)
        center_z = float(row.get("best_gz") or 0.0)

        xs = sorted({round(clamp(center_x + off, x_min, x_max), 3) for off in x_offsets})
        ys = sorted({round(clamp(center_y + off, y_min, y_max), 3) for off in y_offsets})
        zs = sorted({round(clamp(center_z + off, z_min, z_max), 3) for off in z_offsets})

        for gx, gy, gz in itertools.product(xs, ys, zs):
            if (
                math.isclose(gx, center_x, abs_tol=1e-9)
                and math.isclose(gy, center_y, abs_tol=1e-9)
                and math.isclose(gz, center_z, abs_tol=1e-9)
            ):
                continue
            stage = stage_name(gx, gy, gz)
            stage_rows[stage] = (gx, gy, gz)
            writer.writerow(
                [
                    subject,
                    stage,
                    f"{gx:.3f}",
                    f"{gy:.3f}",
                    f"{gz:.3f}",
                    current_stage,
                    f"{current_delta:.6f}",
                    f"{center_x:.3f}",
                    f"{center_y:.3f}",
                    f"{center_z:.3f}",
                ]
            )

with stage_manifest_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["stage", "gx", "gy", "gz"])
    for stage in sorted(stage_rows):
        gx, gy, gz = stage_rows[stage]
        writer.writerow([stage, f"{gx:.3f}", f"{gy:.3f}", f"{gz:.3f}"])
PY

cp "${stage_manifest}" "${all_stage_manifest}"
if [[ ${enable_rigid} == "1" ]]; then
  awk -F'\t' -v OFS='\t' -v suffix="${rigid_suffix}" 'NR>1 {print $1 suffix, $2, $3, $4}' "${stage_manifest}" >> "${all_stage_manifest}"
fi

find_prev_mask() {
  local sid=$1
  local stage=$2
  local path=
  for ext in ".nii.gz" ".nii"; do
    for candidate in "${prev_root}"/chunks/*/"${stage}"/"sn_groupmask_in_${sid}${ext}"; do
      if [[ -e ${candidate} ]]; then
        path=${candidate}
        break 2
      fi
    done
  done
  if [[ -n ${path} ]]; then
    printf "%s\n" "${path}"
    return 0
  fi
  return 1
}

while read -r sid; do
  [[ -z ${sid} ]] && continue
  sid=${sid##* }
  manual_mask="${manual_dir}/sn_groupmask_in_${sid}.nii.gz"
  [[ -f ${manual_mask} ]] || manual_mask="${manual_dir}/sn_groupmask_in_${sid}.nii"
  pre_mask="${pre_dir}/sn_groupmask_in_${sid}.nii.gz"
  [[ -f ${pre_mask} ]] || pre_mask="${pre_dir}/sn_groupmask_in_${sid}.nii"
  [[ -f ${manual_mask} && -f ${pre_mask} ]] || { echo "[WARN] ${sid}: missing manual or pre mask"; continue; }

  pre_out="${out_root}/pre/$(basename "${pre_mask}")"
  ln -sf "${pre_mask}" "${pre_out}"

  current_stage=$(awk -F'\t' -v sid="${sid}" 'NR>1 && $1==sid {print $6; exit}' "${subject_stage_manifest}")
  [[ -n ${current_stage} ]] || { echo "[WARN] ${sid}: missing current stage"; continue; }
  current_mask=
  if [[ -n ${current_dir} ]]; then
    for ext in ".nii.gz" ".nii"; do
      candidate="${current_dir}/sn_groupmask_in_${sid}${ext}"
      if [[ -e ${candidate} ]]; then
        current_mask=${candidate}
        break
      fi
    done
  fi
  if [[ -z ${current_mask} ]]; then
    current_mask=$(find_prev_mask "${sid}" "${current_stage}" || true)
  fi
  [[ -n ${current_mask} ]] || { echo "[WARN] ${sid}: could not locate current mask for stage ${current_stage}"; continue; }
  ln -sf "${current_mask}" "${out_root}/current/$(basename "${current_mask}")"

  read -r dx dy dz <<<"$(calc_delta "${manual_mask}" "${pre_mask}")"
  if [[ ${dx} == "nan" ]]; then
    echo "[WARN] ${sid}: failed COM delta"
    continue
  fi

  awk -F'\t' -v sid="${sid}" 'NR>1 && $1==sid' "${subject_stage_manifest}" | \
  while IFS=$'\t' read -r subject stage gx gy gz _; do
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

    if [[ ${enable_rigid} == "1" ]]; then
      rigid_stage="${stage}${rigid_suffix}"
      rigid_out="${out_root}/${rigid_stage}/sn_groupmask_in_${sid}.nii.gz"
      mkdir -p "${out_root}/${rigid_stage}"
      3dAllineate -overwrite -base "${manual_mask}" -source "${trans_out}" \
        -prefix "${rigid_out}" -final NN -float -warp shift_rotate \
        -cost "${rigid_cost}" -cmass \
        -maxrot "${rigid_maxrot_deg}" -maxshf "${rigid_maxshf_mm}" >/dev/null 2>&1 || {
          rm -f "${rigid_out}"
        }
    fi
  done
done < "${subjects_file}"

python3 Scripts/check_sn_alignment_metrics.py --subjects-file "${subjects_file}" \
  --manual-dir "${manual_dir}" --baseline-dir "${out_root}/current" \
  --output "${out_root}/metrics/current.tsv"

while IFS=$'\t' read -r stage gx gy gz; do
  [[ ${stage} == "stage" ]] && continue
  python3 Scripts/check_sn_alignment_metrics.py --subjects-file "${subjects_file}" \
    --manual-dir "${manual_dir}" --baseline-dir "${out_root}/${stage}" \
    --output "${out_root}/metrics/${stage}.tsv"
done < "${all_stage_manifest}"

python3 - "${out_root}" "${subject_stage_manifest}" "${all_stage_manifest}" <<'PY'
import csv
import sys
from pathlib import Path


root = Path(sys.argv[1])
subject_stage_manifest = Path(sys.argv[2])
all_stage_manifest = Path(sys.argv[3])
metrics_dir = root / "metrics"
best_dir = root / "best"
best_dir.mkdir(exist_ok=True)

subject_info = {}
all_stages = []
stage_gains = {}
with subject_stage_manifest.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        subject_info.setdefault(
            row["subject"],
            {
                "current_stage": row["current_stage"],
                "current_gx": row["current_gx"],
                "current_gy": row["current_gy"],
                "current_gz": row["current_gz"],
            },
        )

with all_stage_manifest.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        stage = row["stage"]
        if stage in stage_gains:
            continue
        all_stages.append(stage)
        stage_gains[stage] = (
            float(row["gx"]),
            float(row["gy"]),
            float(row["gz"]),
        )


def load_table(path: Path) -> dict[str, dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["subject"]: row for row in csv.DictReader(handle, delimiter="\t")}


tables = {"current": load_table(metrics_dir / "current.tsv")}
for stage in all_stages:
    path = metrics_dir / f"{stage}.tsv"
    if path.exists():
        tables[stage] = load_table(path)


summary_rows = []
counts = {}
for subject in sorted(subject_info, key=int):
    current_row = tables["current"].get(subject)
    if current_row is None:
        continue
    try:
        current_delta = float(current_row["delta_mag_mm"])
    except Exception:
        continue

    accepted_stage = subject_info[subject]["current_stage"]
    accepted_delta = current_delta
    accepted_gx = float(subject_info[subject]["current_gx"])
    accepted_gy = float(subject_info[subject]["current_gy"])
    accepted_gz = float(subject_info[subject]["current_gz"])
    best_candidate_stage = accepted_stage
    best_candidate_delta = current_delta
    best_candidate_gx = accepted_gx
    best_candidate_gy = accepted_gy
    best_candidate_gz = accepted_gz

    for stage in all_stages:
        row = tables.get(stage, {}).get(subject)
        if row is None:
            continue
        try:
            delta = float(row["delta_mag_mm"])
        except Exception:
            continue
        if delta < best_candidate_delta:
            best_candidate_delta = delta
            best_candidate_stage = stage
            best_candidate_gx, best_candidate_gy, best_candidate_gz = stage_gains.get(
                stage, (accepted_gx, accepted_gy, accepted_gz)
            )

    if best_candidate_delta < current_delta:
        accepted_stage = best_candidate_stage
        accepted_delta = best_candidate_delta
        accepted_gx = best_candidate_gx
        accepted_gy = best_candidate_gy
        accepted_gz = best_candidate_gz

    counts[accepted_stage] = counts.get(accepted_stage, 0) + 1

    if accepted_stage == subject_info[subject]["current_stage"]:
        src = root / "current" / f"sn_groupmask_in_{subject}.nii.gz"
        if not src.exists():
            src = root / "current" / f"sn_groupmask_in_{subject}.nii"
    else:
        src = root / accepted_stage / f"sn_groupmask_in_{subject}.nii.gz"
        if not src.exists():
            src = root / accepted_stage / f"sn_groupmask_in_{subject}.nii"
    if src.exists():
        dst = best_dir / src.name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        dst.symlink_to(src)

    summary_rows.append(
        [
            subject,
            accepted_stage,
            f"{accepted_delta:.6f}",
            f"{current_delta:.6f}",
            f"{(current_delta - accepted_delta):.6f}",
            subject_info[subject]["current_stage"],
            f"{current_delta:.6f}",
            best_candidate_stage,
            f"{best_candidate_delta:.6f}",
            f"{accepted_gx:.3f}",
            f"{accepted_gy:.3f}",
            f"{accepted_gz:.3f}",
            "1" if accepted_stage != subject_info[subject]["current_stage"] else "0",
        ]
    )

summary_path = root / "best_axis_refine_summary.tsv"
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(
        [
            "subject",
            "accepted_stage",
            "accepted_delta_mm",
            "current_delta_mm",
            "improvement_mm",
            "previous_stage",
            "previous_delta_mm",
            "best_candidate_stage",
            "best_candidate_delta_mm",
            "accepted_gx",
            "accepted_gy",
            "accepted_gz",
            "accepted_improvement",
        ]
    )
    writer.writerows(summary_rows)

counts_path = root / "best_axis_refine_stage_counts.tsv"
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
    f"stayed_current={sum(1 for row in summary_rows if row[-1] == '0')}",
)
PY
