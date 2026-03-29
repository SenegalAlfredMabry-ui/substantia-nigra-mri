#!/usr/bin/env bash
set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}

refine_root=${AXIS_OUTLIER_REFINE_ROOT:-$(
  find "${repo_root}/tmp" -maxdepth 1 -mindepth 1 -type d -name 'native_sn_axis_gain_refine_all108_20*' \
    | sort -r \
    | head -n 1
)}
axis_root=${AXIS_OUTLIER_AXIS_ROOT:-$(
  find "${repo_root}/tmp" -maxdepth 1 -mindepth 1 -type d -name 'native_sn_axis_gain_20*' ! -name 'native_sn_axis_gain_refine_*' \
    | sort -r \
    | head -n 1
)}

[[ -n ${refine_root} && -d ${refine_root} ]] || { echo "[ERROR] Could not resolve refine-all108 root" >&2; exit 1; }
[[ -n ${axis_root} && -d ${axis_root} ]] || { echo "[ERROR] Could not resolve axis root" >&2; exit 1; }

out_root=${AXIS_REFINE_ROOT:-${repo_root}/tmp/native_sn_axis_gain_refine_outliers_$(date -u +%Y%m%dT%H%M%SZ)}
launcher=${repo_root}/Scripts/launch_native_sn_axis_gain_refine_chunks.sh

outlier_subjects_default="222 210 142 128 221 526 551 212 563 227 140"
outlier_subjects=${AXIS_OUTLIER_SUBJECTS:-${outlier_subjects_default}}
outlier_subjects_file=${AXIS_OUTLIER_SUBJECTS_FILE:-}

mkdir -p "${out_root}"

converted_summary="${out_root}/current_best_summary_outliers.tsv"
current_dir="${out_root}/current_best_masks"
subjects_file="${out_root}/subjects_outliers.txt"
selection_audit="${out_root}/selection_audit.tsv"

python3 - "${refine_root}" "${converted_summary}" "${current_dir}" "${subjects_file}" "${selection_audit}" "${outlier_subjects_file}" "${outlier_subjects}" <<'PY'
import csv
import re
import sys
from pathlib import Path


refine_root = Path(sys.argv[1])
summary_out = Path(sys.argv[2])
current_dir = Path(sys.argv[3])
subjects_out = Path(sys.argv[4])
audit_out = Path(sys.argv[5])
subjects_file = Path(sys.argv[6]) if sys.argv[6] else None
subjects_text = sys.argv[7]

current_dir.mkdir(parents=True, exist_ok=True)


def read_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


rows_by_subject = {}
for path in sorted(refine_root.glob("chunks/c*/best_axis_refine_summary.tsv")):
    for row in read_rows(path):
        rows_by_subject[row["subject"]] = row

mask_by_subject = {}
for path in sorted(refine_root.glob("chunks/c*/best/sn_groupmask_in_*.nii.gz")):
    subject = path.name.replace("sn_groupmask_in_", "").replace(".nii.gz", "")
    mask_by_subject[subject] = path
for path in sorted(refine_root.glob("chunks/c*/best/sn_groupmask_in_*.nii")):
    subject = path.name.replace("sn_groupmask_in_", "").replace(".nii", "")
    mask_by_subject.setdefault(subject, path)

requested = []
if subjects_file and subjects_file.exists():
    requested = [line.strip().split()[0] for line in subjects_file.read_text(encoding="utf-8").splitlines() if line.strip()]
else:
    requested = [tok for tok in re.split(r"[\s,]+", subjects_text.strip()) if tok]

requested = [str(int(tok)) for tok in requested]

selected = []
audit_rows = []
for sid in requested:
    has_row = sid in rows_by_subject
    has_mask = sid in mask_by_subject
    status = "selected" if has_row and has_mask else "missing_summary_or_mask"
    audit_rows.append([sid, "1" if has_row else "0", "1" if has_mask else "0", status])
    if has_row and has_mask:
        selected.append(sid)

if not selected:
    raise SystemExit("No selected outlier subjects had both summary rows and best masks.")

fieldnames = [
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

with summary_out.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for sid in selected:
        src = rows_by_subject[sid]
        accepted_delta = float(src["accepted_delta_mm"])
        writer.writerow(
            {
                "subject": sid,
                "accepted_stage": src["accepted_stage"],
                "accepted_delta_mm": f"{accepted_delta:.6f}",
                "pre_delta_mm": f"{accepted_delta:.6f}",
                "improvement_mm": "0.000000",
                "best_candidate_stage": src["accepted_stage"],
                "best_candidate_delta_mm": f"{accepted_delta:.6f}",
                "best_gx": f"{float(src['accepted_gx']):.3f}",
                "best_gy": f"{float(src['accepted_gy']):.3f}",
                "best_gz": f"{float(src['accepted_gz']):.3f}",
                "accepted_improvement": src.get("accepted_improvement", "0"),
            }
        )

subjects_out.write_text("\n".join(selected) + "\n", encoding="utf-8")

for sid in selected:
    src = mask_by_subject[sid]
    dst = current_dir / src.name
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src)

with audit_out.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["subject", "has_refine_summary", "has_best_mask", "status"])
    writer.writerows(audit_rows)

print(summary_out)
print(subjects_out)
print(audit_out)
print(f"requested={len(requested)} selected={len(selected)}")
PY

# Expanded local windows around each subject's current accepted gains.
x_offsets=${AXIS_REFINE_X_OFFSETS:-"-0.75 -0.5 -0.375 -0.25 -0.125 0 0.125 0.25 0.375 0.5 0.75"}
y_offsets=${AXIS_REFINE_Y_OFFSETS:-"-0.75 -0.5 -0.375 -0.25 -0.125 0 0.125 0.25 0.375 0.5 0.75"}
z_offsets=${AXIS_REFINE_Z_OFFSETS:-"-0.75 -0.5 -0.375 -0.25 -0.125 0 0.125 0.25 0.375 0.5 0.75"}

x_min=${AXIS_REFINE_X_MIN:--1.5}
x_max=${AXIS_REFINE_X_MAX:-2.0}
y_min=${AXIS_REFINE_Y_MIN:--1.5}
y_max=${AXIS_REFINE_Y_MAX:-1.5}
z_min=${AXIS_REFINE_Z_MIN:--1.75}
z_max=${AXIS_REFINE_Z_MAX:-1.25}

chunk_count=${AXIS_REFINE_CHUNK_COUNT:-4}
enable_rigid=${AXIS_REFINE_ENABLE_RIGID:-1}
rigid_suffix=${AXIS_REFINE_RIGID_SUFFIX:-__rig}
rigid_maxrot_deg=${AXIS_REFINE_RIGID_MAXROT_DEG:-1.0}
rigid_maxshf_mm=${AXIS_REFINE_RIGID_MAXSHF_MM:-0.5}
rigid_cost=${AXIS_REFINE_RIGID_COST:-mi}

AXIS_REFINE_ROOT="${out_root}" \
AXIS_REFINE_PREV_ROOT="${axis_root}" \
AXIS_REFINE_PREVIOUS_SUMMARY="${converted_summary}" \
AXIS_REFINE_CURRENT_DIR="${current_dir}" \
AXIS_REFINE_SUBJECTS_FILE="${subjects_file}" \
AXIS_REFINE_CHUNK_COUNT="${chunk_count}" \
AXIS_REFINE_X_OFFSETS="${x_offsets}" \
AXIS_REFINE_Y_OFFSETS="${y_offsets}" \
AXIS_REFINE_Z_OFFSETS="${z_offsets}" \
AXIS_REFINE_X_MIN="${x_min}" \
AXIS_REFINE_X_MAX="${x_max}" \
AXIS_REFINE_Y_MIN="${y_min}" \
AXIS_REFINE_Y_MAX="${y_max}" \
AXIS_REFINE_Z_MIN="${z_min}" \
AXIS_REFINE_Z_MAX="${z_max}" \
AXIS_REFINE_ENABLE_RIGID="${enable_rigid}" \
AXIS_REFINE_RIGID_SUFFIX="${rigid_suffix}" \
AXIS_REFINE_RIGID_MAXROT_DEG="${rigid_maxrot_deg}" \
AXIS_REFINE_RIGID_MAXSHF_MM="${rigid_maxshf_mm}" \
AXIS_REFINE_RIGID_COST="${rigid_cost}" \
bash "${launcher}"

