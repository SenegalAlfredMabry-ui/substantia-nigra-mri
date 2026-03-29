#!/usr/bin/env bash
set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}

translation_root=${AXIS_REFINE_ALL108_TRANSLATION_ROOT:-$(
  find "${repo_root}/tmp" -maxdepth 1 -mindepth 1 -type d -name 'native_sn_translation_gain_20*' \
    | sort -r \
    | head -n 1
)}
axis_root=${AXIS_REFINE_ALL108_AXIS_ROOT:-$(
  find "${repo_root}/tmp" -maxdepth 1 -mindepth 1 -type d -name 'native_sn_axis_gain_20*' ! -name 'native_sn_axis_gain_refine_*' \
    | sort -r \
    | head -n 1
)}

[[ -n ${translation_root} && -d ${translation_root} ]] || { echo "[ERROR] Could not resolve translation root" >&2; exit 1; }
[[ -n ${axis_root} && -d ${axis_root} ]] || { echo "[ERROR] Could not resolve axis root" >&2; exit 1; }

out_root=${AXIS_REFINE_ROOT:-${repo_root}/tmp/native_sn_axis_gain_refine_all108_$(date -u +%Y%m%dT%H%M%SZ)}
launcher=${repo_root}/Scripts/launch_native_sn_axis_gain_refine_chunks.sh

mkdir -p "${out_root}"

merged_summary="${out_root}/current_best_summary_all108.tsv"
current_dir="${out_root}/current_best_masks"
subjects_file="${out_root}/subjects_all108.txt"

python3 - "${translation_root}" "${axis_root}" "${merged_summary}" "${current_dir}" "${subjects_file}" <<'PY'
import csv
import re
import sys
from pathlib import Path


translation_root = Path(sys.argv[1])
axis_root = Path(sys.argv[2])
summary_path = Path(sys.argv[3])
current_dir = Path(sys.argv[4])
subjects_path = Path(sys.argv[5])
current_dir.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gain_from_stage(stage: str) -> float:
    if stage == "pre":
        return 0.0
    match = re.fullmatch(r"g(\d+)", stage)
    if not match:
        raise ValueError(f"Unrecognized translation stage: {stage}")
    return int(match.group(1)) / 100.0


def best_mask_paths(root: Path):
    paths = {}
    for path in sorted(root.glob("chunks/c*/best/sn_groupmask_in_*.nii.gz")):
        subject = path.stem.replace("sn_groupmask_in_", "")
        if subject.endswith(".nii"):
            subject = subject[:-4]
        paths[subject] = path
    for path in sorted(root.glob("chunks/c*/best/sn_groupmask_in_*.nii")):
        subject = path.stem.replace("sn_groupmask_in_", "")
        paths.setdefault(subject, path)
    return paths


translation_rows = {}
for path in sorted(translation_root.glob("chunks/c*/best_summary.tsv")):
    for row in read_tsv(path):
        subject = row["subject"]
        stage = row["best_stage"]
        best_delta = float(row["best_delta_mm"])
        pre_delta = float(row["pre_delta_mm"])
        gain = gain_from_stage(stage)
        translation_rows[subject] = {
            "subject": subject,
            "accepted_stage": stage,
            "accepted_delta_mm": f"{best_delta:.6f}",
            "pre_delta_mm": f"{pre_delta:.6f}",
            "improvement_mm": f"{pre_delta - best_delta:.6f}",
            "best_candidate_stage": stage,
            "best_candidate_delta_mm": f"{best_delta:.6f}",
            "best_gx": f"{gain:.3f}",
            "best_gy": f"{gain:.3f}",
            "best_gz": f"{gain:.3f}",
            "accepted_improvement": "1" if stage != "pre" else "0",
        }

axis_rows = {}
for path in sorted(axis_root.glob("chunks/c*/best_axis_gain_summary.tsv")):
    for row in read_tsv(path):
        axis_rows[row["subject"]] = row

translation_masks = best_mask_paths(translation_root)
axis_masks = best_mask_paths(axis_root)

merged = dict(translation_rows)
merged.update(axis_rows)

for subject, row in merged.items():
    src = axis_masks.get(subject) if subject in axis_rows else translation_masks.get(subject)
    if src is None:
        raise FileNotFoundError(f"Missing current best mask for subject {subject}")
    dst = current_dir / src.name
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src)

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

subjects = sorted(merged, key=int)
with summary_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for subject in subjects:
        writer.writerow({key: merged[subject][key] for key in fieldnames})

subjects_path.write_text("\n".join(subjects) + "\n", encoding="utf-8")

print(summary_path)
print(subjects_path)
print(f"count={len(subjects)}")
print(f"translation_only={sum(1 for s in subjects if s not in axis_rows)}")
print(f"axis_override={sum(1 for s in subjects if s in axis_rows)}")
PY

AXIS_REFINE_ROOT="${out_root}" \
AXIS_REFINE_PREV_ROOT="${axis_root}" \
AXIS_REFINE_PREVIOUS_SUMMARY="${merged_summary}" \
AXIS_REFINE_CURRENT_DIR="${current_dir}" \
AXIS_REFINE_SUBJECTS_FILE="${subjects_file}" \
bash "${launcher}"
