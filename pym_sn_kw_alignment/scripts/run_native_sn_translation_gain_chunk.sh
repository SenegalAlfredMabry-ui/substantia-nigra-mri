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

subjects_file=${RIGID_RESCUE_SUBJECTS_FILE:?RIGID_RESCUE_SUBJECTS_FILE required}
manual_dir=${RIGID_RESCUE_MANUAL_DIR:?RIGID_RESCUE_MANUAL_DIR required}
input_dir=${RIGID_RESCUE_INPUT_DIR:?RIGID_RESCUE_INPUT_DIR required}
out_root=${RIGID_RESCUE_OUTPUT_ROOT:?RIGID_RESCUE_OUTPUT_ROOT required}
gain_list=${RIGID_RESCUE_GAINS:-"0.5 1.0 1.5 2.0"}

mkdir -p "${out_root}/pre" "${out_root}/best" "${out_root}/tmp"

calc_delta() {
  python3 - "$1" "$2" <<'PY'
import sys
from pathlib import Path
import nibabel as nib
import numpy as np

def load_mask(p):
    img = nib.load(str(p))
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

stage_name_for_gain() {
  python3 - "$1" <<'PY'
import sys
gain = float(sys.argv[1])
print(f"g{int(round(gain * 100)):03d}")
PY
}

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

  for gain in ${gain_list}; do
    stage=$(stage_name_for_gain "${gain}")
    mkdir -p "${out_root}/${stage}"
    read -r sdx sdy sdz <<<"$(python3 - "${dx}" "${dy}" "${dz}" "${gain}" <<'PY'
import sys
dx, dy, dz, gain = map(float, sys.argv[1:5])
print(f"{dx*gain:.6f} {dy*gain:.6f} {dz*gain:.6f}")
PY
)"
    trans_mat="${out_root}/tmp/translate_${stage}_${sid}.1D"
    printf "1 0 0 %s 0 1 0 %s 0 0 1 %s\n" "${sdx}" "${sdy}" "${sdz}" > "${trans_mat}"
    trans_out="${out_root}/${stage}/sn_groupmask_in_${sid}.nii.gz"
    3dAllineate -overwrite -base "${pre_mask}" -input "${pre_mask}" \
      -1Dmatrix_apply "${trans_mat}" -prefix "${trans_out}" -final NN -float >/dev/null 2>&1
  done
done < "${subjects_file}"

python3 Scripts/check_sn_alignment_metrics.py --subjects-file "${subjects_file}" \
  --manual-dir "${manual_dir}" --baseline-dir "${out_root}/pre" --output "${out_root}/pre_metrics.tsv"

stage_args=()
for gain in ${gain_list}; do
  stage=$(stage_name_for_gain "${gain}")
  python3 Scripts/check_sn_alignment_metrics.py --subjects-file "${subjects_file}" \
    --manual-dir "${manual_dir}" --baseline-dir "${out_root}/${stage}" --output "${out_root}/${stage}_metrics.tsv"
  stage_args+=("${stage}")
done

python3 - "${out_root}" "${stage_args[@]}" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
stages = ["pre", *sys.argv[2:]]
tables = {}
for stage in stages:
    path = root / f"{stage}_metrics.tsv"
    if not path.exists():
        continue
    with path.open(newline="", encoding="utf-8") as f:
        tables[stage] = {r["subject"]: r for r in csv.DictReader(f, delimiter="\t")}

subjects = sorted({s for t in tables.values() for s in t}, key=int)
best_dir = root / "best"
best_dir.mkdir(exist_ok=True)
summary = root / "best_summary.tsv"
with summary.open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["subject", "best_stage", "best_delta_mm", *[f"{stage}_delta_mm" for stage in stages]])
    for sub in subjects:
        best_stage = None
        best_val = None
        vals = {}
        for stage in stages:
            row = tables.get(stage, {}).get(sub)
            if not row:
                continue
            try:
                val = float(row["delta_mag_mm"])
            except Exception:
                continue
            vals[stage] = val
            if best_val is None or val < best_val:
                best_val = val
                best_stage = stage
        if best_stage is None:
            continue
        src = root / best_stage / f"sn_groupmask_in_{sub}.nii.gz"
        if not src.exists():
            src = root / best_stage / f"sn_groupmask_in_{sub}.nii"
        if src.exists():
            dst = best_dir / src.name
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            dst.symlink_to(src)
        row = [sub, best_stage, f"{best_val:.6f}"]
        for stage in stages:
            row.append(f"{vals.get(stage, float('nan')):.6f}")
        w.writerow(row)
print(f"WROTE {summary}")
PY
