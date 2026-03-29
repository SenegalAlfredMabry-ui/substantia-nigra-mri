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

mkdir -p "${out_root}"/{pre,translation,rigid,best,tmp}

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

  trans_mat="${out_root}/tmp/translate_${sid}.1D"
  printf "1 0 0 %s 0 1 0 %s 0 0 1 %s\n" "${dx}" "${dy}" "${dz}" > "${trans_mat}"
  trans_out="${out_root}/translation/sn_groupmask_in_${sid}.nii.gz"
  3dAllineate -overwrite -base "${pre_mask}" -input "${pre_mask}" \
    -1Dmatrix_apply "${trans_mat}" -prefix "${trans_out}" -final NN -float >/dev/null 2>&1

  rigid_out="${out_root}/rigid/sn_groupmask_in_${sid}.nii.gz"
  rigid_mat="${out_root}/tmp/rigid_${sid}.1D"
  3dAllineate -overwrite -base "${manual_mask}" -source "${trans_out}" \
    -prefix "${rigid_out}" -final NN -float -warp shift_rotate -cost mi \
    -cmass -1Dmatrix_save "${rigid_mat}" >/dev/null 2>&1 || {
      echo "[WARN] ${sid}: rigid refine failed; keeping translation"
      rm -f "${rigid_out}"
    }
done < "${subjects_file}"

python3 Scripts/check_sn_alignment_metrics.py --subjects-file "${subjects_file}" \
  --manual-dir "${manual_dir}" --baseline-dir "${out_root}/pre" --output "${out_root}/pre_metrics.tsv"
python3 Scripts/check_sn_alignment_metrics.py --subjects-file "${subjects_file}" \
  --manual-dir "${manual_dir}" --baseline-dir "${out_root}/translation" --output "${out_root}/translation_metrics.tsv"
python3 Scripts/check_sn_alignment_metrics.py --subjects-file "${subjects_file}" \
  --manual-dir "${manual_dir}" --baseline-dir "${out_root}/rigid" --output "${out_root}/rigid_metrics.tsv"

python3 - "${out_root}" <<'PY'
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])
tables = {}
for stage in ("pre", "translation", "rigid"):
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
    w.writerow(["subject", "best_stage", "best_delta_mm", "pre_delta_mm", "translation_delta_mm", "rigid_delta_mm"])
    for sub in subjects:
        best_stage = None
        best_val = None
        vals = {}
        for stage in ("pre", "translation", "rigid"):
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
        w.writerow([sub, best_stage, f"{best_val:.6f}", f"{vals.get('pre', float('nan')):.6f}", f"{vals.get('translation', float('nan')):.6f}", f"{vals.get('rigid', float('nan')):.6f}"])
print(f"WROTE {summary}")
PY
