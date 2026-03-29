#!/bin/bash -l
# antpunk_v1.sh — hybrid sswarper + ANTs warp refinement with LC/SN emphasis

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SNHEART_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)
repo_root_default=${SNHEART_ROOT:-/path/to/SNHEART}
if [[ ! -d ${SNHEART_ROOT}/MRI_data ]]; then
  SNHEART_ROOT=${repo_root_default}
fi

usage() {
  cat <<'USAGE'
Usage: antpunk_v1.sh [options]

Options:
  --subjects <file|list>   Subject roster file or comma/space list (default MRI_data/locations/TSE_SN_ASR_SUBLIST.txt)
  --mode <full|refine_only|qc_only>
  --roi-padding-mm <mm>    ROI bounding-box padding in mm (legacy alias)
  --roi-bbox-pad-mm <mm>   Explicit padding (mm) around the SN ROI box (default ANTPUNK_ROI_BBOX_PAD_MM)
  --roi-bbox-pad-vox <n>   Override padding using template voxels instead of mm
  --roi-erode-mm <mm>      Erode the template-space subject ROI after warping (default 0)
  --roi-erode-vox <n>      Explicit erosion in voxels (overrides mm)
  --roi-strategy <name>    ROI source: atlas or intensity (default atlas)
  --ants-config <path>     JSON/YAML with custom antsRegistration stages (optional)
  --sswarper-root <path>   sswarper warp cache (default MRI_data/WARPS_supervillain)
  --output-root <path>     Destination for antpunk warps (default MRI_data/WARPS_antpunk)
  --log-root <path>        Logs/QC directory (default MRI_data/ants_processing/antpunk_logs)
  --run-sswarper           Auto-run supervillan_wallet.sh when sswarper outputs missing
  --qc-threshold <val>     Minimum Dice vs hand SN mask to accept antpunk (default 0.80)
  --qc-improve <val>       Required Dice improvement over sswarper (default 0.02)
  --roi-only               Regenerate subject ROI masks only (skip antsRegistration)
  --help                   Show this message

Environment overrides mirror the flags (ANTPUNK_SUBJECTS, ANTPUNK_MODE, etc.).
USAGE
}

log() { printf '%s\n' "$*"; }
log_info() { log "[INFO] $*"; }
log_warn() { log "[WARN] $*" >&2; }
die() { log "[ERROR] $*" >&2; exit 1; }

ANTPUNK_DISABLE_ROI_CROP=${ANTPUNK_DISABLE_ROI_CROP:-1}
ANTPUNK_DISABLE_SUBJECT_MASK=${ANTPUNK_DISABLE_SUBJECT_MASK:-1}
disable_roi_crop=${ANTPUNK_DISABLE_ROI_CROP}
disable_subject_mask=${ANTPUNK_DISABLE_SUBJECT_MASK}
log_info "antpunk_v1: disable_roi_crop=${disable_roi_crop} disable_subject_mask=${disable_subject_mask}"

python_site_default="${SNHEART_ROOT}/venv/lib/python3.9/site-packages"
python_site=${ANTPUNK_PYTHON_SITE:-${python_site_default}}
if [[ -d ${python_site} ]]; then
  if [[ -z ${PYTHONPATH:-} ]]; then
    export PYTHONPATH="${python_site}"
  else
    export PYTHONPATH="${python_site}:${PYTHONPATH}"
  fi
else
  log_warn "Python dependencies directory ${python_site} missing; packaging import may fail"
fi

ants_bin_default="${SNHEART_ROOT}/Software/ants-2.6.2/bin"
ants_bin=${ANTPUNK_ANTS_BIN:-${ants_bin_default}}
if [[ -d ${ants_bin} ]]; then
  export ANTSPATH="${ants_bin}"
  export PATH="${ANTSPATH}:${PATH}"
else
  log_warn "ANTs bin ${ants_bin} missing; antsRegistration may be unavailable"
fi

if ! command -v 3dAllineate >/dev/null 2>&1; then
  if command -v module >/dev/null 2>&1; then
    if module load afni >/dev/null 2>&1; then
      log_info "Loaded afni module for ROI helpers"
    else
      log_warn "module load afni failed; AFNI utilities may be unavailable"
    fi
  fi
fi
if ! command -v 3dAllineate >/dev/null 2>&1; then
  log_warn "AFNI binaries (3dAllineate) not found in PATH; ROI prep may fail"
fi

subjects_arg=${ANTPUNK_SUBJECTS:-}
mode=${ANTPUNK_MODE:-full}
roi_padding_mm=${ANTPUNK_ROI_PADDING_MM:-22}
ants_config=${ANTPUNK_ANTS_CONFIG:-}
sswarper_root=${ANTPUNK_SSWARPER_ROOT:-${SNHEART_ROOT}/MRI_data/WARPS_supervillain}
output_root=${ANTPUNK_OUTPUT_ROOT:-${SNHEART_ROOT}/MRI_data/WARPS_antpunk}
log_root=${ANTPUNK_LOG_ROOT:-${SNHEART_ROOT}/MRI_data/ants_processing/antpunk_logs}
template_ref=${ANTPUNK_TEMPLATE:-${SNHEART_ROOT}/MRI_data/MNI152_2009_template_SSW.nii.gz}
sn_template=${ANTPUNK_SN_TEMPLATE:-${SNHEART_ROOT}/MRI_data/ROIs/SN_probabilistic_MNI152/SN_probabilistic.nii}
lc_template=${ANTPUNK_LC_TEMPLATE:-}
manifest_path=${ANTPUNK_MANIFEST:-${output_root}/warp_manifest_antpunk.tsv}
run_sswarper=${ANTPUNK_RUN_SSWARPER:-0}
qc_threshold=${ANTPUNK_QC_THRESHOLD:-0.80}
qc_improvement=${ANTPUNK_QC_IMPROVEMENT:-0.02}
qc_max_delta_mm=${ANTPUNK_QC_MAX_DELTA_MM:-3.0}
force_qc_status=${ANTPUNK_FORCE_QC_STATUS:-}
roi_dilate_mm=${ANTPUNK_ROI_DILATE_MM:-0}
roi_dilate_vox=${ANTPUNK_ROI_DILATE_VOX:-}
roi_bbox_pad_mm=${ANTPUNK_ROI_BBOX_PAD_MM:-${roi_padding_mm}}
roi_bbox_pad_vox=${ANTPUNK_ROI_BBOX_PAD_VOX:-"20 8 22"}
roi_erode_mm=${ANTPUNK_ROI_ERODE_MM:-0}
roi_erode_vox=${ANTPUNK_ROI_ERODE_VOX:-}
roi_only=${ANTPUNK_ROI_ONLY:-0}
roi_target_vox=${ANTPUNK_ROI_TARGET_VOX:-}
roi_offset_json=${ANTPUNK_ROI_OFFSET_JSON:-}
roi_strategy=${ANTPUNK_ROI_STRATEGY:-atlas}
intensity_top_vox=${ANTPUNK_INTENSITY_TOP_VOXELS:-200}
intensity_smooth_sigma=${ANTPUNK_INTENSITY_SMOOTH_SIGMA:-1.0}
shift_roi_script="$SCRIPT_DIR/shift_roi_to_target.py"
if [[ ! -x ${shift_roi_script} && -n ${SLURM_SUBMIT_DIR:-} ]]; then
  alt_shift="${SLURM_SUBMIT_DIR}/Scripts/shift_roi_to_target.py"
  if [[ -x ${alt_shift} ]]; then
    shift_roi_script="${alt_shift}"
  fi
fi
if [[ -n ${roi_target_vox} && ! -x ${shift_roi_script} ]]; then
  log_warn "shift_roi_to_target.py not executable; disabling ROI centroid alignment"
  roi_target_vox=""
fi

motip_root=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data

case ${roi_strategy} in
  atlas|intensity)
    ;;
  *)
    die "Unsupported ROI strategy: ${roi_strategy}"
    ;;
esac
if [[ ${disable_subject_mask} -ne 0 ]]; then
  log_warn "ROI guidance disabled via ANTPUNK_DISABLE_SUBJECT_MASK; ROI strategy ${roi_strategy} will not influence antsRegistration"
fi

while [[ $# -gt 0 ]]; do
  case $1 in
    --subjects)
      subjects_arg=$2; shift 2 ;;
    --mode)
      mode=$2; shift 2 ;;
    --roi-padding-mm)
      roi_padding_mm=$2
      roi_bbox_pad_mm=$2
      shift 2 ;;
    --roi-bbox-pad-mm)
      roi_bbox_pad_mm=$2; shift 2 ;;
    --roi-bbox-pad-vox)
      roi_bbox_pad_vox=$2; shift 2 ;;
    --roi-erode-mm)
      roi_erode_mm=$2; shift 2 ;;
    --roi-erode-vox)
      roi_erode_vox=$2; shift 2 ;;
    --roi-strategy)
      roi_strategy=$2; shift 2 ;;
    --ants-config)
      ants_config=$2; shift 2 ;;
    --sswarper-root)
      sswarper_root=$2; shift 2 ;;
    --output-root)
      output_root=$2; shift 2 ;;
    --log-root)
      log_root=$2; shift 2 ;;
    --run-sswarper)
      run_sswarper=1; shift ;;
    --qc-threshold)
      qc_threshold=$2; shift 2 ;;
    --qc-improve)
      qc_improvement=$2; shift 2 ;;
    --roi-only)
      roi_only=1; shift ;;
    --help|-h)
      usage; exit 0 ;;
    *)
      die "Unknown option: $1" ;;
  esac
done

if [[ -z ${roi_bbox_pad_mm} ]]; then
  roi_bbox_pad_mm=${roi_padding_mm}
fi
if [[ -z ${roi_bbox_pad_mm} ]]; then
  roi_bbox_pad_mm=40
fi
if [[ -z ${roi_dilate_vox} ]]; then
  roi_dilate_vox=$(ANTPUNK_TEMPLATE_PATH="${template_ref}" ANTPUNK_ROI_DILATE_MM="${roi_dilate_mm}" python - <<'PY'
import os, nibabel as nib, numpy as np
tpl=os.environ.get('ANTPUNK_TEMPLATE_PATH')
mm=float(os.environ.get('ANTPUNK_ROI_DILATE_MM','0') or 0)
if not tpl or mm <= 0:
    print(0)
else:
    vox=nib.affines.voxel_sizes(nib.load(tpl).affine)
    mean_vox=float(np.mean(vox)) if vox.size else 1.0
    steps=max(1, int(round(mm / max(mean_vox, 1e-3))))
    print(steps)
PY
)
fi
roi_dilate_vox=${roi_dilate_vox:-0}

if [[ -z ${roi_erode_vox} ]]; then
  roi_erode_vox=$(ANTPUNK_TEMPLATE_PATH="${template_ref}" ANTPUNK_ROI_ERODE_MM="${roi_erode_mm}" python - <<'PY'
import os, nibabel as nib, numpy as np
tpl=os.environ.get('ANTPUNK_TEMPLATE_PATH')
mm=float(os.environ.get('ANTPUNK_ROI_ERODE_MM','0') or 0)
if not tpl or mm <= 0:
    print(0)
else:
    vox=nib.affines.voxel_sizes(nib.load(tpl).affine)
    mean_vox=float(np.mean(vox)) if vox.size else 1.0
    steps=max(1, int(round(mm / max(mean_vox, 1e-3))))
    print(steps)
PY
)
fi
roi_erode_vox=${roi_erode_vox:-0}

subjects_file_default="${SNHEART_ROOT}/MRI_data/locations/TSE_SN_ASR_SUBLIST.txt"

command -v antsRegistration >/dev/null 2>&1 || die "antsRegistration not found in PATH"

[[ -f ${template_ref} ]] || die "Template volume ${template_ref} missing"
[[ -f ${sn_template} ]] || die "SN template ${sn_template} missing"

roi_cache_dir=${output_root}/.antpunk_cache
roi_bounds_file=${roi_cache_dir}/roi_bounds.json
mkdir -p "${output_root}" "${log_root}" "${roi_cache_dir}"
atlas_threshold=${ANTPUNK_ATLAS_THRESHOLD:-0.15}
if [[ ${disable_subject_mask} -ne 0 ]]; then
  log_info "Subject SN mask ROI guidance disabled (ANTPUNK_DISABLE_SUBJECT_MASK=${disable_subject_mask})"
fi

load_subjects() {
  local source=$1
  local -a subs=()
  if [[ -z ${source} ]]; then
    source=${subjects_file_default}
  fi
  if [[ -f ${source} ]]; then
    while IFS= read -r line; do
      [[ -z ${line} || ${line} =~ ^# ]] && continue
      subs+=(${line})
    done <"${source}"
  else
    source=${source//,/ }
    for token in ${source}; do
      [[ -z ${token} ]] && continue
      subs+=(${token})
    done
  fi
  printf '%s\n' "${subs[@]}"
}

subjects=($(load_subjects "${subjects_arg}"))
[[ ${#subjects[@]} -gt 0 ]] || die "No subjects provided"

manifest_init() {
  if [[ ! -f ${manifest_path} ]]; then
    printf 'subject\tmode\troi_strategy\tsswarper_warp\tantpunk_warp\tants_log\tqc_status\tnotes\n' >"${manifest_path}"
  fi
}

run_qc() {
  local sub=$1
  local out_dir="${output_root}/${sub}"
  local metrics_file="${out_dir}/qc_metrics.tsv"
  if [[ ! -f ${metrics_file} ]]; then
    printf 'subject\tdice_intensity_antpunk\tdice_intensity_ss\tdice_kw_antpunk\tdice_kw_ss\n' >"${metrics_file}"
  fi
  local antpunk_sn="${out_dir}/sn_template_inTSE_antpunk.nii.gz"
  local ss_sn="${out_dir}/sn_template_inTSE_sswarper.nii.gz"
  if ! warp_sn_to_tse "${sub}" antpunk "${antpunk_sn}"; then
    log_warn "${sub}: failed to warp SN template with antpunk"
    return
  fi
  if ! warp_sn_to_tse "${sub}" sswarper "${ss_sn}"; then
    log_warn "${sub}: failed to warp SN template with sswarper"
    return
  fi
  local qc_base_mask="${out_dir}/qc_roi_base_in_tse.nii.gz"
  if ! python3 - "${ss_sn}" "${qc_base_mask}" "${atlas_threshold}" <<'PYINNER'
import nibabel as nib, numpy as np, sys
src, dst, thr = sys.argv[1:4]
thr = float(thr)
img = nib.load(src)
data = img.get_fdata()
mask = (data >= thr).astype(np.uint8)
if not mask.any():
    nz = data[data > 0]
    if nz.size:
        cut = float(np.percentile(nz, 95))
        mask = (data >= cut).astype(np.uint8)
if not mask.any():
    sys.exit(1)
nib.Nifti1Image(mask, img.affine, img.header).to_filename(dst)
PYINNER
  then
    log_warn "${sub}: failed to build atlas-derived QC ROI"
    return
  fi
  local qc_intensity_mask=""
  local qc_target_mask="${qc_base_mask}"
  local qc_intensity_path="${out_dir}/qc_roi_intensity_in_tse.nii.gz"
  if build_intensity_roi_mask "${sub}" "${qc_base_mask}" "${qc_intensity_path}"; then
    qc_intensity_mask="${qc_intensity_path}"
    qc_target_mask="${qc_intensity_mask}"
  else
    log_warn "${sub}: intensity-derived QC ROI failed; using atlas mask"
  fi

  local dice_int_antpunk dice_int_ss
  dice_int_antpunk=$(compute_dice "${antpunk_sn}" "${qc_target_mask}")
  dice_int_ss=$(compute_dice "${ss_sn}" "${qc_target_mask}")

  local manual_mask="${motip_root}/TSE/SN/SN_ROI_KW_${sub}_ZP.nii"
  [[ -f ${manual_mask} ]] || manual_mask="${motip_root}/TSE/SN/SN_ROI_KW_${sub}.nii"
  local dice_kw_antpunk="nan"
  local dice_kw_ss="nan"
  if [[ -f ${manual_mask} ]]; then
    dice_kw_antpunk=$(compute_dice "${antpunk_sn}" "${manual_mask}")
    dice_kw_ss=$(compute_dice "${ss_sn}" "${manual_mask}")
  else
    log_warn "${sub}: manual SN ROI missing; skipping KW Dice"
  fi
  printf '%s\t%s\t%s\t%s\t%s\n' "${sub}" "${dice_int_antpunk}" "${dice_int_ss}" "${dice_kw_antpunk}" "${dice_kw_ss}" >>"${metrics_file}"

  local cm_delta_mm="nan"
  if [[ -f ${manual_mask} ]]; then
    cm_delta_mm=$(python3 - "${antpunk_sn}" "${manual_mask}" <<'PY'
import math, sys, nibabel as nib, numpy as np
sn_path, manual_path = sys.argv[1:3]
def centroid_mm(path):
    img = nib.load(path)
    data = img.get_fdata()
    mask = data > 0
    if not mask.any():
        return None
    coords = np.argwhere(mask)
    cm_vox = coords.mean(axis=0)
    return nib.affines.apply_affine(img.affine, cm_vox)
sn = centroid_mm(sn_path)
manual = centroid_mm(manual_path)
if sn is None or manual is None:
    print("nan")
else:
    delta = np.linalg.norm(np.asarray(sn) - np.asarray(manual))
    print(f"{delta:.6f}")
PY
) || cm_delta_mm="nan"
  fi

  local status notes
  notes="dice_intensity_antpunk=${dice_int_antpunk};dice_intensity_ss=${dice_int_ss};cm_delta_mm=${cm_delta_mm}"
  if [[ ${dice_kw_antpunk} != "nan" ]]; then
    notes+=";dice_kw_antpunk=${dice_kw_antpunk};dice_kw_ss=${dice_kw_ss}"
  fi
  status=$(python - "$dice_int_antpunk" "$dice_int_ss" "$qc_threshold" "$qc_improvement" <<'PYINNER'
import sys, math
ant=sys.argv[1]
ss=sys.argv[2]
thr=float(sys.argv[3])
impr=float(sys.argv[4])
try:
    ant_val=float(ant)
except ValueError:
    print('UNKNOWN')
    sys.exit(0)
try:
    ss_val=float(ss)
except ValueError:
    ss_val=float('nan')
if math.isnan(ant_val):
    print('UNKNOWN')
elif ant_val < thr:
    print('SSWARP')
elif not math.isnan(ss_val) and (ant_val - ss_val) < impr:
    print('SSWARP')
else:
    print('ANTPUNK')
PYINNER
)
  if [[ ${status} != "ANTPUNK" && ${cm_delta_mm} != "nan" ]]; then
    if python3 - "${cm_delta_mm}" "${qc_max_delta_mm}" <<'PY'
import sys
delta=float(sys.argv[1])
thr=float(sys.argv[2])
sys.exit(0 if delta <= thr else 1)
PY
    then
      log_info "${sub}: centroid delta ${cm_delta_mm} mm <= ${qc_max_delta_mm} mm; promoting to ANTPUNK"
      status="ANTPUNK"
      notes+=";cm_delta_override=1"
    fi
  fi
  if [[ -n ${force_qc_status} ]]; then
    log_warn "${sub}: overriding QC status ${status} -> ${force_qc_status}"
    status=${force_qc_status}
  fi
  printf '%s' "${status}" >"${out_dir}/qc_status.txt"
  printf '%s' "${notes}" >"${out_dir}/qc_notes.txt"
}

ensure_sswarper_subject() {
  local sub=$1
  local ss_dir="${sswarper_root}/${sub}"
  local anat="${ss_dir}/anatQQ.${sub}.nii"
  if [[ -f ${anat} ]]; then
    return 0
  fi
  if [[ ${run_sswarper} -ne 1 ]]; then
    log_warn "${ss_dir} missing; rerun supervillan_wallet.sh manually or pass --run-sswarper"
    return 1
  fi
  if [[ ! -x ${SCRIPT_DIR}/supervillan_wallet.sh ]]; then
    log_warn "supervillan_wallet.sh not executable; cannot rebuild sswarper cache"
    return 1
  fi
  log_info "Running supervillan_wallet.sh for ${sub}"
  SUPERVILLAIN_SUBJECTS="${sub}" bash "${SCRIPT_DIR}/supervillan_wallet.sh" || {
    log_warn "supervillan_wallet.sh failed for ${sub}"
    return 1
  }
  [[ -f ${anat} ]] || { log_warn "sswarper outputs still missing after rerun"; return 1; }
}

find_tse() {
  local sub=$1
  local patterns=(
    "${motip_root}/TSE/zeropad_tse_${sub}.nii"
    "${motip_root}/TSE/zeropad_tse_${sub}.nii.gz"
    "${motip_root}/TSE/${sub}_TSE_zeropad.nii"
    "${motip_root}/TSE/${sub}_TSE_zeropad.nii.gz"
    "${motip_root}/TSE/${sub}_TSE.nii"
    "${motip_root}/TSE/${sub}_TSE.nii.gz"
  )
  for p in "${patterns[@]}"; do
    if [[ -f ${p} ]]; then
      printf '%s' "${p}"
      return 0
    fi
  done
  return 1
}

find_corrected_tse() {
  local sub=$1
  local patterns=(
    "${motip_root}/TSE/SN/corrected_tse_DC_${sub}_ptfix.nii"
    "${motip_root}/TSE/SN/corrected_tse_DC_${sub}_ptfix.nii.gz"
    "${motip_root}/TSE/SN/corrected_tse_DC_${sub}.nii"
    "${motip_root}/TSE/SN/corrected_tse_DC_${sub}.nii.gz"
    "${motip_root}/TSE/SN/corrected_tse_${sub}.nii"
    "${motip_root}/TSE/SN/corrected_tse_${sub}.nii.gz"
  )
  for p in "${patterns[@]}"; do
    if [[ -f ${p} ]]; then
      printf '%s' "${p}"
      return 0
    fi
  done
  return 1
}

compute_roi_bounds() {
  local bounds_file=$1
  python - "$sn_template" "$bounds_file" "$roi_bbox_pad_mm" <<'PY'
import itertools, json, math, sys
import nibabel as nib
import numpy as np
sn_path, out_path, pad_mm = sys.argv[1], sys.argv[2], float(sys.argv[3])
img = nib.load(sn_path)
data = img.get_fdata()
mask = data > 0.2
if not mask.any():
    raise SystemExit("SN template mask empty")
coords = mask.nonzero()
mins = [int(c.min()) for c in coords]
maxs = [int(c.max())+1 for c in coords]
vox_sizes = nib.affines.voxel_sizes(img.affine)
pad = [int(math.ceil(pad_mm / vs)) for vs in vox_sizes]
shape = data.shape
start = [max(0, mins[i]-pad[i]) for i in range(3)]
stop = [min(shape[i], maxs[i]+pad[i]) for i in range(3)]
aff = img.affine
corner_indices = []
for i in (start[0], max(start[0], stop[0]-1)):
    for j in (start[1], max(start[1], stop[1]-1)):
        for k in (start[2], max(start[2], stop[2]-1)):
            corner_indices.append((i, j, k))
corner_indices = np.array(corner_indices, dtype=float)
corner_world = nib.affines.apply_affine(aff, corner_indices)
world_min = corner_world.min(axis=0)
world_max = corner_world.max(axis=0)
with open(out_path, 'w') as f:
    json.dump({
        'start': start,
        'stop': stop,
        'world_min': world_min.tolist(),
        'world_max': world_max.tolist()
    }, f)
PY
}

crop_volume() {
  local input_vol=$1
  local output_vol=$2
  local bounds_file=$3
  python - "$input_vol" "$output_vol" "$bounds_file" <<'PY'
import itertools, json, math, sys
import nibabel as nib
import numpy as np
inp, outp, bounds_path = sys.argv[1:]
with open(bounds_path) as f:
    bounds = json.load(f)
img = nib.load(inp)
data = img.get_fdata()
while data.ndim > 3:
    data = data[..., 0]
shape = data.shape
aff = img.affine
if 'world_min' in bounds and 'world_max' in bounds:
    world_min = np.array(bounds['world_min'], dtype=float)
    world_max = np.array(bounds['world_max'], dtype=float)
    inv_aff = np.linalg.inv(aff)
    corners = []
    for x in (world_min[0], world_max[0]):
        for y in (world_min[1], world_max[1]):
            for z in (world_min[2], world_max[2]):
                corners.append(nib.affines.apply_affine(inv_aff, (x, y, z)))
    corners = np.array(corners)
    starts = np.floor(corners.min(axis=0)).astype(int)
    stops = np.ceil(corners.max(axis=0)).astype(int) + 1
    start = [max(0, s) for s in starts[:3]]
    stop = [min(shape[i], stops[i]) for i in range(3)]
else:
    start = bounds['start']
    stop = bounds['stop']
slices = tuple(slice(start[i], stop[i]) for i in range(3))
sub = data[slices]
new_aff = aff.copy()
offset = new_aff[:3,:3] @ start + new_aff[:3,3]
new_aff[:3,3] = offset
nib.Nifti1Image(sub, new_aff).to_filename(outp)
PY
}

embed_warp() {
  local roi_warp=$1
  local out_warp=$2
  antsApplyTransforms -d 3 -o "[${out_warp},1]" -t "${roi_warp}" -r "${template_ref}" >/dev/null 2>&1 || return 1
  return 0
}

embed_roi_volume() {
  local roi_vol=$1
  local ref_vol=$2
  local bounds_file=$3
  local out_vol=$4
  python - "${roi_vol}" "${ref_vol}" "${bounds_file}" "${out_vol}" <<'PY'
import json, sys
import nibabel as nib
import numpy as np
roi_path, ref_path, bounds_path, out_path = sys.argv[1:]
roi_img = nib.load(roi_path)
roi_data = roi_img.get_fdata()
while roi_data.ndim > 3:
    roi_data = roi_data[..., 0]
ref_img = nib.load(ref_path)
with open(bounds_path) as f:
    bounds = json.load(f)
aff = ref_img.affine
shape = ref_img.shape[:3]
if 'world_min' in bounds and 'world_max' in bounds:
    world_min = np.array(bounds['world_min'], dtype=float)
    world_max = np.array(bounds['world_max'], dtype=float)
    inv_aff = np.linalg.inv(aff)
    corners = []
    for x in (world_min[0], world_max[0]):
        for y in (world_min[1], world_max[1]):
            for z in (world_min[2], world_max[2]):
                corners.append(nib.affines.apply_affine(inv_aff, (x, y, z)))
    corners = np.array(corners)
    starts = np.floor(corners.min(axis=0)).astype(int)
    stops = np.ceil(corners.max(axis=0)).astype(int) + 1
    start = [max(0, s) for s in starts[:3]]
    stop = [min(shape[i], stops[i]) for i in range(3)]
else:
    start = bounds['start']
    stop = bounds['stop']
expected = tuple(max(0, stop[i]-start[i]) for i in range(3))
embed = np.zeros(expected, dtype=roi_data.dtype)
min_shape = tuple(min(expected[i], roi_data.shape[i]) for i in range(3))
embed_slices = tuple(slice(0, m) for m in min_shape)
roi_slices = tuple(slice(0, m) for m in min_shape)
embed[embed_slices] = roi_data[roi_slices]
full = np.zeros(shape, dtype=roi_data.dtype)
x, y, z = start
x2 = min(shape[0], x + embed.shape[0])
y2 = min(shape[1], y + embed.shape[1])
z2 = min(shape[2], z + embed.shape[2])
full[x:x2, y:y2, z:z2] = embed[:(x2-x), :(y2-y), :(z2-z)]
nib.Nifti1Image(full, ref_img.affine).to_filename(out_path)
PY
}

ensure_template_bbox_mask() {
  local pad_mm_tag=${roi_bbox_pad_mm//[^0-9A-Za-z.+-]/}
  pad_mm_tag=${pad_mm_tag//./p}
  [[ -z ${pad_mm_tag} ]] && pad_mm_tag=auto
  local vox_tag=${roi_bbox_pad_vox//,/x}
  vox_tag=${vox_tag//[^0-9xX+-]/}
  if [[ -z ${vox_tag} ]]; then
    vox_tag=auto
  else
    vox_tag=${vox_tag,,}
  fi
  local mask_path="${roi_cache_dir}/sn_bbox_mask_padmm${pad_mm_tag}_padvox${vox_tag}.nii.gz"
  if [[ ! -f ${mask_path} ]]; then
    log_info "Building template ROI mask (pad=${roi_bbox_pad_mm} mm, vox_override=${roi_bbox_pad_vox:-auto})" >&2
    python - "$sn_template" "$mask_path" "${roi_bbox_pad_mm}" "${roi_bbox_pad_vox:-}" <<'PY'
import math, sys
import nibabel as nib
import numpy as np
tpl_path, out_path, pad_mm, pad_vox_raw = sys.argv[1:5]
img = nib.load(tpl_path)
data = img.get_fdata()
mask = data > 0.2
if not mask.any():
    raise SystemExit("SN probabilistic mask is empty")
coords = mask.nonzero()
mins = [int(c.min()) for c in coords]
maxs = [int(c.max()) + 1 for c in coords]
vox_sizes = nib.affines.voxel_sizes(img.affine)[:3]
pad_vox = None
raw = (pad_vox_raw or '').replace('X', 'x').replace(';', ',')
if raw.strip():
    tokens = [tok.strip() for tok in raw.replace('x', ',').split(',') if tok.strip()]
    if len(tokens) == 1:
        val = max(0, int(round(float(tokens[0]))))
        pad_vox = [val, val, val]
    elif len(tokens) >= 3:
        pad_vox = [max(0, int(round(float(tokens[i])))) for i in range(3)]
    elif len(tokens) == 2:
        vals = [max(0, int(round(float(tok)))) for tok in tokens]
        pad_vox = [vals[0], vals[1], vals[1]]
if pad_vox is None:
    try:
        pad_mm_val = float(pad_mm)
    except ValueError:
        pad_mm_val = 0.0
    pad_vox = [max(0, int(math.ceil(pad_mm_val / max(v, 1e-3)))) for v in vox_sizes]
shape = data.shape[:3]
start = [max(0, mins[i] - pad_vox[i]) for i in range(3)]
stop = [min(shape[i], maxs[i] + pad_vox[i]) for i in range(3)]
bbox = np.zeros(shape, dtype=np.uint8)
bbox[start[0]:stop[0], start[1]:stop[1], start[2]:stop[2]] = 1
nib.Nifti1Image(bbox, img.affine).to_filename(out_path)
PY
  fi
  printf '%s' "${mask_path}"
}

setup_roi_cache() {
  mkdir -p "${roi_cache_dir}"
  if [[ ! -f ${roi_bounds_file} ]]; then
    compute_roi_bounds "${roi_bounds_file}"
  fi
}
setup_roi_cache

template_bbox_mask=""
if [[ ${disable_roi_crop} -eq 0 ]]; then
  if template_bbox_mask=$(ensure_template_bbox_mask); then
    log_info "Template ROI mask: ${template_bbox_mask}"
  else
    template_bbox_mask=""
    log_warn "Failed to prepare template ROI mask; antsRegistration will run without fixed-mask guidance"
  fi
else
  log_info "ROI bbox cropping disabled; antsRegistration will use full template/anat volumes"
fi



ensure_subject_roi_template() {
  local sub=$1
  local out_dir="${output_root}/${sub}"
  mkdir -p "${out_dir}"
  local subj_roi_template="${out_dir}/subject_roi_in_template.nii.gz"
  local subj_roi_anat="${out_dir}/subject_roi_in_anat.nii.gz"
  if [[ -f ${subj_roi_template} && -f ${subj_roi_anat} ]]; then
    printf '%s\n%s\n' "${subj_roi_template}" "${subj_roi_anat}"
    return 0
  fi
  local ss_dir="${sswarper_root}/${sub}"
  local anat_img="${ss_dir}/anatSS.${sub}.nii"
  local anat_to_tse="${motip_root}/TSE/${sub}_anat_toTSE_mat.aff12.1D"
  local warp_field="${ss_dir}/anatQQ.${sub}_WARP.nii"
  local warp_affine="${ss_dir}/anatQQ.${sub}.aff12.1D"
  if [[ ! -f ${anat_img} || ! -f ${anat_to_tse} || ! -f ${warp_field} || ! -f ${warp_affine} ]]; then
    log_warn "${sub}: missing sswarper inputs for ROI generation"
    return 1
  fi
  local tmp_dir
  if [[ -n ${ANTPUNK_ROI_DEBUG_DIR:-} ]]; then
    tmp_dir="${ANTPUNK_ROI_DEBUG_DIR}/tmp_subject_roi_${sub}" && rm -rf "${tmp_dir}"
    mkdir -p "${tmp_dir}"
  else
    tmp_dir=$(mktemp -d "${out_dir}/tmp_subject_roi_${sub}.XXXX")
  fi
  local tse_to_anat_mat="${tmp_dir}/tse_to_anat_mat.aff12.1D"
  local cat_matvec_bin
  cat_matvec_bin=$(command -v cat_matvec || true)
  if [[ -n ${cat_matvec_bin} ]]; then
    if ! ${cat_matvec_bin} ${anat_to_tse} -I > "${tse_to_anat_mat}"; then
      log_warn "${sub}: cat_matvec failed for ${anat_to_tse}"
      [[ -z ${ANTPUNK_ROI_DEBUG_DIR:-} ]] && rm -rf "${tmp_dir}"
      return 1
    fi
  else
    if ! python3 - "${anat_to_tse}" "${tse_to_anat_mat}" <<'PYINNER'
import numpy as np, sys
src, dst = sys.argv[1], sys.argv[2]
vals = []
with open(src) as fh:
    for line in fh:
        vals.extend(line.split())
if len(vals) != 12:
    raise SystemExit(1)
mat = np.array(vals, dtype=float).reshape(3, 4)
rot = mat[:, :3]
trans = mat[:, 3]
try:
    rot_inv = np.linalg.inv(rot)
except np.linalg.LinAlgError:
    raise SystemExit(1)
shift = -rot_inv.dot(trans)
full = np.hstack([rot_inv, shift[:, None]]).reshape(-1)
with open(dst, 'w') as out:
    out.write('	'.join(f"{val: .6f}" for val in full))
    out.write('
')
PYINNER
    then
      log_warn "${sub}: python fallback failed to invert ${anat_to_tse}"
      [[ -z ${ANTPUNK_ROI_DEBUG_DIR:-} ]] && rm -rf "${tmp_dir}"
      return 1
    fi
  fi

  local base_tse="${tmp_dir}/sn_roi_base_in_tse.nii.gz"
  if ! warp_sn_to_tse "${sub}" sswarper "${base_tse}"; then
    log_warn "${sub}: failed to warp SN template into TSE space"
    [[ -z ${ANTPUNK_ROI_DEBUG_DIR:-} ]] && rm -rf "${tmp_dir}"
    return 1
  fi
  local base_mask_tse="${tmp_dir}/sn_roi_mask_in_tse.nii.gz"
  if ! python3 - "${base_tse}" "${base_mask_tse}" "${atlas_threshold}" <<'PYINNER'
import nibabel as nib, numpy as np, sys
src, dst, thr = sys.argv[1:4]
thr = float(thr)
img = nib.load(src)
data = (img.get_fdata() >= thr).astype(np.uint8)
if not data.any():
    flat = img.get_fdata().ravel()
    nz = flat[flat > 0]
    if nz.size:
        cut = float(np.percentile(nz, 95))
        data = (img.get_fdata() >= cut).astype(np.uint8)
img.__class__(data, img.affine, img.header).to_filename(dst)
PYINNER
  then
    log_warn "${sub}: failed to threshold warped SN atlas"
    [[ -z ${ANTPUNK_ROI_DEBUG_DIR:-} ]] && rm -rf "${tmp_dir}"
    return 1
  fi
  local subject_roi_tse="${base_mask_tse}"
  if [[ ${roi_strategy} == "intensity" ]]; then
    local intensity_mask="${out_dir}/roi_intensity_in_tse.nii.gz"
    if build_intensity_roi_mask "${sub}" "${base_mask_tse}" "${intensity_mask}"; then
      subject_roi_tse="${intensity_mask}"
      log_info "${sub}: built intensity-derived ROI (${intensity_top_vox} voxels, sigma=${intensity_smooth_sigma})"
    else
      log_warn "${sub}: intensity ROI build failed; falling back to atlas mask"
    fi
  fi

  local subj_roi_anat_tmp="${tmp_dir}/subject_roi_in_anat"
  if ! 3dAllineate -overwrite -master "${anat_img}" -base "${anat_img}"     -1Dmatrix_apply "${tse_to_anat_mat}" -input "${subject_roi_tse}"     -prefix "${subj_roi_anat_tmp}" -final NN -float >/dev/null 2>&1; then
    log_warn "${sub}: 3dAllineate failed when warping ROI into anat space"
    rm -rf "${tmp_dir}"
    return 1
  fi
  if [[ -f ${subj_roi_anat_tmp}.nii.gz ]]; then
    mv "${subj_roi_anat_tmp}.nii.gz" "${subj_roi_anat}"
  elif [[ -f ${subj_roi_anat_tmp}.nii ]]; then
    mv "${subj_roi_anat_tmp}.nii" "${subj_roi_anat}"
  else
    local afni_dataset=""
    if [[ -f ${subj_roi_anat_tmp}+orig.HEAD ]]; then
      afni_dataset="${subj_roi_anat_tmp}+orig"
    elif [[ -f ${subj_roi_anat_tmp}+tlrc.HEAD ]]; then
      afni_dataset="${subj_roi_anat_tmp}+tlrc"
    fi
    if [[ -n ${afni_dataset} ]]; then
      if ! 3dAFNItoNIFTI -overwrite -prefix "${subj_roi_anat}" "${afni_dataset}" >/dev/null 2>&1; then
        log_warn "${sub}: 3dAFNItoNIFTI failed for subject ROI"
        rm -rf "${tmp_dir}"
        return 1
      fi
      rm -f "${afni_dataset}.HEAD" "${afni_dataset}.BRIK"
    else
      log_warn "${sub}: subject ROI anat file missing after 3dAllineate"
      rm -rf "${tmp_dir}"
      return 1
    fi
  fi

  if (( roi_dilate_vox > 0 )); then
    if command -v 3dmask_tool >/dev/null 2>&1; then
      local dil_tmp="${tmp_dir}/dilated_${sub}_anat.nii.gz"
      if 3dmask_tool -input "${subj_roi_anat}" -dilate_input ${roi_dilate_vox} -prefix "${dil_tmp}" >/dev/null 2>&1; then
        mv "${dil_tmp}" "${subj_roi_anat}"
      else
        log_warn "${sub}: 3dmask_tool dilation failed; keeping original ROI"
      fi
    else
      log_warn "${sub}: 3dmask_tool not found; skipping ROI dilation"
    fi
  fi

  local align_tmp="${tmp_dir}/aligned_${sub}_anat.nii.gz"
  local target_native_vox=""
  local -a target_components=()
  local roi_shift_applied=0
  if [[ -n ${roi_offset_json} ]]; then
    if python3 "${shift_roi_script}" "${subj_roi_anat}" "${align_tmp}"         --offset-json "${roi_offset_json}" --subject "${sub}" >/dev/null 2>&1; then
      mv "${align_tmp}" "${subj_roi_anat}"
      roi_shift_applied=1
      log_info "${sub}: aligned anat ROI using per-subject offsets from ${roi_offset_json}"
    else
      rm -f "${align_tmp}"
      log_warn "${sub}: failed to apply ROI offset from ${roi_offset_json}; keeping native centroid"
    fi
  fi

  if [[ ${roi_shift_applied} -eq 0 && -n ${roi_target_vox} ]]; then
    target_components=(${roi_target_vox//,/ })
    if [[ ${#target_components[@]} -ge 3 ]]; then
      local template_target_mask="${tmp_dir}/template_target_mask_${sub}.nii.gz"
      if python3 - "${template_ref}" "${template_target_mask}"         "${target_components[0]}" "${target_components[1]}" "${target_components[2]}" <<'PYINNER'
import nibabel as nib, numpy as np, sys
tpl=nib.load(sys.argv[1])
out=sys.argv[2]
x=float(sys.argv[3])
y=float(sys.argv[4])
z=float(sys.argv[5])
ijk=[int(round(v)) for v in (x,y,z)]
if any(i<0 or i>=tpl.shape[idx] for idx,i in enumerate(ijk)):
    sys.exit(1)
data=np.zeros(tpl.shape, dtype=np.uint8)
data[tuple(ijk)]=1
nib.Nifti1Image(data, tpl.affine, tpl.header).to_filename(out)
PYINNER
      then
        local template_target_src="${template_target_mask}"
        if command -v 3dmask_tool >/dev/null 2>&1; then
          local template_target_dil="${tmp_dir}/template_target_mask_${sub}_dil.nii.gz"
          if 3dmask_tool -input "${template_target_mask}" -dilate_input 1             -prefix "${template_target_dil}" >/dev/null 2>&1; then
            template_target_src="${template_target_dil}"
          fi
        fi
        local target_point_anat="${tmp_dir}/roi_target_in_anat_${sub}.nii.gz"
        if 3dNwarpApply -overwrite -master "${anat_img}"           -nwarp "INV(${warp_field}) INV(${warp_affine})"           -source "${template_target_src}" -ainterp NN -prefix "${target_point_anat}" >/dev/null 2>&1; then
          if target_native_vox=$(python3 - "${target_point_anat}" <<'PYINNER'
import nibabel as nib, numpy as np, sys
img=nib.load(sys.argv[1])
data=img.get_fdata()
coords=np.argwhere(data>0.5)
if coords.size==0:
    sys.exit(1)
centroid=coords.mean(axis=0)
print(f"{centroid[0]:.6f} {centroid[1]:.6f} {centroid[2]:.6f}")
PYINNER
); then
            :
          else
            target_native_vox=""
          fi
        else
          target_native_vox=""
        fi
      else
        target_native_vox=""
      fi
    fi
  fi
  if [[ ${roi_shift_applied} -eq 0 && -n ${target_native_vox} ]]; then
    local target_native_components=(${target_native_vox})
    if python3 "${shift_roi_script}" "${subj_roi_anat}" "${align_tmp}"       "${target_native_components[0]}" "${target_native_components[1]}" "${target_native_components[2]}" >/dev/null 2>&1; then
      mv "${align_tmp}" "${subj_roi_anat}"
      log_info "${sub}: aligned anat ROI centroid to template target (${roi_target_vox}) via subject coords (${target_native_vox})"
    else
      log_warn "${sub}: failed to align anat ROI centroid to converted target (${target_native_vox})"
      rm -f "${align_tmp}"
    fi
  elif [[ ${roi_shift_applied} -eq 0 && -n ${roi_target_vox} ]]; then
    log_warn "${sub}: unable to convert ROI target (${roi_target_vox}) into subject voxels; keeping native centroid"
  fi

  if ! 3dNwarpApply -overwrite -master "${template_ref}" -nwarp "${warp_field} ${warp_affine}"     -source "${subj_roi_anat}" -ainterp NN -prefix "${subj_roi_template}" >/dev/null 2>&1; then
    log_warn "${sub}: 3dNwarpApply failed when pushing ROI into template space"
    if [[ -z ${ANTPUNK_ROI_DEBUG_DIR:-} ]]; then
      rm -rf "${tmp_dir}"
    else
      log_info "${sub}: ROI debug artifacts left in ${tmp_dir}"
    fi
    return 1
  fi
  if [[ ${disable_roi_crop} -eq 0 && -n ${template_bbox_mask} && -f ${template_bbox_mask} ]]; then
    local subj_roi_cropped="${tmp_dir}/subject_roi_template_cropped.nii.gz"
    if 3dcalc -a "${subj_roi_template}" -b "${template_bbox_mask}" -expr 'step(a)*step(b)'       -prefix "${subj_roi_cropped}" >/dev/null 2>&1; then
      mv "${subj_roi_cropped}" "${subj_roi_template}"
    else
      log_warn "${sub}: failed to crop ROI with template bbox mask"
    fi
  fi
  if (( roi_erode_vox > 0 )); then
    if command -v 3dmask_tool >/dev/null 2>&1; then
      local subj_roi_eroded="${tmp_dir}/subject_roi_template_eroded.nii.gz"
      if 3dmask_tool -input "${subj_roi_template}" -dilate_input -${roi_erode_vox}         -prefix "${subj_roi_eroded}" >/dev/null 2>&1; then
        mv "${subj_roi_eroded}" "${subj_roi_template}"
      else
        log_warn "${sub}: 3dmask_tool erosion failed; keeping uncropped ROI"
      fi
    else
      log_warn "${sub}: 3dmask_tool not found; skipping ROI erosion"
    fi
  fi
  local subj_roi_from_template="${tmp_dir}/subject_roi_in_anat_from_template.nii.gz"
  if 3dNwarpApply -overwrite -master "${anat_img}"     -nwarp "INV(${warp_field}) INV(${warp_affine})"     -source "${subj_roi_template}" -ainterp NN -prefix "${subj_roi_from_template}" >/dev/null 2>&1; then
    mv "${subj_roi_from_template}" "${subj_roi_anat}"
  else
    log_warn "${sub}: failed to map template ROI back into anat space"
  fi

  if [[ -n ${ANTPUNK_ROI_DEBUG_DIR:-} ]]; then
    local debug_dir="${ANTPUNK_ROI_DEBUG_DIR}/tmp_subject_roi_${sub}"
    mkdir -p "${debug_dir}"
    cp -f "${subj_roi_anat}" "${debug_dir}/subject_roi_in_anat_final.nii.gz"
    cp -f "${subj_roi_template}" "${debug_dir}/subject_roi_in_template_final.nii.gz"
    [[ -f ${template_bbox_mask} ]] && cp -f "${template_bbox_mask}" "${debug_dir}/template_bbox_mask.nii.gz"
    cp -f "${anat_img}" "${debug_dir}/anatSS.${sub}.nii" >/dev/null 2>&1 || true
  fi
  if [[ -z ${ANTPUNK_ROI_DEBUG_DIR:-} ]]; then
    rm -rf "${tmp_dir}"
  else
    log_info "${sub}: ROI debug artifacts left in ${tmp_dir}"
  fi
  printf '%s\n%s\n' "${subj_roi_template}" "${subj_roi_anat}"
  return 0
}

warp_sn_to_tse() {
  local sub=$1 mode=$2 out_file=$3
  local tse_vol
  tse_vol=$(find_tse "${sub}") || return 1
  local anat_to_tse=${motip_root}/TSE/${sub}_anat_toTSE_mat.aff12.1D
  [[ -f ${anat_to_tse} ]] || return 1

  if [[ ${mode} == antpunk ]]; then
    local out_dir="${output_root}/${sub}"
    local roi_inv="${out_dir}/antpunk_${sub}_roi_1InverseWarp.nii.gz"
    local aff_mat="${out_dir}/anat_to_template_antpunk_0GenericAffine.mat"
    local anat_img="${sswarper_root}/${sub}/anatQQ.${sub}.nii"
    [[ -f ${roi_inv} && -f ${aff_mat} && -f ${anat_img} ]] || return 1
    local tmp_sn_dir
    tmp_sn_dir=$(mktemp -d "${out_dir}/tmp_snwarp_${mode}_${sub}.XXXX")
    local sn_anat="${tmp_sn_dir}/sn_template_in_anat.nii.gz"
    if ! antsApplyTransforms -d 3 -i "${sn_template}" -r "${anat_img}" \
      -t "${roi_inv}" -t "[${aff_mat},1]" -n Linear -o "${sn_anat}" >/dev/null 2>&1; then
      rm -rf "${tmp_sn_dir}"
      return 1
    fi
    if ! 3dAllineate -overwrite -base "${tse_vol}" -1Dmatrix_apply "${anat_to_tse}" \
      -input "${sn_anat}" -prefix "${out_file}" -final wsinc5 -float >/dev/null 2>&1; then
      rm -rf "${tmp_sn_dir}"
      return 1
    fi
    rm -rf "${tmp_sn_dir}"
    return 0
  fi

  local warp_field="${sswarper_root}/${sub}/anatQQ.${sub}_WARP.nii"
  local warp_mat="${sswarper_root}/${sub}/anatQQ.${sub}.aff12.1D"
  [[ -f ${warp_field} && -f ${warp_mat} ]] || return 1
  local nwarp="INV(${warp_field}) INV(${warp_mat}) ${anat_to_tse}"
  3dNwarpApply -overwrite -master "${tse_vol}" -nwarp "${nwarp}" \
    -source "${sn_template}" -ainterp wsinc5 -prefix "${out_file}" >/dev/null 2>&1 || return 1
  return 0
}

compute_dice() {
  local mask_a=$1 mask_b=$2
  python - "$mask_a" "$mask_b" <<'PY'
import nibabel as nib, numpy as np, sys
a,b=sys.argv[1],sys.argv[2]
try:
    da=nib.load(a).get_fdata()
    db=nib.load(b).get_fdata()
except Exception:
    print('nan')
    sys.exit(0)
ma=da>0
mb=db>0
intersection=np.logical_and(ma,mb).sum()
vol_a=ma.sum()
vol_b=mb.sum()
if vol_a+vol_b==0:
    print('nan')
else:
    dice=(2.0*intersection)/(vol_a+vol_b)
    print(f"{dice:.5f}")
PY
}

mask_voxel_count() {
  local mask_path=$1
  python - "$mask_path" <<'PY'
import nibabel as nib, numpy as np, sys
path=sys.argv[1]
try:
    data=nib.load(path).get_fdata()
except Exception:
    print(0)
    sys.exit(0)
print(int(np.count_nonzero(data)))
PY
}

build_intensity_roi_mask() {
  local sub=$1
  local base_roi=$2
  local out_mask=$3
  local corrected_path
  if ! corrected_path=$(find_corrected_tse "${sub}"); then
    log_warn "${sub}: corrected TSE not found; cannot build intensity ROI"
    return 1
  fi
  if [[ ! -f ${base_roi} ]]; then
    log_warn "${sub}: base ROI ${base_roi} missing; cannot build intensity ROI"
    return 1
  fi
  python3 - "${corrected_path}" "${base_roi}" "${out_mask}" \
    "${intensity_top_vox}" "${intensity_smooth_sigma}" <<'PY'
import math, sys
import nibabel as nib
import numpy as np
tse_path, roi_path, out_path, top_str, sigma_str = sys.argv[1:6]
top_n = max(1, int(float(top_str)))
sigma = float(sigma_str)
tse_img = nib.load(tse_path)
tse_data = np.asarray(tse_img.get_fdata(), dtype=np.float32)
roi_data = nib.load(roi_path).get_fdata()
mask = roi_data > 0.1
if not mask.any():
    raise SystemExit('ROI mask empty')
try:
    from scipy.ndimage import gaussian_filter
except Exception:
    gaussian_filter = None
if sigma > 0 and gaussian_filter is not None:
    blurred = gaussian_filter(tse_data, sigma=sigma)
else:
    blurred = tse_data
masked_vals = blurred[mask]
if masked_vals.size == 0:
    raise SystemExit('No voxels inside ROI to rank')
if masked_vals.size <= top_n:
    thresh = float(np.min(masked_vals))
else:
    idx = np.argpartition(masked_vals, -top_n)[-top_n]
    thresh = float(np.sort(masked_vals)[-top_n])
new_mask = np.zeros_like(tse_data, dtype=np.uint8)
sel = (blurred >= thresh) & mask
if not np.any(sel):
    sel = mask
new_mask[sel] = 1
nib.Nifti1Image(new_mask, tse_img.affine, tse_img.header).to_filename(out_path)
PY
}

if [[ ${roi_only} -eq 1 ]]; then
  for sub in "${subjects[@]}"; do
    log_info "Rebuilding subject ROI masks for ${sub}"
    subj_roi_template="${output_root}/${sub}/subject_roi_in_template.nii.gz"
    subj_roi_anat="${output_root}/${sub}/subject_roi_in_anat.nii.gz"
    rm -f "${subj_roi_template}" "${subj_roi_anat}"
    if ensure_subject_roi_template "${sub}" >/dev/null; then
      log_info "${sub}: wrote ${subj_roi_template} and ${subj_roi_anat}"
    else
      log_warn "${sub}: ROI rebuild failed"
    fi
  done
  log_info "ROI-only mode complete"
  exit 0
fi

manifest_append() {
  local sub=$1 warp=$2 log_file=$3 status=$4 notes=$5
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${sub}" "${mode}" "${roi_strategy}" "${sswarper_root}/${sub}/anatQQ.${sub}_WARP.nii" \
    "${warp}" "${log_file}" "${status}" "${notes}" >>"${manifest_path}"
}

run_antpunk() {
  local sub=$1
  local ss_dir="${sswarper_root}/${sub}"
  ensure_sswarper_subject "${sub}" || return 1
  local anat="${ss_dir}/anatQQ.${sub}.nii"
  local warp_hint="${ss_dir}/anatQQ.${sub}_WARP.nii"
  [[ -f ${warp_hint} ]] || log_warn "${warp_hint} missing; proceeding"
  local out_dir="${output_root}/${sub}"
  mkdir -p "${out_dir}"
  local subject_roi_template=""
  local subject_roi_anat=""
  if [[ ${disable_subject_mask} -eq 0 ]]; then
    local -a roi_paths=()
    if mapfile -t roi_paths < <(ensure_subject_roi_template "${sub}" 2>/dev/null); then
      if [[ ${#roi_paths[@]} -ge 2 ]]; then
        subject_roi_template=${roi_paths[0]}
        subject_roi_anat=${roi_paths[1]}
        local template_vox anat_vox
        template_vox=$(mask_voxel_count "${subject_roi_template}")
        anat_vox=$(mask_voxel_count "${subject_roi_anat}")
        if (( template_vox == 0 || anat_vox == 0 )); then
          log_warn "${sub}: subject ROI mask empty (template_vox=${template_vox}, anat_vox=${anat_vox}); disabling ROI guidance"
          subject_roi_template=""
          subject_roi_anat=""
        else
          log_info "${sub}: subject ROI voxels template=${template_vox}, anat=${anat_vox}"
        fi
      else
        log_warn "${sub}: subject ROI helper returned incomplete paths"
      fi
    else
      log_warn "${sub}: failed to build subject ROI masks"
    fi
  fi
  local tmp_dir
  tmp_dir=$(mktemp -d "${out_dir}/tmp_antpunk_${sub}.XXXX")
  local log_stamp
  if [[ -n ${SLURM_JOB_ID:-} ]]; then
    log_stamp="${SLURM_JOB_ID}_$(date +%Y%m%dT%H%M%S)"
  else
    log_stamp="$(date +%Y%m%dT%H%M%S)"
  fi
  local log_file="${log_root}/${sub}_${log_stamp}.log"
  log_info "${sub}: logging antsRegistration output to ${log_file}"
  local template_input="${template_ref}"
  local anat_roi="${anat}"
  if [[ ${disable_roi_crop} -ne 0 ]]; then
    log_info "ROI cropping disabled; using full template/anat volumes"
  fi
  log_info "antsRegistration: ${sub}"
  {
    printf 'antsRegistration starting at %s\n' "$(date)"
    printf 'Input: %s\n' "${anat}"
  } >"${log_file}"

  local -a cmd=(antsRegistration --dimensionality 3 --float 1 \
    --output "[${tmp_dir}/antpunk_${sub}_,${tmp_dir}/antpunk_${sub}_Warped.nii.gz,${tmp_dir}/antpunk_${sub}_InverseWarped.nii.gz]" \
    --interpolation Linear --winsorize-image-intensities [0.005,0.995] \
    --use-histogram-matching 0 --initial-moving-transform "[${template_ref},${anat},1]" \
    --transform Rigid[0.1] --metric MI[${template_ref},${anat},1,64,Regular,0.25] \
    --convergence [1000x500x250,1e-8,10] --shrink-factors 8x4x2 --smoothing-sigmas 3x2x1vox \
    --transform Affine[0.05] --metric MI[${template_ref},${anat},1,64,Regular,0.25] \
    --convergence [500x250x100,1e-8,10] --shrink-factors 4x2x1 --smoothing-sigmas 2x1x0vox \
    --transform SyN[0.05,3,0] --metric CC[${template_ref},${anat},1,4] \
    --convergence [200x200x100,1e-8,10] --shrink-factors 3x2x1 --smoothing-sigmas 2x1x0vox)
  if [[ ${disable_roi_crop} -eq 0 && -n ${template_bbox_mask} ]]; then
    if [[ -n ${subject_roi_anat} ]]; then
      cmd+=( --x "[${template_bbox_mask},${subject_roi_anat}]" )
    else
      log_warn "${sub}: subject ROI missing; using rectangular template mask only"
      cmd+=( --x "[${template_bbox_mask}]" )
    fi
  elif [[ -n ${subject_roi_template} && -n ${subject_roi_anat} ]]; then
    log_info "${sub}: rectangular ROI mask unavailable; falling back to subject ROI pair"
    cmd+=( --x "[${subject_roi_template},${subject_roi_anat}]" )
  else
    log_info "${sub}: proceeding without antsRegistration ROI mask"
  fi

  if [[ -n ${ants_config} ]]; then
    log_warn "Custom ANTs config ${ants_config} not yet parsed; using built-in stages"
  fi

  {
    printf 'Command: %s\n' "${cmd[*]}"
    "${cmd[@]}"
    printf 'antsRegistration finished at %s\n' "$(date)"
  } >>"${log_file}" 2>&1

  local warp_out="${out_dir}/anat_to_template_antpunk_1Warp.nii.gz"
  local inv_out="${out_dir}/template_to_anat_antpunk_1InverseWarp.nii.gz"
  local aff_out="${out_dir}/anat_to_template_antpunk_0GenericAffine.mat"
  cp "${tmp_dir}/antpunk_${sub}_1Warp.nii.gz" "${out_dir}/antpunk_${sub}_roi_1Warp.nii.gz"
  cp "${tmp_dir}/antpunk_${sub}_1InverseWarp.nii.gz" "${out_dir}/antpunk_${sub}_roi_1InverseWarp.nii.gz"
  embed_warp "${tmp_dir}/antpunk_${sub}_1Warp.nii.gz" "${warp_out}"
  embed_warp "${tmp_dir}/antpunk_${sub}_1InverseWarp.nii.gz" "${inv_out}"
  mv "${tmp_dir}/antpunk_${sub}_0GenericAffine.mat" "${aff_out}"
  if [[ -z ${ANTPUNK_ROI_DEBUG_DIR:-} ]]; then
    rm -rf "${tmp_dir}"
  else
    log_info "${sub}: ROI debug artifacts left in ${tmp_dir}"
  fi

  local components_file="${out_dir}/template_to_tse_components.txt"
  local anat_to_tse=${motip_root}/TSE/${sub}_anat_toTSE_mat.aff12.1D
  if [[ -f ${anat_to_tse} ]]; then
    printf '%s %s\n' "${inv_out}" "${anat_to_tse}" >"${components_file}"
  fi

  printf '%s\n%s' "${warp_out}" "${log_file}"
}

manifest_init

log_info "antpunk_v1: ${#subjects[@]} subjects (mode=${mode}, ROI strategy=${roi_strategy}, ROI bbox pad=${roi_bbox_pad_mm} mm; vox override=${roi_bbox_pad_vox:-auto})"

for sub in "${subjects[@]}"; do
  log_info "Processing ${sub}"
  local_warp=""
  local_log=""
  case ${mode} in
    full|refine_only)
      mapfile -t run_output < <(run_antpunk "${sub}" || true)
      if [[ ${#run_output[@]} -lt 2 ]]; then
        log_warn "${sub}: antpunk run failed"
        continue
      fi
      local_warp=${run_output[0]}
      local_log=${run_output[1]}
      run_qc "${sub}"
      ;;
    qc_only)
      local_warp="${output_root}/${sub}/anat_to_template_antpunk_1Warp.nii.gz"
      [[ -f ${local_warp} ]] || { log_warn "${sub}: antpunk warp missing"; continue; }
      run_qc "${sub}"
      ;;
    *)
      die "Unsupported mode ${mode}"
      ;;
  esac
  status_file="${output_root}/${sub}/qc_status.txt"
  notes_file="${output_root}/${sub}/qc_notes.txt"
  qc_status="PENDING"
  qc_notes=""
  [[ -f ${status_file} ]] && qc_status=$(<"${status_file}")
  [[ -f ${notes_file} ]] && qc_notes=$(<"${notes_file}")
  manifest_append "${sub}" "${local_warp}" "${local_log}" "${qc_status}" "${qc_notes}"
done

log_info "Manifest updated at ${manifest_path}"
