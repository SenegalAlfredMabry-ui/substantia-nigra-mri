#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: open_roi_variant_freeview.sh <subject_id> [variant ...]
Examples:
  open_roi_variant_freeview.sh 239
  open_roi_variant_freeview.sh 239 vA vC vE
Variants can be specified as variantA/variantB or shorthand (vA, a, C, etc.).
Set FREEVIEW to override the Freeview binary path. Set UNDERLAY to force a custom underlay volume.
USAGE
}

if [[ $# -lt 1 ]]; then
  usage >&2
  exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
VARIANT_ROOT="$REPO_ROOT/MRI_data/WARPS_antpunk/variants"
SUBJECT=$1
shift
SUBJECT_ROOT="$REPO_ROOT/MRI_data/WARPS_antpunk/$SUBJECT"

normalize_variant() {
  local raw=${1^^}
  case "$raw" in
    VARIANTA|VA|A) echo "variantA" ;;
    VARIANTB|VB|B) echo "variantB" ;;
    VARIANTC|VC|C) echo "variantC" ;;
    VARIANTD|VD|D) echo "variantD" ;;
    VARIANTE|VE|E) echo "variantE" ;;
    *) return 1 ;;
  esac
}

variant_letter() {
  local variant=$1
  echo "${variant#variant}"
}

if [[ $# -gt 0 ]]; then
  declare -a requested=()
  for arg in "$@"; do
    if ! norm=$(normalize_variant "$arg"); then
      echo "[ERROR] Unknown variant label '$arg'" >&2
      exit 1
    fi
    requested+=("$norm")
  done
else
  requested=(variantA variantB variantC variantD variantE)
fi

# Remove duplicate entries while preserving order
variants=()
for v in "${requested[@]}"; do
  skip=false
  for existing in "${variants[@]}"; do
    if [[ "$existing" == "$v" ]]; then
      skip=true
      break
    fi
  done
  [[ "$skip" == true ]] || variants+=("$v")
done

FREEVIEW_BIN=${FREEVIEW:-freeview}
if ! command -v "$FREEVIEW_BIN" >/dev/null 2>&1; then
  echo "[ERROR] Could not find Freeview binary '$FREEVIEW_BIN'" >&2
  exit 1
fi

color_for_variant() {
  case "$1" in
    A) echo "1,0.2,0.2" ;;
    B) echo "0.2,0.8,0.2" ;;
    C) echo "1,0.85,0.2" ;;
    D) echo "0.2,0.8,1" ;;
    E) echo "0.85,0.3,0.85" ;;
    *) echo "1,1,1" ;;
  esac
}

resolve_variant_anat() {
  local variant=$1
  local subj_dir="$VARIANT_ROOT/$variant/$SUBJECT"
  shopt -s nullglob
  for candidate in "$subj_dir"/anatQQ.*.nii; do
    echo "$candidate"
    shopt -u nullglob
    return 0
  done
  shopt -u nullglob
  return 1
}

resolve_template_volume() {
  local label=${1,,}
  local path=""
  case "$label" in
    template|antpunk)
      path="$SUBJECT_ROOT/sn_template_inTSE_antpunk.nii.gz"
      ;;
    sswarper|ss)
      path="$SUBJECT_ROOT/sn_template_inTSE_sswarper.nii.gz"
      ;;
    *)
      return 1
      ;;
  esac
  if [[ -f "$path" ]]; then
    echo "$path"
    return 0
  fi
  return 1
}

resolve_underlay() {
  local override=${UNDERLAY:-}
  if [[ -n "$override" ]]; then
    if [[ -f "$override" ]]; then
      echo "$override"
      return 0
    fi
    if norm=$(normalize_variant "$override" 2>/dev/null); then
      if path=$(resolve_variant_anat "$norm" 2>/dev/null); then
        echo "$path"
        return 0
      fi
    fi
    if path=$(resolve_template_volume "$override" 2>/dev/null); then
      echo "$path"
      return 0
    fi
    echo "[ERROR] UNDERLAY '$override' is neither a readable file/variant nor a recognized template alias (template|antpunk|sswarper)." >&2
    return 1
  fi

  for variant in "${variants[@]}"; do
    if path=$(resolve_variant_anat "$variant" 2>/dev/null); then
      echo "$path"
      return 0
    fi
  done

  for label in template sswarper; do
    if path=$(resolve_template_volume "$label" 2>/dev/null); then
      echo "$path"
      return 0
    fi
  done

  echo "[ERROR] Could not find an underlay volume for subject $SUBJECT. Set UNDERLAY to override." >&2
  return 1
}

underlay=$(resolve_underlay) || exit 1
echo "[INFO] Using underlay: $underlay" >&2

declare -a overlays=()
for variant in "${variants[@]}"; do
  letter=$(variant_letter "$variant")
  subj_dir="$VARIANT_ROOT/$variant/$SUBJECT"
  if [[ ! -d "$subj_dir" ]]; then
    echo "[WARN] No directory for $variant / subject $SUBJECT (expected $subj_dir)" >&2
    continue
  fi
  mask="$subj_dir/subject_roi_in_anat_final_v${letter}.nii.gz"
  if [[ ! -f "$mask" ]]; then
    mask="$subj_dir/subject_roi_in_anat_v${letter}.nii.gz"
  fi
  if [[ ! -f "$mask" ]]; then
    echo "[WARN] Missing ROI for $variant / subject $SUBJECT" >&2
    continue
  fi
  color=$(color_for_variant "$letter")
  mask_name=$(basename "$mask")
  overlays+=('-v' "${mask}:name=${mask_name}:color=${color}:opacity=0.45")
  echo "[INFO] Added overlay ($variant): $mask" >&2
  if anat=$(resolve_variant_anat "$variant" 2>/dev/null); then
    anat_name=$(basename "$anat")
    overlays+=('-v' "${anat}:name=${anat_name}:visible=0")
    echo "[INFO] Added reference volume ($variant): $anat" >&2
  fi
done

for label in template sswarper; do
  if path=$(resolve_template_volume "$label" 2>/dev/null); then
    name=$(basename "$path")
    overlays+=('-v' "${path}:name=${name}:visible=0")
    echo "[INFO] Added template reference ($label): $path" >&2
  fi
done

if [[ ${#overlays[@]} -eq 0 ]]; then
  echo "[ERROR] No ROI overlays available for subject $SUBJECT" >&2
  exit 1
fi

echo "[INFO] Launching Freeview with underlay $underlay" >&2
exec "$FREEVIEW_BIN" -v "$underlay" "${overlays[@]}"
