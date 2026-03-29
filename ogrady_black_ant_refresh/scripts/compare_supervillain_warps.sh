#!/usr/bin/env bash
# Compare refreshed SSwarper outputs against legacy warps and emit displacement stats.
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: compare_supervillain_warps.sh [-m manifest] [-n new_root] [-o old_root] [-t output]
  -m  Path to warp_manifest.tsv (default MRI_data/WARPS_supervillain/warp_manifest.tsv)
  -n  Root directory holding refreshed warps (default MRI_data/WARPS_supervillain)
  -o  Root directory with legacy warps (default ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/WARPS)
  -t  Output TSV for stats (default warp_diff_stats.tsv)
  -h  Show this help
Environment:
  AFNI_BIN (optional) points to directory containing 3dcalc/3dBrickStat.
USAGE
}

manifest="MRI_data/WARPS_supervillain/warp_manifest.tsv"
new_root="MRI_data/WARPS_supervillain"
old_root="${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/WARPS"
output_tsv="warp_diff_stats.tsv"

while getopts "m:n:o:t:h" opt; do
  case "$opt" in
    m) manifest="$OPTARG" ;;
    n) new_root="$OPTARG" ;;
    o) old_root="$OPTARG" ;;
    t) output_tsv="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

AFNI_BIN=${AFNI_BIN:-afni}
if [[ ! -x "$AFNI_BIN/3dcalc" || ! -x "$AFNI_BIN/3dBrickStat" ]]; then
  echo "ERROR: AFNI tools not found under $AFNI_BIN" >&2
  exit 1
fi

if [[ ! -f "$manifest" ]]; then
  echo "ERROR: manifest $manifest not found" >&2
  exit 1
fi

printf 'subject\tmean_mm\tmax_mm\n' > "$output_tsv"
expr_str='sqrt((a-b)^2+(c-d)^2+(e-f)^2)'

tail -n +2 "$manifest" | cut -f1 | while read -r subj; do
  [[ -z "$subj" ]] && continue
  new_warp="$new_root/$subj/anatQQ.${subj}_WARP.nii"
  old_warp="$old_root/$subj/anatQQ.${subj}_WARP.nii"
  if [[ ! -f "$new_warp" || ! -f "$old_warp" ]]; then
    printf '%s\tNA\tNA\n' "$subj" >> "$output_tsv"
    echo "WARN: missing warp for $subj (new: $new_warp, old: $old_warp)" >&2
    continue
  fi
  tmp=$(mktemp /tmp/warp_diff_${subj}_XXXX)
  tmp_nii="${tmp}.nii.gz"
  "$AFNI_BIN/3dcalc" -overwrite \
    -a "$new_warp"'[0]' -b "$old_warp"'[0]' \
    -c "$new_warp"'[1]' -d "$old_warp"'[1]' \
    -e "$new_warp"'[2]' -f "$old_warp"'[2]' \
    -expr "$expr_str" -prefix "$tmp_nii" >/dev/null
  mean=$("$AFNI_BIN/3dBrickStat" -mean "$tmp_nii")
  max=$("$AFNI_BIN/3dBrickStat" -max "$tmp_nii")
  printf '%s\t%s\t%s\n' "$subj" "$mean" "$max" >> "$output_tsv"
  rm -f "$tmp_nii"
done

echo "Wrote $output_tsv"
