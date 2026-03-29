#!/bin/bash -l
#SBATCH -J irredeemable_v12_kwsn
#SBATCH -p cloud
#SBATCH -c 4
#SBATCH -o sbatch-%j.out
#SBATCH -e sbatch-%j.err

set -euo pipefail

if command -v module >/dev/null 2>&1; then
  module load afni >/dev/null 2>&1 || true
fi
if ! command -v 3dcalc >/dev/null 2>&1; then
  echo "[ERROR] AFNI/3dcalc not available on PATH; did module load afni run?" >&2
  exit 1
fi

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
shared_root=${SNHEART_ROOT:-/path/to/SNHEART}
cd "$repo_root"

MATLAB_ROOT=${MATLAB_ROOT:-"$repo_root/Software"}
if [[ ! -x "$MATLAB_ROOT/bin/matlab" ]]; then
  echo "[ERROR] MATLAB binary not found at $MATLAB_ROOT/bin/matlab" >&2
  exit 1
fi

run_dir=${KWBACKUP_RUN_DIR:-$(ls -dt "${shared_root}"/tmp/raw_kw_loop_* 2>/dev/null | head -n 1 || true)}
if [[ -z ${run_dir} || ! -d ${run_dir} ]]; then
  echo "[ERROR] Could not resolve KWBACKUP_RUN_DIR" >&2
  exit 1
fi

resolve_raw_kw_subjects() {
  local hh=${run_dir}/highhat/highhat_rawkw_subjects.txt
  local cr=${run_dir}/crash/crash_rawkw_subjects.txt
  [[ -f ${hh} && -f ${cr} ]] || return 1
  awk 'NF{print $1}' "${hh}" "${cr}" | sort -u -n | paste -sd' ' -
}

subjects=${IRREDEEMABLE_SUBJECTS:-}
if [[ -z ${subjects} ]]; then
  subjects=$(resolve_raw_kw_subjects || true)
fi
if [[ -z ${subjects} ]]; then
  echo "[ERROR] Could not resolve raw-KW subject list" >&2
  exit 1
fi

kw_backup_root=${KWBACKUP_OUTPUT_ROOT:-${repo_root}/MRI_data/TSE/probabilistic_masks_kwsn_v12}
analysis_root=${IRREDEEMABLE_ANALYSIS_ROOT:-${repo_root}/MRI_data/analysis/tests/v12_kwsn}
brainstem_dir=${IRREDEEMABLE_BRAINSTEM_DIR:-${repo_root}/MRI_data/TSE/brainstem_masks_v10}
chunk_tag=${IRREDEEMABLE_CHUNK_TAG:-}
subject_expect=${IRREDEEMABLE_EXPECT_COUNT:-$(wc -w <<<"${subjects}")}
mkdir -p "$analysis_root"

export KWBACKUP_SOURCE_ROOT=${KWBACKUP_SOURCE_ROOT:-${repo_root}/MRI_data/TSE/probabilistic_masks_v10}
export KWBACKUP_OUTPUT_ROOT="${kw_backup_root}"
export KWBACKUP_PT_VARIANTS="${KWBACKUP_PT_VARIANTS:-KW}"
export KWBACKUP_SUBJECTS="${subjects}"

"${repo_root}/Scripts/prepare_irredeemable_kwsn_mask_root.sh"

echo "[INFO] Using KW-backed mask root ${kw_backup_root}"
echo "[INFO] Subjects: $(wc -w <<<"${subjects}")"
if [[ -n ${chunk_tag} ]]; then
  echo "[INFO] Chunk tag: ${chunk_tag}"
fi

if [[ -n ${IRREDEEMABLE_PT_VARIANTS:-} ]]; then
  read -r -a pt_variants <<<"${IRREDEEMABLE_PT_VARIANTS}"
else
  pt_variants=(KW)
fi
variant_dirs=()
for variant in "${pt_variants[@]}"; do
  label=${variant^^}
  dir="${kw_backup_root}/PT_${label}"
  if [[ -d $dir ]]; then
    variant_dirs+=("$dir")
  else
    echo "[WARN] PT variant $label missing under $kw_backup_root; skipping" >&2
  fi
done
if [[ ${#variant_dirs[@]} -eq 0 ]]; then
  echo "[ERROR] No PT variant directories available under $kw_backup_root" >&2
  exit 1
fi

if [[ -n ${IRREDEEMABLE_CORR_MODES:-} ]]; then
  read -r -a corr_modes <<<"${IRREDEEMABLE_CORR_MODES}"
else
  corr_modes=(legacy)
fi
watch_list=(142 503 505 522 548 602 614)
prob_slug=15

run_idx=0
for variant_dir in "${variant_dirs[@]}"; do
  variant_label=${variant_dir##*/PT_}
  [[ -z ${variant_label} ]] && continue
  for corr_mode in "${corr_modes[@]}"; do
    case ${corr_mode} in
      legacy)
        use_corrected=1
        corr_suffix=""
        corr_tag="legacy"
        ;;
      raw)
        use_corrected=0
        corr_suffix=""
        corr_tag="raw"
        ;;
      *)
        echo "[WARN] Unknown corr mode ${corr_mode}; skipping" >&2
        continue
        ;;
    esac

    clamp_flag=1
    clamp_tag="clamp"
    run_label=${variant_label}_${corr_tag}_${clamp_tag}
    if [[ -n ${chunk_tag} ]]; then
      run_label=${run_label}_${chunk_tag}
    fi
    run_idx=$((run_idx + 1))
    analysis_dir="${analysis_root}/${variant_label}/${corr_tag}/${clamp_tag}"
    mkdir -p "$analysis_dir/logs"

    export IRREDEEMABLE_SUBJECTS="${subjects}"
    export IRREDEEMABLE_MASK_DIR="${variant_dir}"
    export IRREDEEMABLE_ANALYSIS_DIR="${analysis_dir}"
    export IRREDEEMABLE_OUTPUT_BASENAME="SN_nmdata_irredeemable_v12_kwsn_${run_label}"
    export IRREDEEMABLE_OUTPUT_TAG="kwsn_${run_label}"
    export IRREDEEMABLE_PT_CLAMP_SN_SLICES="$clamp_flag"
    export IRREDEEMABLE_USE_CORRECTED="$use_corrected"
    export IRREDEEMABLE_SN_PROB_FIXED=0.15
    export IRREDEEMABLE_ENABLE_NIGROSOMES=0
    export IRREDEEMABLE_PT_USE_SUPERIOR=0
    export IRREDEEMABLE_LOG="$analysis_dir/logs/irredeemable_v12_kwsn_${run_label}.log"
    if [[ $use_corrected -eq 1 && -n ${corr_suffix} ]]; then
      export IRREDEEMABLE_CORR_SUFFIX="$corr_suffix"
    else
      unset IRREDEEMABLE_CORR_SUFFIX
    fi

    echo "[INFO] (${run_idx}) Variant=${variant_label} Corr=${corr_tag} Clamp=${clamp_tag}"
    "$MATLAB_ROOT/bin/matlab" -nodisplay -nosplash -nodesktop \
      -r "addpath('$repo_root/Scripts'); run('Scripts/irredeemable_black_ant_V9.m'); exit;"

    output_file="${analysis_dir}/SN_nmdata_irredeemable_v12_kwsn_${run_label}_prob${prob_slug}_kwsn_${run_label}.txt"
    if [[ ! -f $output_file ]]; then
      echo "[WARN] Expected output $output_file missing" >&2
      continue
    fi

    log_path="${analysis_dir}/logs/check_sn_metrics_${run_label}.log"
    python3 Scripts/check_sn_metrics.py \
      --main "$output_file" \
      --expect "${subject_expect}" \
      --watch "${watch_list[@]}" \
      --check-overlap \
      --mask-dir "$variant_dir" \
      --brainstem-dir "$brainstem_dir" \
      2>&1 | tee "$log_path"
    echo "[INFO] Completed ${run_label}"
  done
done

unset IRREDEEMABLE_SUBJECTS IRREDEEMABLE_MASK_DIR IRREDEEMABLE_ANALYSIS_DIR \
  IRREDEEMABLE_OUTPUT_BASENAME IRREDEEMABLE_OUTPUT_TAG IRREDEEMABLE_PT_CLAMP_SN_SLICES \
  IRREDEEMABLE_USE_CORRECTED IRREDEEMABLE_SN_PROB_FIXED IRREDEEMABLE_ENABLE_NIGROSOMES \
  IRREDEEMABLE_PT_USE_SUPERIOR IRREDEEMABLE_CORR_SUFFIX IRREDEEMABLE_LOG \
  IRREDEEMABLE_CHUNK_TAG IRREDEEMABLE_EXPECT_COUNT
