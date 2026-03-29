#!/bin/bash -l
# Analysis_Generator.sh
# Refreshes the low-signal manifest and renders SN_ASR_DRAFT twice—once for the
# Black ANT V1 metrics table and once for the Black ANT V2 table—so we can compare
# outputs side-by-side without running the pipeline multiple times.

set -euo pipefail

R_CMD="${SN_ASR_R_BIN:-R}"

setup_r_environment() {
  local requested_module=${SN_ASR_R_MODULE:-R/4.2.1}
  if ! type module >/dev/null 2>&1; then
    for init_script in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
      if [[ -f "$init_script" ]]; then
        # shellcheck disable=SC1091
        source "$init_script"
        break
      fi
    done
  fi

  if type module >/dev/null 2>&1; then
    local loaded_module=""
    local candidates=("${requested_module}" "R")
    for candidate in "${candidates[@]}"; do
      [[ -z "$candidate" ]] && continue
      echo "Attempting to load module $candidate"
      if module load "$candidate"; then
        loaded_module="$candidate"
        break
      fi
    done
    if [[ -z "$loaded_module" ]]; then
      echo "[WARN] Unable to load requested R module(s); falling back to PATH"
    fi
  else
    echo "[WARN] module command unavailable; assuming R is already on PATH"
  fi

  local pandoc_root=${SN_ASR_PANDOC_DIR:-$PWD/Software/pandoc-3.1.12.3}
  local pandoc_bin="${pandoc_root}/bin"
  if [[ -x "${pandoc_bin}/pandoc" ]]; then
    export PATH="${pandoc_bin}:$PATH"
    export RSTUDIO_PANDOC="${pandoc_bin}"
  else
    echo "[WARN] Pandoc not found under ${pandoc_root}; relying on system installation"
  fi

  export R_LIBS_USER=${R_LIBS_USER:-$HOME/R/local/lib/R/library}
  export PATH=$HOME/R/local/bin:$PATH

  R_CMD=${SN_ASR_R_BIN:-${R_CMD:-R}}
  if ! command -v "$R_CMD" >/dev/null 2>&1; then
    echo "ERROR: R executable not found after module/Pandoc setup. Aborting." >&2
    exit 1
  fi
}

refresh_manifest() {
  echo "Refreshing low-signal manifest..."
  python3 Scripts/refresh_low_signal_manifest.py
}

run_render() {
  local dataset_path="$1"
  local label="$2"

  if [[ ! -f "$dataset_path" ]]; then
    echo "ERROR: dataset for ${label} not found at ${dataset_path}" >&2
    return 1
  fi

  local slug
  slug=$(echo "$label" | tr -cd '[:alnum:]')
  if [[ -z "$slug" ]]; then
    slug="SN"
  fi

  echo "Rendering SN_ASR_DRAFT for ${label} using ${dataset_path}"
  SN_ASR_SN_TABLE="$dataset_path" \
  SN_ASR_SN_LABEL="$label" \
  "$R_CMD" -e "rmarkdown::render('Scripts/SN_ASR_DRAFT_analysis.Rmd', output_dir='MRI_data/analysis', output_file=sprintf('SN_ASR_DRAFT_analysis_%s_%s.html', '${slug}', format(Sys.time(), '%Y%m%d')) )"
}

main() {
  setup_r_environment

  refresh_manifest

  local sn_v1="${SN_ASR_SN_TABLE_V1:-${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data/analysis/SN_nmdata_autosaved.txt}"
  local sn_v2="${SN_ASR_SN_TABLE_V2:-${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data/analysis/SN_nmdata_autosaved_2.txt}"
  local sn_v3="${SN_ASR_SN_TABLE_V3:-${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data/analysis/SN_nmdata_autosaved_20251124_thr80.txt}"

  run_render "$sn_v1" "${SN_ASR_SN_LABEL_V1:-BlackANT_V1}"
  run_render "$sn_v2" "${SN_ASR_SN_LABEL_V2:-BlackANT_V2}"

  if [[ -f "$sn_v3" ]]; then
    run_render "$sn_v3" "${SN_ASR_SN_LABEL_V3:-BlackANT_V3_thr40}"
  else
    echo "Skipping V3 render; dataset not found at ${sn_v3}" >&2
  fi
}

main "$@"
