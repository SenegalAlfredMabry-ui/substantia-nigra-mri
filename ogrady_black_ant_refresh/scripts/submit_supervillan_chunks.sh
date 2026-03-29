#!/bin/bash -l
# submit_supervillan_chunks.sh
# Split the SNHEART roster into chunks of outstanding subjects and queue
# supervillan_wallet.sh runs with 2-node requests per chunk.
set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
cd "$repo_root"

roster_file=$repo_root/MRI_data/locations/TSE_SN_ASR_SUBLIST.txt
if [[ ! -f $roster_file ]]; then
  roster_file=$repo_root/MRI_data/locations/TSE_subjects.txt
fi
manifest=$repo_root/MRI_data/WARPS_supervillain/warp_manifest.tsv

if [[ ! -f $roster_file ]]; then
  echo "ERROR: no roster file found" >&2
  exit 1
fi

chunk_size=${SUPERVILLAIN_CHUNK_SIZE:-30}
partition=${SUPERVILLAIN_PARTITION:-cloud}
nodes=${SUPERVILLAIN_NODES:-2}
ntasks=${SUPERVILLAIN_TASKS:-2}
cpus=${SUPERVILLAIN_CPUS:-4}
submit_flag=${SUPERVILLAIN_SUBMIT:-0}
max_chunks=${SUPERVILLAIN_MAX_CHUNKS:-0}

# Map of completed subjects from the manifest.
declare -A done_subs=()
if [[ -f $manifest ]]; then
  tail -n +2 "$manifest" | awk -F'\t' '{print $1}' | sort -u | while read -r sid; do
    [[ -z $sid ]] && continue
    sid=$(printf '%03d' "$sid")
    done_subs[$sid]=1
  done || true
fi

# Build list of outstanding subjects.
remaining=()
while read -r line; do
  line=${line%%#*}
  line=$(echo "$line" | xargs || true)
  [[ -z $line ]] && continue
  if [[ $line =~ ([0-9]+)$ ]]; then
    sid=$(printf '%03d' "${BASH_REMATCH[1]}")
    if [[ -z ${done_subs[$sid]:-} ]]; then
      remaining+=("$sid")
    fi
  fi
done < "$roster_file"

if [[ ${#remaining[@]} -eq 0 ]]; then
  echo "[INFO] No outstanding subjects detected; manifest already covers roster."
  exit 0
fi

echo "[INFO] Outstanding subjects (${#remaining[@]} total): ${remaining[*]}"

chunk_id=0
for ((idx=0; idx<${#remaining[@]}; idx+=chunk_size)); do
  ((chunk_id++))
  if [[ ${max_chunks} -gt 0 && ${chunk_id} > ${max_chunks} ]]; then
    echo "[INFO] Reached SUPERVILLAIN_MAX_CHUNKS=${max_chunks}; stopping after chunk $((chunk_id-1))."
    break
  fi
  chunk=("${remaining[@]:idx:chunk_size}")
  chunk_str=$(printf '%s ' "${chunk[@]}")
  chunk_str=${chunk_str%% }
  job_name=${SUPERVILLAIN_JOB_PREFIX:-svw}_chunk${chunk_id}
  cmd=(sbatch -N "$nodes" -n "$ntasks" -c "$cpus" -p "$partition" -J "$job_name" Scripts/supervillan_wallet.sh)
  if [[ $submit_flag -eq 1 ]]; then
    echo "[INFO] Submitting $job_name for subjects: $chunk_str"
    SUPERVILLAIN_SUBJECTS="$chunk_str" "${cmd[@]}"
  else
    echo "[DRY-RUN] SUPERVILLAIN_SUBJECTS=\"$chunk_str\" ${cmd[*]}"
  fi
done

echo "[INFO] Set SUPERVILLAIN_SUBMIT=1 to actually queue the jobs."
