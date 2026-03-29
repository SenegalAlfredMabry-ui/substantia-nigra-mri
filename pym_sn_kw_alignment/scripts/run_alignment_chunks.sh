#!/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: $0 <type> [chunk-size] [mode]
  type       = highhat | crash
  chunk-size = subjects per chunk (default 3)
  mode       = serial | parallel (default serial)
Example: ./Scripts/run_alignment_chunks.sh highhat 3 parallel
USAGE
}

if (( $# < 1 )); then
  usage
  exit 1
fi

type=$1
chunk_size=${2:-3}
mode=${3:-serial}

case $type in
  highhat)
    roster=tmp/highhat_with_kw_subjects.txt
    prefix=highhat
    ;;
  crash)
    roster=tmp/crash_offsets_subjects.txt
    prefix=crash
    ;;
  *)
    echo "Unknown type '$type'" >&2
    usage
    exit 1
    ;;
esac

if [[ ! -f $roster ]]; then
  echo "Roster file '$roster' not found" >&2
  exit 1
fi

subjects=()
while read -r line; do
  [[ -z $line ]] && continue
  subjects+=($line)
done < "$roster"

total=${#subjects[@]}
if (( total == 0 )); then
  echo "No subjects listed in $roster" >&2
  exit 1
fi

last_job=""
chunk=0
for ((i=0; i<total; i+=chunk_size)); do
  chunk_id=$((chunk+1))
  chunk_file=tmp/${prefix}_chunk${chunk_id}_subjects.txt
  chunk_output=tmp/kw_masks_in_anat_${prefix}_chunk${chunk_id}
  baseline_dir=tmp/${prefix}_chunk${chunk_id}_baseline
  summary=tmp/${prefix}_chunk${chunk_id}_alignment.tsv

  printf "%s\n" "${subjects[@]:i:chunk_size}" > "$chunk_file"

  jobname="align_${prefix}_chunk${chunk_id}"
  log=tmp/${jobname}.log
  dep_option=""
  if [[ $mode == "serial" && -n $last_job ]]; then
    dep_option="--dependency=afterok:$last_job"
  fi

  jobid=$(sbatch --parsable --job-name="$jobname" --output="$log" $dep_option \
    --wrap="${PWD}/Scripts/run_alignment_chunk_job.sh '$chunk_file' '$chunk_output' '$baseline_dir' '$summary'" )

  echo "Submitted chunk $chunk_id ($jobname) -> job $jobid"
  last_job=$jobid
  chunk=$((chunk+1))
done

echo "All chunks submitted; monitor tmp/${prefix}_chunk*_alignment.tsv for results."
