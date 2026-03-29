#!/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: $0 <type> [chunk-size] [--nodes node1,node2,...] [--dedicated]
  type       = highhat | crash
  chunk-size = subjects per chunk (default 3)
  --nodes    = comma-delimited list of nodes to rotate through (default: c-0..c-6)
  --dedicated=reserve a specific node for each chunk (adds --nodelist/--exclusive)
Example: ./Scripts/submit_antpunk_chunks.sh crash 3 --nodes cluster-node-7,cluster-node-8 --dedicated
USAGE
}

if (( $# < 1 )); then
  usage
  exit 1
fi

type=$1
shift
chunk_size=3
if (( $# > 0 )) && [[ $1 != --* ]]; then
  chunk_size=$1
  shift
fi

nodes_arg=""
dedicated=false
while (( $# > 0 )); do
  case "$1" in
    --nodes)
      shift
      nodes_arg="$1"
      shift
      ;;
    --dedicated)
      dedicated=true
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

case $type in
  highhat)
    roster=tmp/highhat_with_kw_subjects.txt
    json=tmp/highhat_roi_offsets_hybrid.json
    prefix=highhat
    ;;
  crash)
    roster=tmp/crash_offsets_subjects.txt
    json=tmp/crash_roi_offsets_hybrid.json
    prefix=crash
    ;;
  *)
    echo "Unknown type '$type'" >&2
    usage
    exit 1
    ;;
esac

if [[ ! -f $roster || ! -f $json ]]; then
  echo "Missing roster ($roster) or offset JSON ($json)" >&2
  exit 1
fi

declare -a subjects
while read -r line; do
  [[ -z $line ]] && continue
  subjects+=("$line")
done < "$roster"

total=${#subjects[@]}
if (( total == 0 )); then
  echo "No subjects listed in $roster" >&2
  exit 1
fi

if [[ -n "$nodes_arg" ]]; then
  IFS=',' read -r -a nodes <<< "$nodes_arg"
else
  nodes=(cluster-node-0 cluster-node-1 cluster-node-2 cluster-node-3 cluster-node-4 cluster-node-5 cluster-node-6)
fi

chunk_count=$(((total + chunk_size - 1) / chunk_size))
if [[ "$dedicated" == true && ${#nodes[@]} -lt chunk_count ]]; then
  echo "Not enough nodes (${#nodes[@]}) for $chunk_count chunks when using --dedicated" >&2
  exit 1
fi

for ((chunk_id=1; chunk_id<=chunk_count; chunk_id++)); do
  start=$(( (chunk_id - 1) * chunk_size ))
  chunk_subjects="${subjects[@]:start:chunk_size}"
  jobname="antpunk_${prefix}_chunk${chunk_id}"
  log=tmp/${jobname}.log
  cmd="export ANTPUNK_ROI_OFFSET_JSON=$json; export ANTPUNK_DISABLE_SUBJECT_MASK=0; export ANTPUNK_DISABLE_ROI_CROP=1; export ANTPUNK_SUBJECTS='$chunk_subjects'; cd ${SNHEART_ROOT:-/path/to/SNHEART}; Scripts/antpunk_v1.sh --mode full"
  sbatch_args=(--parsable --job-name="$jobname" --output="$log")
  if [[ "$dedicated" == true ]]; then
    node="${nodes[chunk_id-1]}"
    sbatch_args+=(--nodes=1 --nodelist="$node" --exclusive)
    extra_info=" node=$node"
  else
    index=$(( (chunk_id - 1) % ${#nodes[@]} ))
    extra_info=" node=${nodes[index]}"
  fi
  sbatch "${sbatch_args[@]}" --wrap="$cmd"
  echo "Submitted $jobname -> $(date)$extra_info subjects=$chunk_subjects"
done
