#!/usr/bin/env bash
# run_black_ant_v5_pipeline.sh
# -----------------------------------------------------------------------------
# Wrapper that submits the v5 mask builder to SLURM, waits for completion,
# checks its final state via sacct, and if successful launches the
# Irredeemable V5 metrics job. Mirrors the manual cadence described in
# ANT_FARM_RUNS but automates the monitoring.
# -----------------------------------------------------------------------------
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_SCRIPT="${ROOT_DIR}/Scripts/run_build_black_ant_masks_v5.sh"
IRRED_SCRIPT="${ROOT_DIR}/Scripts/run_irredeemable_black_ant_V5.sh"
POLL_INTERVAL=60

submit_job() {
  local script=$1
  if [[ ! -x ${script} ]]; then
    echo "[ERROR] ${script} not found or not executable" >&2
    exit 1
  fi
  sbatch --parsable "${script}"
}

wait_for_completion() {
  local job_id=$1
  echo "[INFO] Waiting for job ${job_id} to finish (poll every ${POLL_INTERVAL}s)..."
  while true; do
    sleep "${POLL_INTERVAL}"
    local state
    state=$(sacct -X -j "${job_id}" --format=State --noheader 2>/dev/null | head -n1 | tr -d ' ')
    if [[ -z ${state} ]]; then
      echo "[INFO] sacct data not ready for job ${job_id}; continuing to poll."
      continue
    fi
    case ${state} in
      COMPLETED)
        echo "[INFO] Job ${job_id} completed successfully."
        return 0
        ;;
      FAILED|CANCELLED|TIMEOUT|PREEMPTED|OUT_OF_MEMORY)
        echo "[ERROR] Job ${job_id} ended with state: ${state}" >&2
        return 1
        ;;
      *)
        echo "[INFO] Job ${job_id} state=${state}; rechecking soon."
        ;;
    esac
  done
}

echo "[INFO] Submitting build_black_ant_masks_v5..."
build_job=$(submit_job "${BUILD_SCRIPT}")
echo "[INFO] build job id: ${build_job}"

if wait_for_completion "${build_job}"; then
  echo "[INFO] Launching Irredeemable V5 now that masks are ready..."
  irred_job=$(submit_job "${IRRED_SCRIPT}")
  echo "[INFO] Irredeemable job id: ${irred_job}"
else
  echo "[ERROR] Skipping Irredeemable launch because build job failed." >&2
  exit 1
fi

echo "[INFO] Pipeline submission complete. Monitor ${build_job} and ${irred_job} with sacct/squeue for details."
