#!/bin/bash -l
#
# Cleanup-targeted SegmentAAN wrapper: generate arousalNetworkLabels.v10.mgz for
# the 22 late HUMAN_EBR IDs only when the FreeSurfer outputs exist and labels are
# missing. Based on segment_brainstem copy.sh.
#
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -J segmentaan_cleanup
#SBATCH -o sbatch-segmentaan-cleanup-%j.out
#SBATCH -e sbatch-segmentaan-cleanup-%j.err

set -euo pipefail

# Public GitHub copy: use an env override or replace these example IDs.
default_subjects=(100 101 102)
if [[ -n ${CLEANUP_SUBJECTS:-} ]]; then
  read -r -a subjects <<<"${CLEANUP_SUBJECTS}"
else
  subjects=("${default_subjects[@]}")
fi
seg_root=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/segmentation

echo "SLURM_JOB_ID = ${SLURM_JOB_ID:-unknown}"
echo "Cleanup SegmentAAN subjects: ${subjects[*]}"

FREESURFER_HOME=${FREESURFER_HOME:-/path/to/freesurfer}
source "$FREESURFER_HOME/SetUpFreeSurfer.sh"
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=8

for sub in "${subjects[@]}"; do
  printf '\n== HUMAN_EBR-%s ==\n' "$sub"
  subj_dir="${seg_root}/${sub}"
  if [[ ! -d $subj_dir ]]; then
    echo "  !! Missing FreeSurfer directory; run recon-all first."
    continue
  fi
  labels="${subj_dir}/mri/arousalNetworkLabels.v10.mgz"
  if [[ -f $labels ]]; then
    echo "  -> Arousal labels already exist; skipping."
    continue
  fi
  echo "  -> Running SegmentAAN.sh ${sub} ${seg_root}"
  SegmentAAN.sh "$sub" "$seg_root"
done

echo "SegmentAAN cleanup batch complete."
