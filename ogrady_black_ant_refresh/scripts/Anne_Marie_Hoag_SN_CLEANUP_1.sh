#!/bin/bash -l
#
# FreeSurfer recon-all wrapper for the 22 cleanup HUMAN_EBR IDs
# Cloned from freesurfer_preprocess copy.sh but scoped to the late cohort so we
# can rerun recon-all only when a segmentation folder is missing.
#
#SBATCH -n 8
#SBATCH -c 4
#SBATCH -J freesurfer_cleanup
#SBATCH -o sbatch-freesurfer-cleanup-%j.out
#SBATCH -e sbatch-freesurfer-cleanup-%j.err

set -euo pipefail

# Public GitHub copy: use an env override or replace these example IDs.
default_subjects=(100 101 102)
if [[ -n ${CLEANUP_SUBJECTS:-} ]]; then
  read -r -a subjects <<<"${CLEANUP_SUBJECTS}"
else
  subjects=("${default_subjects[@]}")
fi

mprage_sub_file=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/locations/MPRAGE_subjects.txt
mprage_loc_file=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/locations/MPRAGE_loc.txt
seg_root=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/segmentation

echo "SLURM_JOB_ID = ${SLURM_JOB_ID:-unknown}"
echo "Cleanup subjects: ${subjects[*]}"

module load afni
FREESURFER_HOME=${FREESURFER_HOME:-/path/to/freesurfer}
source "$FREESURFER_HOME/SetUpFreeSurfer.sh"

# Build associative maps from the master subject/loc files so we can reuse the
# same lookup logic as the original script.
if [[ ! -f $mprage_sub_file || ! -f $mprage_loc_file ]]; then
  echo "Missing MPRAGE subject/location files; aborting." >&2
  exit 1
fi

declare -A loc_map=()
numsub=$(wc -l < "$mprage_sub_file")
for i in $(seq 1 "$numsub"); do
  sub_folder=$(sed -n "${i}p" "$mprage_sub_file")
  sub=${sub_folder: -3}
  loc=$(sed -n "${i}p" "$mprage_loc_file")
  loc_map[$sub]=$loc
  done

find_mprage() {
  local sub=$1
  local loc=${loc_map[$sub]:-}
  [[ -n $loc ]] || return 1
  local base="${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/HUMAN_EBR-${sub}/${loc}"
  local candidates=(
    "$base/mprage.nii"
    "$base/MPRAGE.nii"
    "$base/mprage_scan_0001.nii"
  )
  for cand in "${candidates[@]}"; do
    [[ -f $cand ]] && { echo "$cand"; return 0; }
  done
  return 1
}

for sub in "${subjects[@]}"; do
  printf '\n== HUMAN_EBR-%s ==\n' "$sub"
  if [[ -d ${seg_root}/${sub} ]]; then
    echo "  -> Segmentation already exists; skipping recon-all."
    continue
  fi
  mprage=$(find_mprage "$sub") || {
    echo "  !! No MPRAGE located for ${sub}; skipping." >&2
    continue
  }
  echo "  MPRAGE: $mprage"
  echo "  Deobliquing..."
  3drefit -deoblique "$mprage"
  echo "  Running recon-all..."
  recon-all -s "$sub" -openmp 4 -i "$mprage" -all
  echo "  Copying outputs into ${seg_root}/${sub}"
  cp -r "$FREESURFER_HOME/subjects/${sub}" "$seg_root/" || true
  if [[ ! -d ${seg_root}/${sub} ]]; then
    echo "  !! recon-all output missing for ${sub}" >&2
  fi
done

# Post-run: ensure aparc.a2009s+aseg gets converted to NIfTI where missing
for sub in "${subjects[@]}"; do
  sub_dir="${seg_root}/${sub}"
  [[ -d $sub_dir ]] || continue
  mgz="${sub_dir}/mri/aparc.a2009s+aseg.mgz"
  nii="${sub_dir}/mri/aparc.a2009s+aseg.nii"
  if [[ -f $mgz && ! -f $nii ]]; then
    echo "Converting aparc for ${sub}"
    mri_convert "$mgz" "$nii"
  fi
done

echo "Cleanup FreeSurfer batch complete."
