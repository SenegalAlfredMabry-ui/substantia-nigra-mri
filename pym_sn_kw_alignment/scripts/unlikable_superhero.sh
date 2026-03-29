#!/bin/bash -l
# unlikable_superhero.sh
# -----------------------------------------------------------------------------
# PURPOSE
#   Regenerate the anat->TSE affine matrices (`${sub}_anat_toTSE_mat.aff12.1D`)
#   for the MOTIP IDs currently blocked from build_black_ant_masks.sh. This is
#   step 2 of the legacy neuromelanin pipeline but scoped to a curated roster
#   so we only rebuild what's missing.
# OUTPUTS
#   - Writes `${motip_root}/MRI_data/TSE/<ID>_anat_toTSE_mat.aff12.1D` per
#     subject. That is the exact filename the mask builder and
#     S.H.I.E.L.D_SN_Warp_v1.sh look for when they project template masks back
#     into raw TSE space.
# DEPENDENCIES
#   - AFNI (align_epi_anat.py + 3drefit) on PATH.
#   - MOTIP basic_anat tree plus raw TSE volumes. Script auto-decompresses the
#     `.nii.gz` when only gzip copies exist and mirrors the legacy parameter
#     split (<513 uses ginormous_move, ≥513 uses align_centers+giant_move).
# RELIES ON
#   - build_black_ant_masks.sh, irredeemable_black_ant_V1.m, and downstream
#     MATLAB/R runs. Without these affines, the probabilistic masks cannot land
#     in subject TSE space, so rerun this wrapper whenever those jobs warn about
#     missing `${sub}_anat_toTSE_mat.aff12.1D`.
# -----------------------------------------------------------------------------
#SBATCH -J unlikable_superhero
#SBATCH -p cloud
#SBATCH -c 4
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o sbatch-wrap-align-%j.out
#SBATCH -e sbatch-wrap-align-%j.err
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=FAIL

set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
cd "$repo_root"

module load afni

motip_root=${MOTIP_ROOT:-/path/to/MOTIP2018}
base_dir=${motip_root}/MRI_data
basic_anat_dir=${base_dir}/basic_anat
raw_tse_dir=${base_dir}/TSE
process_dir=${raw_tse_dir}/process

# Public GitHub copy: the real late-cohort roster was replaced with examples.
default_subjects=(100 101 102)

if [[ -n ${UNLIKABLE_SUBJECTS:-} ]]; then
  read -r -a subjects <<<"${UNLIKABLE_SUBJECTS}"
  echo "[INFO] Using UNLIKABLE_SUBJECTS override: ${subjects[*]}"
else
  subjects=("${default_subjects[@]}")
fi

echo "[INFO] Refreshing anat->TSE affines for: ${subjects[*]}"

mkdir -p "${process_dir}"

for sub in "${subjects[@]}"; do
  printf '\n=== HUMAN_EBR-%03d ===\n' "${sub}"
  matrix_out=${raw_tse_dir}/${sub}_anat_toTSE_mat.aff12.1D
  if [[ -f ${matrix_out} ]]; then
    echo "  -> ${matrix_out} already exists; skipping"
    continue
  fi

  anat_base=${basic_anat_dir}/HUMAN_EBR-${sub}/anat.nii
  raw_tse=${raw_tse_dir}/${sub}_raw_tse.nii
  raw_tse_gz=${raw_tse_dir}/${sub}_raw_tse.nii.gz

  if [[ -f ${raw_tse_gz} && ! -f ${raw_tse} ]]; then
    echo "  -> Decompressing ${raw_tse_gz}"
    gunzip -c "${raw_tse_gz}" >"${raw_tse}"
  fi

  if [[ ! -f ${raw_tse} ]]; then
    echo "  !! Missing raw TSE ${raw_tse}; skipping"
    continue
  fi
  if [[ ! -f ${anat_base} ]]; then
    echo "  !! Missing anatomy ${anat_base}; skipping"
    continue
  fi

  echo "  -> Deobliquing ${raw_tse}"
  3drefit -deoblique "${raw_tse}"

  align_args=(
    -overwrite
    -output_dir "${process_dir}"
    -dset1 "${anat_base}"
    -dset2 "${raw_tse}"
    -master_dset1 "${anat_base}"
    -dset1to2
    -dset1_strip None
    -dset2_strip 3dAutomask
    -suffix "${sub}_anat_toTSE"
    -cost lpa)

  if (( sub < 513 )); then
    align_args+=(-ginormous_move -deoblique off)
  else
    align_args+=(-align_centers yes -giant_move)
  fi

  echo "  -> Running align_epi_anat.py"
  align_epi_anat.py "${align_args[@]}"

done

# Move results back into TSE directory and clean temp AFNI outputs
shopt -s nullglob
for aff in "${process_dir}"/*_anat_toTSE_mat.aff12.1D; do
  base=$(basename "${aff}")
  canonical="${base#anat}"
  if [[ ! ${canonical} =~ ^[0-9]+_anat_toTSE_mat\.aff12\.1D$ ]]; then
    canonical="${base}"
  fi
  mv -f "${aff}" "${raw_tse_dir}/${canonical}"
  echo "[INFO] Installed ${raw_tse_dir}/${canonical}"
done

rm -f "${process_dir}"/*.BRIK "${process_dir}"/*.HEAD

echo "[INFO] Done refreshing anat->TSE affines."
