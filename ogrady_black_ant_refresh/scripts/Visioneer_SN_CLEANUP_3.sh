#!/bin/bash -l
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -J tse_cleanup_masks
#SBATCH -o sbatch-tse-cleanup-%j.out
#SBATCH -e sbatch-tse-cleanup-%j.err

module load afni
export FREESURFER_HOME=${FREESURFER_HOME:-/path/to/freesurfer}
export SUBJECTS_DIR=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/segmentation
export FUNCTIONALS_DIR=$SUBJECTS_DIR
set +u
source $FREESURFER_HOME/SetUpFreeSurfer.sh
set -u
#
# Generate the MOTIP cleanup prerequisites for correct_TSE_forSN_phase4.m
# Combines tse_phase1 + tse_phase3 steps for a fixed set of HUMAN_EBR IDs:
#   1) copy/convert raw TSE into the shared TSE folder
#   2) extract brainstem/DC/cerebellum masks from FreeSurfer
#   3) zero-pad and align the masks + MPRAGE into TSE space (writes *_inTSE)
#

set -euo pipefail

if [[ -n ${VISIONEER_CLEANUP_IDS:-} ]]; then
  read -r -a cleanup_ids <<<"${VISIONEER_CLEANUP_IDS}"
  echo "[INFO] Visioneer cleanup restricted to: ${cleanup_ids[*]}"
else
  cleanup_ids=(530)
fi

motip_root=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data
snheart_root=${SNHEART_ROOT:-/path/to/SNHEART}
sn_mask_dir=${snheart_root}/MRI_data/TSE/probabilistic_masks
sub_file=${motip_root}/locations/TSE_subjects.txt
tse_loc_file=${motip_root}/locations/TSE_loc.txt
tmp_mask_dir=${motip_root}/TSE/tmp_sn_cleanup_masks
mkdir -p "$tmp_mask_dir"


# Helper to resolve the long folder + TSE session name for a subject
lookup_entry() {
  local sub=$1
  local line
  line=$(grep -n "${sub}$" "$sub_file" | head -n1 | cut -d: -f1 || true)
  [[ -n $line ]] || return 1
  main_folder=$(sed -n "${line}p" "$sub_file")
  tse_folder=$(sed -n "${line}p" "$tse_loc_file")
  return 0
}

for sub in "${cleanup_ids[@]}"; do
  echo "=== HUMAN_EBR-${sub} ==="
  if ! lookup_entry "$sub"; then
    echo "  !! Could not locate entry in ${sub_file}; skipping." >&2
    continue
  fi

  subj_root="${motip_root}/${main_folder}/${tse_folder}"
  raw_dst="${motip_root}/TSE/${sub}_raw_tse.nii"

  if [[ -f $raw_dst ]]; then
    echo "  -> Raw TSE already copied."
  else
    if [[ $sub -lt 500 ]]; then
      src_glob="${subj_root}/frfseopt*"
      converted="${subj_root}/frfseopt.nii"
    else
      src_glob="${subj_root}/TSE*"
      converted="${subj_root}/TSE.nii"
    fi
    echo "  -> Converting ${src_glob}"
    3dAFNItoNIFTI -overwrite -prefix "$converted" $src_glob
    echo "  -> Copying to ${raw_dst}"
    cp "$converted" "$raw_dst"
  fi

  # Ensure anatomical masks exist (brainstem/DC/cerebellum)
  seg_dir="${motip_root}/segmentation/${sub}/mri"
  aparc_mgz="${seg_dir}/aparc.a2009s+aseg.mgz"
  aparc_nii="${seg_dir}/aparc.a2009s+aseg.nii"
  [[ -f $aparc_nii ]] || mri_convert "$aparc_mgz" "$aparc_nii"

  if [[ ! -f ${motip_root}/TSE/brainstem_4v_mask_${sub}.nii ]]; then
    echo "  -> Extracting brainstem/4V masks"
    mri_extract_label "$aparc_nii" 15 16 "${motip_root}/TSE/brainstem_4v_mask_${sub}.nii"
    mri_extract_label "$aparc_nii" 16 "${motip_root}/TSE/brainstem_mask_${sub}.nii"
  fi
  if [[ ! -f ${motip_root}/TSE/SN/DC_mask_${sub}.nii ]]; then
    echo "  -> Extracting DC mask"
    mri_extract_label "$aparc_nii" 28 60 "${motip_root}/TSE/SN/DC_mask_${sub}.nii"
  fi
  if [[ ! -f ${motip_root}/TSE/cerebellum/cerebellum_mask_${sub}.nii ]]; then
    echo "  -> Extracting cerebellum control mask"
    mkdir -p "${motip_root}/TSE/cerebellum"
    mri_extract_label "$aparc_nii" 8 47 "${motip_root}/TSE/cerebellum/cerebellum_mask_${sub}.nii"
  fi

  # Align masks into zeropadded TSE space (phase 3 equivalent)
  mat="${motip_root}/TSE/${sub}_anat_toTSE_mat.aff12.1D"
  if [[ ! -f $mat ]]; then
    echo "  !! Missing ${mat}; run align_mprage_to_tse_phase2 for ${sub}." >&2
    continue
  fi

  raw_tse="${motip_root}/TSE/${sub}_raw_tse.nii"
  zeropad="${motip_root}/TSE/zeropad_tse_${sub}.nii"
  if [[ -f $zeropad ]]; then
    echo "  -> Zeropad volume already exists."
  else
    3dZeropad -I 40 -S 40 -A 40 -P 40 -L 40 -R -40 -Z 40 -prefix "$zeropad" "$raw_tse"
  fi

  mprage="${motip_root}/basic_anat/${main_folder}/anat.nii"
  anat_in_tse="${motip_root}/TSE/ANATinTSE_${sub}.nii"
  if [[ ! -f $anat_in_tse ]]; then
    3dAllineate -final wsinc5 -master "$zeropad" -1Dmatrix_apply "$mat" \
      -input "$mprage" -prefix "$anat_in_tse"
  else
    echo "  -> ANATinTSE already exists."
  fi

  brainstem4v_in_tse="${motip_root}/TSE/brainstem_4v_mask_inTSE_${sub}.nii"
  if [[ ! -f $brainstem4v_in_tse ]]; then
    3dAllineate -final wsinc5 -master "$zeropad" -1Dmatrix_apply "$mat" \
      -input "${motip_root}/TSE/brainstem_4v_mask_${sub}.nii" \
      -prefix "$brainstem4v_in_tse"
  else
    echo "  -> brainstem_4v_mask_inTSE already exists."
  fi

  brainstem_in_tse="${motip_root}/TSE/brainstem_mask_inTSE_${sub}.nii"
  if [[ ! -f $brainstem_in_tse ]]; then
    3dAllineate -final wsinc5 -master "$zeropad" -1Dmatrix_apply "$mat" \
      -input "${motip_root}/TSE/brainstem_mask_${sub}.nii" \
      -prefix "$brainstem_in_tse"
  else
    echo "  -> brainstem_mask_inTSE already exists."
  fi

  mkdir -p "${motip_root}/TSE/SN"
  dc_in_tse="${motip_root}/TSE/SN/DC_mask_inTSE_${sub}.nii"
  sn_mask=$(find "${sn_mask_dir}" -maxdepth 1 -name "sn_groupmask_in_${sub}.nii*" | head -n1 || true)
  if [[ -n $sn_mask ]]; then
    echo "  -> Rebuilding DC_mask_inTSE from SN group mask with AFNI dilation"
    tmp_mask_base="${tmp_mask_dir}/sn_mask_${sub}.nii"
    3dmask_tool -input "$sn_mask" -dilate_input 1 -prefix "${tmp_mask_base}" >/dev/null 2>&1 || true
    if [[ -f ${tmp_mask_base} ]]; then
      mv "${tmp_mask_base}" "$dc_in_tse"
    else
      echo "  !! Failed to create dilated SN mask for ${sub}; falling back to legacy DC mask warp." >&2
    fi
  fi
  if [[ ! -f $dc_in_tse ]]; then
    echo "  -> (Fallback) warping FreeSurfer DC mask"
    3dAllineate -final wsinc5 -master "$zeropad" -1Dmatrix_apply "$mat" \
      -input "${motip_root}/TSE/SN/DC_mask_${sub}.nii" \
      -prefix "$dc_in_tse"
  fi

  mkdir -p "${motip_root}/TSE/cerebellum"
  cereb_in_tse="${motip_root}/TSE/cerebellum/cerebellum_mask_inTSE_${sub}.nii"
  if [[ ! -f $cereb_in_tse ]]; then
    3dAllineate -final wsinc5 -master "$zeropad" -1Dmatrix_apply "$mat" \
      -input "${motip_root}/TSE/cerebellum/cerebellum_mask_${sub}.nii" \
      -prefix "$cereb_in_tse"
  else
    echo "  -> cerebellum_mask_inTSE already exists."
  fi

done

echo "Cleanup TSE masks complete."
