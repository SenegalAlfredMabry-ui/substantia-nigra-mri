#!/bin/bash -l
# supervillan_wallet.sh
# -----------------------------------------------------------------------------
# PURPOSE
#   Rebuild the missing SSwarper outputs (subject MPRAGE -> MNI152_2009 warp
#   field + affine) for the MOTIP IDs that keep blocking
#   build_black_ant_masks.sh. Running this script from the SNHEART repo means we
#   don't have to cd into 33_MOTIP2018 manually.
# OUTPUTS
#   - For each requested subject, writes `${motip_root}/MRI_data/WARPS/<ID>/`
#     contents exactly as sswarper2 normally produces (anatQQ.<ID>_WARP.nii,
#     anatQQ.<ID>_aff12.1D, QC snapshots, etc.). build_black_ant_masks.sh and
#     S.H.I.E.L.D_SN_Warp_v1.sh already look in that location, so no path edits
#     are required downstream.
# DEPENDENCIES
#   - Legacy AFNI virtualenv (provides sswarper2) and the MOTIP data tree.
#   - `MPRAGE_subjects.txt` + `MPRAGE_loc.txt` to map HUMAN_EBR labels to scan
#     folders. The wrapper mirrors warper.sh but constrains the subject roster
#     to the 22 IDs we know are missing warps.
# RELIES ON
#   - Later scripts (build_black_ant_masks.sh, S.H.I.E.L.D_SN_Warp_v1.sh,
#     irredeemable_black_ant_V1.m) require these warp fields to exist so they
#     can move atlases between spaces. Run this wrapper whenever those scripts
#     complain about `anat.un.aff.qw_WARP.nii` being absent.
# -----------------------------------------------------------------------------
#SBATCH -J supervillan_wallet
#SBATCH -p cloud
#SBATCH -c 4
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o sbatch-wrap-warper-%j.out
#SBATCH -e sbatch-wrap-warper-%j.err
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=FAIL

set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
cd "$repo_root"

source ${MOTIP_ROOT:-/path/to/MOTIP2018}/afnienv/bin/activate
module load afni

echo '#####################################################'
echo 'supervillan_wallet: SSWarper refresh for selected MOTIP subjects'
echo '#####################################################'

motip_root=${MOTIP_ROOT:-/path/to/MOTIP2018}
base_folder=${motip_root}
data_folder=${motip_root}/MRI_data
mni_template=${motip_root}/MNI152_2009_template_SSW.nii.gz

snheart_repo=${SNHEART_ROOT:-/path/to/SNHEART}
snheart_root=${snheart_repo}/MRI_data
default_roster_file=${snheart_root}/locations/TSE_SN_ASR_SUBLIST.txt
if [[ ! -f ${default_roster_file} ]]; then
  default_roster_file=${snheart_root}/locations/TSE_subjects.txt
fi

warp_root=${SUPERVILLAIN_WARP_ROOT:-${repo_root}/MRI_data/WARPS_supervillain}
mkdir -p "${warp_root}"

anat_subject_file=${SUPERVILLAIN_MPRAGE_SUBJECTS:-${data_folder}/locations/MPRAGE_subjects.txt}
anat_scan_loc=${SUPERVILLAIN_MPRAGE_LOC:-${data_folder}/locations/MPRAGE_loc.txt}

if [[ -f ${default_roster_file} ]]; then
  mapfile -t default_subjects < <(grep -Ev '^[[:space:]]*(#|$)' "${default_roster_file}")
else
  # Public GitHub copy: use an env/file override or replace these example IDs.
  default_subjects=(100 101 102)
fi
if [[ -n ${SUPERVILLAIN_SUBJECTS:-} ]]; then
  read -r -a subjects <<<"${SUPERVILLAIN_SUBJECTS}"
elif [[ -n ${SUPERVILLAIN_SUBJECT_FILE:-} && -f ${SUPERVILLAIN_SUBJECT_FILE} ]]; then
  mapfile -t subjects <"${SUPERVILLAIN_SUBJECT_FILE}"
else
  subjects=(${default_subjects[@]})
fi

echo "[INFO] Refreshing SSWarper outputs for: ${subjects[*]}"

overwrite=${SUPERVILLAIN_OVERWRITE:-0}
manifest_path=${warp_root}/warp_manifest.tsv
if [[ ! -f ${manifest_path} ]]; then
  printf "subject\twarp_dir\tmprage\tstatus\ttimestamp\n" >"${manifest_path}"
fi

missing_subjects=0

for entry in "${subjects[@]}"; do
  entry_trimmed=$(echo "${entry}" | xargs)
  if [[ ${entry_trimmed} =~ ^[0-9]+$ ]]; then
    sub=$(printf '%03d' "${entry_trimmed}")
    subject_label=$(printf 'HUMAN_EBR-%03d' "${entry_trimmed}")
  else
    subject_label=${entry_trimmed}
    sub=${entry_trimmed##*-}
  fi
  sub=${sub##0}
  printf -v sub "%03d" "${sub}"
  echo "\n=== ${subject_label} ==="
  match=$(grep -n "^${subject_label}$" "${anat_subject_file}" || true)
  if [[ -z ${match} ]]; then
    echo "  !! ${subject_label} not found in ${anat_subject_file}; skipping"
    ((missing_subjects+=1))
    continue
  fi
  line=${match%%:*}
  scan_num=$(sed -n "${line}p" "${anat_scan_loc}")
  if [[ -z ${scan_num} ]]; then
    echo "  !! Missing scan number entry for line ${line}; skipping"
    ((missing_subjects+=1))
    continue
  fi

  subject_dir=${data_folder}/${subject_label}/${scan_num}
  if [[ ! -d ${subject_dir} ]]; then
    echo "  !! Missing directory ${subject_dir}; skipping"
    ((missing_subjects+=1))
    continue
  fi

  if [[ -f ${subject_dir}/anat.nii && ! -f ${subject_dir}/MPRAGE.nii ]]; then
    cp "${subject_dir}/anat.nii" "${subject_dir}/MPRAGE.nii"
  fi
  if [[ -f ${subject_dir}/mprage.nii && ! -f ${subject_dir}/MPRAGE.nii ]]; then
    cp "${subject_dir}/mprage.nii" "${subject_dir}/MPRAGE.nii"
  fi

  mprage=${subject_dir}/MPRAGE.nii
  if [[ ! -f ${mprage} ]]; then
    echo "  !! No MPRAGE.nii found in ${subject_dir}; skipping"
    ((missing_subjects+=1))
    continue
  fi

  warp_dir=${warp_root}/${sub}
  warp_field=${warp_dir}/anatQQ.${sub}_WARP.nii
  if [[ -f ${warp_field} && ${overwrite} -eq 0 ]]; then
    echo "  -> Warp already exists at ${warp_field}; skipping"
    printf "%s\t%s\t%s\tskipped_existing\t%s\n" \
      "${sub}" "${warp_dir}" "${mprage}" "$(date -u +%FT%TZ)" >>"${manifest_path}"
    continue
  fi

  mkdir -p "${warp_dir}"
  echo "  -> Running sswarper2 (output ${warp_dir})"
  if ! sswarper2 -input "${mprage}" \
      -base "${mni_template}" \
      -subid "${sub}" \
      -odir "${warp_dir}"; then
    echo "  !! sswarper2 failed for ${sub}" >&2
    printf "%s\t%s\t%s\tfailure\t%s\n" \
      "${sub}" "${warp_dir}" "${mprage}" "$(date -u +%FT%TZ)" >>"${manifest_path}"
    ((missing_subjects+=1))
    continue
  fi

  for core in "${warp_dir}/anatQQ.${sub}" "${warp_dir}/anatQQ.${sub}_WARP"; do
    if [[ -f ${core}.nii.gz && ! -f ${core}.nii ]]; then
      gunzip -c "${core}.nii.gz" >"${core}.nii"
    fi
  done

  if [[ ! -f ${warp_dir}/anatQQ.${sub}_WARP.nii || ! -f ${warp_dir}/anatQQ.${sub}.aff12.1D ]]; then
    echo "  !! Missing warp outputs for ${sub} even after sswarper" >&2
    printf "%s\t%s\t%s\tmissing_outputs\t%s\n" \
      "${sub}" "${warp_dir}" "${mprage}" "$(date -u +%FT%TZ)" >>"${manifest_path}"
    ((missing_subjects+=1))
    continue
  fi

  printf "%s\t%s\t%s\tsuccess\t%s\n" \
    "${sub}" "${warp_dir}" "${mprage}" "$(date -u +%FT%TZ)" >>"${manifest_path}"

done

echo "[INFO] Wrapper complete (missing subjects: ${missing_subjects})"
