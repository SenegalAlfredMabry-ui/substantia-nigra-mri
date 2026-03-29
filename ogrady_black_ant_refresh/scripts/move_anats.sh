#!/bin/bash -l

#
#
#
# TO RUN THIS SCRIPT, you just type "./move_anats.sh" or "sh move_anats.sh" from within the correct folder, 
# This script attempts to make convenient anat files from **any subject in the MPRAGE_subjects.txt file***
# So make sure you have run make_scan_loc_files.txt before running this, or it may not process your newest subjects
# Written by Lissa Riley
#
#

data_dir=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data
sub_file=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/locations/MPRAGE_subjects.txt
mprage_loc_file=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/locations/MPRAGE_loc.txt

# Optional override: space-delimited numeric IDs, e.g. "647 648 650"
override_ids=${MOVE_ANATS_SUBJECTS:-}

cd ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/segmentation

numsub=`wc -l < ${sub_file}`

for i in $(seq 1 $numsub); do

sub_folder=`sed "${i}q;d" ${sub_file}`
echo $sub_folder

sub=${sub_folder: -3}
if [[ -n "$override_ids" ]]; then
	match=0
	for wanted in $override_ids; do
		wanted_fmt=$(printf "%03d" "$wanted")
		if [[ "$sub" == "$wanted_fmt" ]]; then
			match=1
			break
		fi
	done
	if [[ $match -eq 0 ]]; then
		continue
	fi
fi

anat_dir=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder

if ! [[ -f $anat_dir/anat.nii ]]
then
	mkdir -p $anat_dir
	mprage_folder=`sed "${i}q;d" $mprage_loc_file`
	echo $mprage_folder
	mprage_file1=${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/$sub_folder/$mprage_folder/mprage.nii
	mprage_file2=`ls ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/$sub_folder/$mprage_folder/mprage_scan*.nii`
	mprage_file3=`ls ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/$sub_folder/$mprage_folder/mprage_scan*.BRIK`
	mprage_file4=`ls ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/$sub_folder/$mprage_folder/MPRAGE.nii`
	
	if [[ -f $mprage_file1 ]]
	then
		echo "First format file found..."
		cp ${mprage_file1} ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder
		mv ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder/mprage.nii ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder/anat.nii
	elif [[ -f $mprage_file2 ]]
	then
		echo "Second format file found..."
		cp ${mprage_file2} ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder
		mv ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder/mprage_scan*.nii ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder/anat.nii
	elif [[ -f $mprage_file4 ]]
	then
		echo "Second format file found..."
		cp ${mprage_file4} ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder
		mv ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder/MPRAGE.nii ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder/anat.nii
	elif [[ -f $mprage_file3 ]]
	then
		echo "Only BRIK/HEAD format currently exists..."
		3dAFNItoNIFTI -prefix ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/$sub_folder/$mprage_folder/mprage.nii ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/$sub_folder/$mprage_folder/mprage_scan*.BRIK
		cp ${mprage_file1} ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder
		mv ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder/mprage.nii ${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/basic_anat/$sub_folder/anat.nii		
	else
		echo "Could not find an appropriate MPRAGE file to copy to basic_anat for ${sub_folder}"
	fi
else	
	echo "anat.nii already exists for ${sub_folder}"
fi



done
