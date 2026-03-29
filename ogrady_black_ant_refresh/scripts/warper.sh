#!/bin/bash -l
#This script calculates warps from individual space to standard space for everybody with oddball data

#Written by Lissa Riley
#the following commands tell the Slurm scheduler how many nodes and CPUs we want, what to call our job, etc.
#SBATCH -n 8
#SBATCH -N 8
#SBATCH -c 2
#SBATCH -J warper
#SBATCH -o sbatch-%j.out
#SBATCH -e sbatch-%j.err
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=ALL
#SBATCH --oversubscribe
#SBATCH --ntasks-per-core=1
#SBATCH -exclude=cluster-node-0


source ${MOTIP_ROOT:-/path/to/MOTIP2018}/afnienv/bin/activate


module load afni

echo '#####################################################'
echo 'SSWarper2'
echo '#####################################################'

# directories that will de-clutter code below
base_folder=${MOTIP_ROOT:-/path/to/MOTIP2018}
prepro_folder=${base_folder}/scripts/preprocess
data_folder=${base_folder}/MRI_data

# load scan location file for $task and $run
subject_folder=${data_folder}/locations/oddball_subjects.txt

# get scan location file for subject's MPRAGE
anat_scan_loc=${data_folder}/locations/MPRAGE_loc.txt
anat_subject_folder=${data_folder}/locations/MPRAGE_subjects.txt


subjects=`cat ${anat_subject_folder} `
echo $subjects


# loop through all subjects that you input above
for subject in $subjects
do
    #echo 'I AM IN THE MATRIX'
    sub_position_in_anat_scan_file=`grep -n ${subject} $anat_subject_folder | cut -d : -f 1`
    anat_scan_num=$(sed -n ${sub_position_in_anat_scan_file}p $anat_scan_loc)

    echo $subject
    sub=${subject: -3}


    if [[ -f ${data_folder}/HUMAN_EBR-${sub}/${anat_scan_num}/anat.nii ]];
    then
        cp ${data_folder}/HUMAN_EBR-${sub}/${anat_scan_num}/anat.nii ${data_folder}/HUMAN_EBR-${sub}/${anat_scan_num}/MPRAGE.nii
    fi

    if [[ -f ${data_folder}/HUMAN_EBR-${sub}/${anat_scan_num}/mprage.nii ]];
    then
        cp ${data_folder}/HUMAN_EBR-${sub}/${anat_scan_num}/mprage.nii ${data_folder}/HUMAN_EBR-${sub}/${anat_scan_num}/MPRAGE.nii
    fi
    

    #check if afni_proc script already exists for that subject, task, and run
    if [[ -f ${data_folder}/WARPS/${sub}/anatQQ.${sub}_WARP.nii ]];
    then
        echo "Warp already exists"
    else
        echo "Running SSWarper2..."
        sswarper2 -input ${data_folder}/HUMAN_EBR-${sub}/${anat_scan_num}/MPRAGE.nii \
        -base ${base_folder}/MNI152_2009_template_SSW.nii.gz \
        -subid ${sub} \
        -odir ${data_folder}/WARPS/${sub}
    fi
done
