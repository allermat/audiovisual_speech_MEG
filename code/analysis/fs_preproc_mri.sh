#!/bin/sh
#SBATCH --nodes 1-1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 48:0:00
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user mate.aller@mrc-cbu.cam.ac.uk

# set up freesurfer
export FSVER='6.0.0'
export FSDIR=${FSROOT}/${FSVER}
export FREESURFER_HOME=/imaging/local/software/freesurfer/${FSVER}/`arch`

echo $FREESURFER_HOME
echo $SUBJECTS_DIR

source $FREESURFER_HOME/FreeSurferEnv.sh

# root directory for processed MRI data
export SUBJECTS_DIR=/imaging/ma09/Projects/AVSpeechMEG/data/derivatives/anat
echo $SUBJECTS_DIR

# list of subject names here
for SUBJECT in sub-01 sub-02 sub-03 sub-04 sub-05 sub-06 sub-07 sub-08 sub-09 sub-10 sub-11 sub-12 sub-13 sub-14
 do
    RAW_DIR=/imaging/ma09/Projects/AVSpeechMEG/data/rawdata/${SUBJECT}/anat
    my_NIfTI=${RAW_DIR}/${SUBJECT}_T1w.nii
    recon-all -i $my_NIfTI -s $SUBJECT -all

done
