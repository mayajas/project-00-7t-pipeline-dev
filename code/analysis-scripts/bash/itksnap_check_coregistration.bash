SUBJECT=sub-04
cd ~/scratch/project-00-7t-pipeline-dev/
itksnap -g ./output/anat/_subject_id_$SUBJECT/UNI_corrected.nii -o ./output/func/coreg_calcarine/_subject_id_$SUBJECT/registered_Warped.nii.gz ./output/func/coreg_midocc/_subject_id_$SUBJECT/registered_Warped.nii.gz ./output/func/coreg_occipital/_subject_id_$SUBJECT/registered_Warped.nii.gz -s ./manualcorr/func/_subject_id_$SUBJECT/occipitalMask.nii
