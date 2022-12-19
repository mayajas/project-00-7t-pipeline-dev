SUBJECT_ID=sub-02
PROJECT_ID=project-00-7t-pipeline-dev
UNI=/home/mayajas/scratch/$PROJECT_ID/output/anat/_subject_id_$SUBJECT_ID/UNI_corrected.nii
manualcorr_dir=/home/mayajas/scratch/project-00-7t-pipeline-dev/manualcorr/freesurfer_externalbrainmask/_subject_id_$SUBJECT_ID

freeview -v $UNI \
$manualcorr_dir/brain.finalsurfs.mgz \
-f $manualcorr_dir/lh.white:edgecolor=yellow \
$manualcorr_dir/lh.pial:edgecolor=red \
$manualcorr_dir/rh.white:edgecolor=yellow \
$manualcorr_dir/rh.pial:edgecolor=red
