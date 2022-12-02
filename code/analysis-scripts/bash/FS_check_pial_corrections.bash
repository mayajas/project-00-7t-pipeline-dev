SUBJECT_ID=sub-04
PROJECT_ID=project-00-7t-pipeline-dev
UNI=/home/mayajas/scratch/$PROJECT_ID/output/anat/_subject_id_$SUBJECT_ID/UNI_corrected.nii
FS_orig=/home/mayajas/scratch/$PROJECT_ID/output/freesurfer_externalbrainmask/_subject_id_$SUBJECT_ID
FS_corr=/home/mayajas/scratch/$PROJECT_ID/output/freesurfer_regenerated/_subject_id_$SUBJECT_ID

freeview -v $UNI \
-f $FS_orig/lh.white:edgecolor=red \
$FS_orig/lh.pial:edgecolor=red \
$FS_orig/rh.white:edgecolor=red \
$FS_orig/rh.pial:edgecolor=red \
$FS_corr/lh.white:edgecolor=yellow \
$FS_corr/lh.pial:edgecolor=yellow \
$FS_corr/rh.white:edgecolor=yellow \
$FS_corr/rh.pial:edgecolor=yellow
