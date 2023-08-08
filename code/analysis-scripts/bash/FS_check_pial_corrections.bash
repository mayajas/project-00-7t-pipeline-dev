SUBJECT_ID=sub-03
PROJECT_ID=project-00-7t-pipeline-dev
UNI=/home/mayajas/scratch/$PROJECT_ID/output/anat/_subject_id_$SUBJECT_ID/UNI_corrected.nii
T1=/home/mayajas/scratch/$PROJECT_ID/output/func/anat/_subject_id_$SUBJECT_ID/T1_out.nii
# FS_orig=/home/mayajas/scratch/$PROJECT_ID/output/freesurfer_externalbrainmask/_subject_id_$SUBJECT_ID
# FS_corr=/home/mayajas/scratch/$PROJECT_ID/output/freesurfer_regenerated/_subject_id_$SUBJECT_ID
FS_orig=/home/mayajas/scratch/$PROJECT_ID/output/freesurfer_externalbrainmask/_subject_id_$SUBJECT_ID
FS_corr=/home/mayajas/scratch/$PROJECT_ID/output/freesurfer_regenerated/_subject_id_$SUBJECT_ID

freeview -v $T1 \
-f $FS_orig/lh.white:edgecolor=red \
$FS_orig/lh.pial:edgecolor=red \
$FS_orig/rh.white:edgecolor=red \
$FS_orig/rh.pial:edgecolor=red \
$FS_corr/lh.white:edgecolor=blue \
$FS_corr/lh.pial:edgecolor=yellow \
$FS_corr/rh.white:edgecolor=blue \
$FS_corr/rh.pial:edgecolor=yellow
