SUBJECT_ID=sub-03
PROJECT_ID=project-00-7t-pipeline-dev

UNI=/home/mayajas/scratch/$PROJECT_ID/output/anat/_subject_id_$SUBJECT_ID/UNI_corrected.nii
manualcorr_dir=/home/mayajas/scratch/project-00-7t-pipeline-dev/manualcorr/freesurfer_externalbrainmask/_subject_id_$SUBJECT_ID
# FS_orig=/home/mayajas/scratch/$PROJECT_ID/output/freesurfer_externalbrainmask/_subject_id_$SUBJECT_ID
# FS_corr=/home/mayajas/scratch/$PROJECT_ID/output/freesurfer_regenerated/_subject_id_$SUBJECT_ID
FS_orig=/home/mayajas/scratch/$PROJECT_ID/output/freesurfer_regenerated/_subject_id_$SUBJECT_ID
FS_corr=/home/mayajas/scratch/$PROJECT_ID/output/freesurfer_reregenerated/_subject_id_$SUBJECT_ID

wm_mgz=$manualcorr_dir/wm.mgz


[ -f "$wm_mgz" ] || cp $FS_orig/wm.mgz $wm_mgz


freeview -v $UNI \
$manualcorr_dir/brain.finalsurfs.mgz \
$wm_mgz:colormap=heat:opacity=0.7 \
-f $FS_corr/lh.white:edgecolor=blue \
$FS_corr/lh.pial:edgecolor=red \
$FS_corr/rh.white:edgecolor=blue \
$FS_corr/rh.pial:edgecolor=red
