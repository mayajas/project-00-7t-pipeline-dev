SUBJECT_ID=sub-02
PROJECT_ID=project-00-7t-pipeline-dev

PROJECT_DIR=/home/mayajas/scratch/$PROJECT_ID

if [ "$SUBJECT_ID" = "sub-04" ]; then
    FS_DIR=$PROJECT_DIR/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUBJECT_ID/autorecon_pial
else
    FS_DIR=$PROJECT_DIR/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial
fi

#UNI=/home/mayajas/scratch/$PROJECT_ID/output/anat/_subject_id_$SUBJECT_ID/UNI_corrected.nii
manualcorr_dir=/home/mayajas/scratch/$PROJECT_ID/manualcorr/freesurfer_externalbrainmask/_subject_id_$SUBJECT_ID

freeview -v $FS_DIR/$SUBJECT_ID/mri/T1.mgz \
$FS_DIR/$SUBJECT_ID/mri/brain.finalsurfs.mgz \
$PROJECT_DIR/output/prfpy_prepped_inputs/$SUBJECT_ID/reg_meanFunc.nii \
-f $FS_DIR/$SUBJECT_ID/surf/lh.white:edgecolor=yellow \
$FS_DIR/$SUBJECT_ID/surf/lh.pial:edgecolor=red \
$FS_DIR/$SUBJECT_ID/surf/rh.white:edgecolor=yellow \
$FS_DIR/$SUBJECT_ID/surf/rh.pial:edgecolor=red
