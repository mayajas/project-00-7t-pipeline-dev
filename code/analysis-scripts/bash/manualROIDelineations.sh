SUBJECT_ID="sub-0$1"
HEMI="$2h"
DELVER="$3"

PRF_MODEL=prfpy_fit_hrf_True_start_from_avg_True
#prfpy_vol_fit_hrf_True

if [ "$SUBJECT_ID" = "sub-04" ]; then
    FS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUBJECT_ID/autorecon_pial/$SUBJECT_ID
else
    FS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial/$SUBJECT_ID
fi
PRFPY_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/${PRF_MODEL}

if [ "$DELVER" = "rois" ]; then
    freeview -f $FS_DIR/surf/$HEMI.inflated:overlay=$PRFPY_DIR/$SUBJECT_ID/$HEMI.ecc.mgh:overlay_color='colorwheel':overlay=$PRFPY_DIR/$SUBJECT_ID/$HEMI.pol.mgh:overlay_custom=$PRFPY_DIR/${HEMI}_pol:label=$FS_DIR/label/$HEMI.V1_exvivo.thresh.label:label_outline=True
elif [ "$DELVER" = "ecc2" ]; then
    freeview -f $FS_DIR/surf/$HEMI.inflated:overlay=$PRFPY_DIR/$SUBJECT_ID/$HEMI.ecc.mgh:overlay_custom=$PRFPY_DIR/ecc2band:overlay=$PRFPY_DIR/$SUBJECT_ID/$HEMI.pol.mgh:overlay_custom=$PRFPY_DIR/${HEMI}_pol:label=$PRFPY_DIR/$SUBJECT_ID/$HEMI.V1.label:label_outline=True:label=$PRFPY_DIR/$SUBJECT_ID/$HEMI.V2d.label:label_outline=True:label=$PRFPY_DIR/$SUBJECT_ID/$HEMI.V3d.label:label_outline=True:label=$PRFPY_DIR/$SUBJECT_ID/$HEMI.V2v.label:label_outline=True:label=$PRFPY_DIR/$SUBJECT_ID/$HEMI.V3v.label:label_outline=True:label=$PRFPY_DIR/$SUBJECT_ID/$HEMI.V4.label:label_outline=True
elif [ "$DELVER" = "dv" ]; then
    freeview -f $FS_DIR/surf/$HEMI.inflated:overlay=$PRFPY_DIR/$SUBJECT_ID/$HEMI.y.mgh:overlay=$PRFPY_DIR/$SUBJECT_ID/$HEMI.ecc.mgh:overlay_color='colorwheel':overlay=$PRFPY_DIR/$SUBJECT_ID/$HEMI.pol.mgh:overlay_custom=$PRFPY_DIR/${HEMI}_pol:label=$PRFPY_DIR/$SUBJECT_ID/$HEMI.V2d.label:label_outline=True:label=$PRFPY_DIR/$SUBJECT_ID/$HEMI.V3d.label:label_outline=True:label=$PRFPY_DIR/$SUBJECT_ID/$HEMI.V2v.label:label_outline=True:label=$PRFPY_DIR/$SUBJECT_ID/$HEMI.V3v.label:label_outline=True:label=$PRFPY_DIR/$SUBJECT_ID/$HEMI.V4.label:label_outline=True
fi