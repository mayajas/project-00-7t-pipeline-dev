SUBJECT_ID="sub-0$1"
HEMI="$2h"

if [ "$SUBJECT_ID" = "sub-04" ]; then
    FS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUBJECT_ID/autorecon_pial_rerun/$SUBJECT_ID
else
    FS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial_rerun/$SUBJECT_ID
fi
PRFPY_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/prfpy

freeview -f $FS_DIR/surf/$HEMI.inflated:overlay=$PRFPY_DIR/$SUBJECT_ID/$HEMI.ecc.mgh:overlay_color='colorwheel':overlay=$PRFPY_DIR/$SUBJECT_ID/$HEMI.pol.mgh:overlay_custom=$PRFPY_DIR/${HEMI}_pol:label=$FS_DIR/label/$HEMI.V1_exvivo.thresh.label:label_outline=True