SUBJECT_ID=sub-01
FS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial_rerun
COREG_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_laminar_fMRI_func_pRF/_subject_id_$SUBJECT_ID/coreg
OUT_DIR=$FS_DIR/$SUBJECT_ID/mri

MOVING_IMG=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_laminar_fMRI_func_pRF/_subject_id_sub-01/meanFunc/merged_func_mcf.nii_mean_reg.nii

rGM_val=42
rWM_val=41
lGM_val=3
lWM_val=2

# right hemisphere
FIXED_IMG=$FS_DIR/$SUBJECT_ID/mri/rGM.nii

mri_binarize --i $FS_DIR/$SUBJECT_ID/mri/ribbon.mgz --o $FS_DIR/$SUBJECT_ID/mri/rGM.mgz --min $rGM_val --max $rGM_val --dilate 1

mri_convert $FS_DIR/$SUBJECT_ID/mri/rGM.mgz $FIXED_IMG

${ANTSPATH}antsApplyTransforms \
  -d 3 \
  -i $FIXED_IMG \
  -r $MOVING_IMG \
  -t [ $COREG_DIR/registered_0GenericAffine.mat, 1 ] \
  -t $COREG_DIR/registered_1InverseWarp.nii.gz \
  -o $OUT_DIR/rGM_funcSpace.nii \
  --interpolation NearestNeighbor

# left hemisphere
FIXED_IMG=$FS_DIR/$SUBJECT_ID/mri/lGM.nii

mri_binarize --i $FS_DIR/$SUBJECT_ID/mri/ribbon.mgz --o $FS_DIR/$SUBJECT_ID/mri/lGM.mgz --min $lGM_val --max $lGM_val --dilate 1

mri_convert $FS_DIR/$SUBJECT_ID/mri/lGM.mgz $FIXED_IMG

${ANTSPATH}antsApplyTransforms \
  -d 3 \
  -i $FIXED_IMG \
  -r $MOVING_IMG \
  -t [ $COREG_DIR/registered_0GenericAffine.mat, 1 ] \
  -t $COREG_DIR/registered_1InverseWarp.nii.gz \
  -o $OUT_DIR/lGM_funcSpace.nii \
  --interpolation NearestNeighbor