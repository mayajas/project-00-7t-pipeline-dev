SUBJECT_ID=sub-04

if [ "$SUBJECT_ID" = "sub-04" ]; then
    FS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUBJECT_ID/autorecon_pial
else
    FS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial_rerun
fi



COREG_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/manualcorr/coreg/_subject_id_$SUBJECT_ID
OUT_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/prfpy/$SUBJECT_ID

MOVING_IMG=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_laminar_fMRI_func_pRF/_subject_id_$SUBJECT_ID/meanFunc/merged_func_mcf.nii_mean_reg.nii
ANAT_IMG=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/anat/_subject_id_$SUBJECT_ID/UNI_corrected.nii

rGM_val=42
rWM_val=41
lGM_val=3
lWM_val=2

cp $FS_DIR/$SUBJECT_ID/mri/ribbon.mgz $OUT_DIR/ribbon.mgz
cp $MOVING_IMG $OUT_DIR/meanFunc.nii
cp $ANAT_IMG $OUT_DIR/UNI_corrected.nii


fslreorient2std $OUT_DIR/meanFunc.nii $OUT_DIR/meanFunc_fsl.nii

# right hemisphere
FIXED_IMG=$OUT_DIR/rh_GM.nii

mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/rh_GM.mgz --min $rGM_val --max $rGM_val 
#--dilate 1

mri_convert $OUT_DIR/rh_GM.mgz $FIXED_IMG

# ${ANTSPATH}antsApplyTransforms \
#   -d 3 \
#   -i $FIXED_IMG \
#   -r $MOVING_IMG \
#   -t [ $COREG_DIR/registered_0GenericAffine.mat, 1 ] \
#   -t $COREG_DIR/registered_1InverseWarp.nii.gz \
#   -o $OUT_DIR/rh_GM_funcSpace.nii.nii \
#   --interpolation NearestNeighbor
${ANTSPATH}antsApplyTransforms \
  -d 3 \
  -i $FIXED_IMG \
  -r $MOVING_IMG \
  -t [ $COREG_DIR/registered_0GenericAffine.mat, 1 ] \
  -o $OUT_DIR/rh_GM_funcSpace.nii.nii \
  --interpolation NearestNeighbor

# left hemisphere
FIXED_IMG=$OUT_DIR/lh_GM.nii

mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/lh_GM.mgz --min $lGM_val --max $lGM_val 
#--dilate 1

mri_convert $OUT_DIR/lh_GM.mgz $FIXED_IMG

# ${ANTSPATH}antsApplyTransforms \
#   -d 3 \
#   -i $FIXED_IMG \
#   -r $MOVING_IMG \
#   -t [ $COREG_DIR/registered_0GenericAffine.mat, 1 ] \
#   -t $COREG_DIR/registered_1InverseWarp.nii.gz \
#   -o $OUT_DIR/lh_GM_funcSpace.nii.nii \
#   --interpolation NearestNeighbor
${ANTSPATH}antsApplyTransforms \
  -d 3 \
  -i $FIXED_IMG \
  -r $MOVING_IMG \
  -t [ $COREG_DIR/registered_0GenericAffine.mat, 1 ] \
  -o $OUT_DIR/lh_GM_funcSpace.nii.nii \
  --interpolation NearestNeighbor


# apply inverse transform to UNI image
${ANTSPATH}antsApplyTransforms \
  -d 3 \
  -i $ANAT_IMG \
  -r $MOVING_IMG \
  -t [ $COREG_DIR/registered_0GenericAffine.mat, 1 ] \
  -o $OUT_DIR/UNI_funcSpace.nii \
  --interpolation BSpline[5]

  fslreorient2std $OUT_DIR/UNI_funcSpace.nii $OUT_DIR/UNI_funcSpace_fsl.nii