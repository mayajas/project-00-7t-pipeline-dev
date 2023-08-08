SUBJECT_ID="sub-0$1"
PROJECT_ID=project-00-7t-pipeline-dev
PRF_MODEL=prfpy_Iso2DGaussianModel
#_LinearDeveining

echo "Making cortical masks for subject $SUBJECT_ID"

if [ "$SUBJECT_ID" = "sub-04" ]; then
    FS_DIR=/home/mayajas/scratch/$PROJECT_ID/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUBJECT_ID/autorecon_pial_rerun
else
    FS_DIR=/home/mayajas/scratch/$PROJECT_ID/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial_rerun
fi

COREG_DIR=/home/mayajas/scratch/$PROJECT_ID/manualcorr/coreg/_subject_id_$SUBJECT_ID/struct2func
OUT_DIR=/home/mayajas/scratch/$PROJECT_ID/output/${PRF_MODEL}/$SUBJECT_ID

# make prfpy output directory
if [ ! -d "/home/mayajas/scratch/$PROJECT_ID/output/${PRF_MODEL}" ]; then
  mkdir /home/mayajas/scratch/$PROJECT_ID/output/${PRF_MODEL}
fi

# make prfpy output directory
if [ ! -d "$OUT_DIR" ]; then
  mkdir $OUT_DIR
  echo "Making directory ${OUT_DIR}..."
else
  echo "Outputs will be saved to ${OUT_DIR}..."
fi

REF_IMG=/home/mayajas/scratch/$PROJECT_ID/derivatives/wf_laminar_fMRI_func_pRF/_subject_id_$SUBJECT_ID/meanFunc/merged_func_mcf.nii_mean_reg.nii
#MOVING_IMG=$OUT_DIR/meanFunc.nii
UNI_IMG=/home/mayajas/scratch/$PROJECT_ID/output/anat/_subject_id_$SUBJECT_ID/UNI_corrected.nii

rGM_val=42
rWM_val=41
lGM_val=3
lWM_val=2

cp $FS_DIR/$SUBJECT_ID/mri/ribbon.mgz $OUT_DIR/ribbon.mgz
cp $REF_IMG $OUT_DIR/meanFunc.nii
cp $UNI_IMG $OUT_DIR/UNI_corrected.nii

fslreorient2std $OUT_DIR/meanFunc.nii $OUT_DIR/meanFunc_fsl.nii

## WM
# right hemisphere
MOVING_IMG=$OUT_DIR/rh_WM.nii

mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/rh_WM.mgz --match $rWM_val 
mri_convert $OUT_DIR/rh_WM.mgz $MOVING_IMG

antsApplyTransforms --default-value 0 \
    -d 3 \
    --input $MOVING_IMG \
    --interpolation NearestNeighbor \
    --output $OUT_DIR/rh_WM_funcSpace.nii \
    --reference-image $REF_IMG \
    --transform $COREG_DIR/regT12UNI_Composite.h5 \
    --transform $COREG_DIR/regUNI2FUNC_Composite.h5

# left hemisphere
MOVING_IMG=$OUT_DIR/lh_WM.nii

mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/lh_WM.mgz --match $lWM_val 
mri_convert $OUT_DIR/lh_WM.mgz $MOVING_IMG

antsApplyTransforms --default-value 0 \
    -d 3 \
    --input $MOVING_IMG \
    --interpolation NearestNeighbor \
    --output $OUT_DIR/lh_WM_funcSpace.nii \
    --reference-image $REF_IMG \
    --transform $COREG_DIR/regT12UNI_Composite.h5 \
    --transform $COREG_DIR/regUNI2FUNC_Composite.h5

## GM - inflated into WM
# right hemisphere
MOVING_IMG=$OUT_DIR/rh_GM.nii

# binarize and dilate into WM
mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/rh_GM.mgz --match $rGM_val 
mri_convert $OUT_DIR/rh_GM.mgz $MOVING_IMG

mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/rh_GM_dilated.mgz --match $rGM_val --dilate 2
mri_convert $OUT_DIR/rh_GM_dilated.mgz $OUT_DIR/rh_GM_dilated.nii

fslmaths $OUT_DIR/rh_GM_dilated.nii -sub $MOVING_IMG -mul $OUT_DIR/rh_WM.nii -add $OUT_DIR/rh_GM.nii $OUT_DIR/rh_GM_submulladd.nii.gz
mri_convert $OUT_DIR/rh_GM_submulladd.nii.gz $MOVING_IMG

antsApplyTransforms --default-value 0 \
    -d 3 \
    --input $MOVING_IMG \
    --interpolation NearestNeighbor \
    --output $OUT_DIR/rh_GM_funcSpace.nii \
    --reference-image $REF_IMG \
    --transform $COREG_DIR/regT12UNI_Composite.h5 \
    --transform $COREG_DIR/regUNI2FUNC_Composite.h5

# left hemisphere
MOVING_IMG=$OUT_DIR/lh_GM.nii

# binarize and dilate into WM
mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/lh_GM.mgz --match $lGM_val 
mri_convert $OUT_DIR/lh_GM.mgz $MOVING_IMG

mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/lh_GM_dilated.mgz --match $lGM_val --dilate 2
mri_convert $OUT_DIR/lh_GM_dilated.mgz $OUT_DIR/lh_GM_dilated.nii

fslmaths $OUT_DIR/lh_GM_dilated.nii -sub $MOVING_IMG -mul $OUT_DIR/lh_WM.nii -add $OUT_DIR/lh_GM.nii $OUT_DIR/lh_GM_submulladd.nii.gz
mri_convert $OUT_DIR/lh_GM_submulladd.nii.gz $MOVING_IMG

antsApplyTransforms --default-value 0 \
    -d 3 \
    --input $MOVING_IMG \
    --interpolation NearestNeighbor \
    --output $OUT_DIR/lh_GM_funcSpace.nii \
    --reference-image $REF_IMG \
    --transform $COREG_DIR/regT12UNI_Composite.h5 \
    --transform $COREG_DIR/regUNI2FUNC_Composite.h5

## GM - un-inflated
# right hemisphere
MOVING_IMG=$OUT_DIR/rh_GM.nii

# binarize and dilate into WM
mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/rh_GM.mgz --match $rGM_val 
mri_convert $OUT_DIR/rh_GM.mgz $MOVING_IMG

antsApplyTransforms --default-value 0 \
    -d 3 \
    --input $MOVING_IMG \
    --interpolation NearestNeighbor \
    --output $OUT_DIR/rh_GM_funcSpace_uninflated.nii \
    --reference-image $REF_IMG \
    --transform $COREG_DIR/regT12UNI_Composite.h5 \
    --transform $COREG_DIR/regUNI2FUNC_Composite.h5

# left hemisphere
MOVING_IMG=$OUT_DIR/lh_GM.nii

# binarize and dilate into WM
mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/lh_GM.mgz --match $lGM_val 
mri_convert $OUT_DIR/lh_GM.mgz $MOVING_IMG

antsApplyTransforms --default-value 0 \
    -d 3 \
    --input $MOVING_IMG \
    --interpolation NearestNeighbor \
    --output $OUT_DIR/lh_GM_funcSpace_uninflated.nii \
    --reference-image $REF_IMG \
    --transform $COREG_DIR/regT12UNI_Composite.h5 \
    --transform $COREG_DIR/regUNI2FUNC_Composite.h5



## apply transform to UNI image
antsApplyTransforms --default-value 0 \
    -d 3 \
    --input $UNI_IMG \
    --interpolation Linear \
    --output $OUT_DIR/UNI_funcSpace.nii \
    --reference-image $REF_IMG \
    --transform $COREG_DIR/regUNI2FUNC_Composite.h5

fslreorient2std $OUT_DIR/UNI_funcSpace.nii $OUT_DIR/UNI_funcSpace_fsl.nii