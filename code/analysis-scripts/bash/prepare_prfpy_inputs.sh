SUBJECT_ID="sub-0$1"
PROJECT_ID=project-00-7t-pipeline-dev
PRF_MODEL=prfpy_vol

echo "Preparing prfpy inputs for subject $SUBJECT_ID"

###########################################################################################################################
## directories, filenames, important variables
PROJECT_DIR=/home/mayajas/scratch/$PROJECT_ID
#PROJECT_DIR=/scratch/mayaaj90/$PROJECT_ID
if [ "$SUBJECT_ID" = "sub-04" ]; then
    FS_DIR=$PROJECT_DIR/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUBJECT_ID/autorecon_pial
else
    FS_DIR=$PROJECT_DIR/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial
fi
SUBDIR=$PROJECT_DIR/manualcorr/coreg/_subject_id_$SUBJECT_ID
COREG_DIR=$SUBDIR/func2struct
OUT_DIR=$PROJECT_DIR/output/${PRF_MODEL}/$SUBJECT_ID

MEANFUNC_IMG=$PROJECT_DIR/output/func/meanFunc/_subject_id_$SUBJECT_ID/merged_func_mcf.nii_mean_reg.nii
BAR1_IMG=$PROJECT_DIR/output/func/sliceTimeCorr/_subject_id_$SUBJECT_ID/_sess_id_task-bar_run-01_sess_nr_0_sess_nvol_124/atask-bar_run-01_roi_warp4D_unwarped.nii
BAR2_IMG=$PROJECT_DIR/output/func/sliceTimeCorr/_subject_id_$SUBJECT_ID/_sess_id_task-bar_run-02_sess_nr_1_sess_nvol_124/atask-bar_run-02_roi_warp4D_unwarped.nii
T1_IMG=$PROJECT_DIR/output/func/anat/_subject_id_$SUBJECT_ID/T1_out.nii

rWM_val=41
lWM_val=2
rGM_val=42
lGM_val=3

###########################################################################################################################
## make output directories, copy needed files
# make prfpy general output directory
if [ ! -d "$PROJECT_DIR/output/${PRF_MODEL}" ]; then
  mkdir $PROJECT_DIR/output/${PRF_MODEL}
fi

# make prfpy subject output directory
if [ ! -d "$OUT_DIR" ]; then
  mkdir $OUT_DIR
  echo "Making directory ${OUT_DIR}..."
else
  echo "Outputs will be saved to ${OUT_DIR}..."
fi

if [ ! -f "$OUT_DIR/ribbon.mgz" ]; then
  cp $FS_DIR/$SUBJECT_ID/mri/ribbon.mgz $OUT_DIR/ribbon.mgz
fi
if [ ! -f "$OUT_DIR/meanFunc.nii" ]; then
  cp $MEANFUNC_IMG $OUT_DIR/meanFunc.nii
fi
if [ ! -f "$OUT_DIR/T1_out.nii" ]; then
  cp $T1_IMG $OUT_DIR/T1_out.nii
fi
if [ ! -f "$OUT_DIR/bar1.nii" ]; then
  cp $BAR1_IMG $OUT_DIR/bar1.nii
fi
if [ ! -f "$OUT_DIR/bar2.nii" ]; then
  cp $BAR2_IMG $OUT_DIR/bar2.nii
fi
###########################################################################################################################
# coregister mean functional & bar run volumes to T1 resolution
if [ ! -f "$OUT_DIR/reg_meanFunc.nii" ]; then
  antsApplyTransforms --default-value 0 \
      --input $OUT_DIR/meanFunc.nii \
      --interpolation BSpline[5] \
      --output $OUT_DIR/reg_meanFunc.nii \
      --reference-image $OUT_DIR/meanFunc.nii \
      --transform $COREG_DIR/regFUNC2UNI_Composite.h5 \
      --transform $COREG_DIR/regUNI2T1_Composite.h5
fi

if [ ! -f "$OUT_DIR/reg_bar1.nii" ]; then
  antsApplyTransforms --default-value 0 \
      -e 3 \
      --input $OUT_DIR/bar1.nii \
      --interpolation BSpline[5] \
      --output $OUT_DIR/reg_bar1.nii \
      --reference-image $OUT_DIR/meanFunc.nii \
      --transform $COREG_DIR/regFUNC2UNI_Composite.h5 \
      --transform $COREG_DIR/regUNI2T1_Composite.h5
fi

if [ ! -f "$OUT_DIR/reg_bar2.nii" ]; then
  antsApplyTransforms --default-value 0 \
      -e 3 \
      --input $OUT_DIR/bar2.nii \
      --interpolation BSpline[5] \
      --output $OUT_DIR/reg_bar2.nii \
      --reference-image $OUT_DIR/meanFunc.nii \
      --transform $COREG_DIR/regFUNC2UNI_Composite.h5 \
      --transform $COREG_DIR/regUNI2T1_Composite.h5   
fi


if [ "$PRF_MODEL" = "prfpy_vol" ]; then
  ###########################################################################################################################
  ## make WM mask
  # right hemisphere
  if [ ! -f "$OUT_DIR/rh_WM.nii" ]; then
    MOVING_IMG=$OUT_DIR/rh_WM.nii

    mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/rh_WM.mgz --match $rWM_val 
    mri_convert $OUT_DIR/rh_WM.mgz $MOVING_IMG
  fi

  # left hemisphere
  if [ ! -f "$OUT_DIR/lh_WM.nii" ]; then
    MOVING_IMG=$OUT_DIR/lh_WM.nii

    mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/lh_WM.mgz --match $lWM_val 
    mri_convert $OUT_DIR/lh_WM.mgz $MOVING_IMG
  fi

  ###########################################################################################################################
  ## make GM mask
  ## right hemisphere
  if [ ! -f "$OUT_DIR/rh_GM_inflated_ds.nii" ]; then
    MOVING_IMG=$OUT_DIR/rh_GM_inflated.nii

    # binarize and dilate into WM
    mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/rh_GM.mgz --match $rGM_val 
    mri_convert $OUT_DIR/rh_GM.mgz $OUT_DIR/rh_GM.nii

    mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/rh_GM_dilated.mgz --match $rGM_val --dilate 2
    mri_convert $OUT_DIR/rh_GM_dilated.mgz $OUT_DIR/rh_GM_dilated.nii

    fslmaths $OUT_DIR/rh_GM_dilated.nii -sub $OUT_DIR/rh_GM.nii -mul $OUT_DIR/rh_WM.nii -add $OUT_DIR/rh_GM.nii $MOVING_IMG

    # reslice to functional dimensions
    antsApplyTransforms --default-value 0 \
        --input $OUT_DIR/rh_GM_inflated.nii \
        --interpolation NearestNeighbor \
        --output $OUT_DIR/rh_GM_inflated_ds.nii \
        --reference-image $OUT_DIR/meanFunc.nii \
        --transform identity
  fi

  ## left hemisphere
  if [ ! -f "$OUT_DIR/lh_GM_inflated_ds.nii" ]; then
    MOVING_IMG=$OUT_DIR/lh_GM_inflated.nii

    # binarize and dilate into WM
    mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/lh_GM.mgz --match $lGM_val 
    mri_convert $OUT_DIR/lh_GM.mgz $OUT_DIR/lh_GM.nii

    mri_binarize --i $OUT_DIR/ribbon.mgz --o $OUT_DIR/lh_GM_dilated.mgz --match $lGM_val --dilate 2
    mri_convert $OUT_DIR/lh_GM_dilated.mgz $OUT_DIR/lh_GM_dilated.nii

    fslmaths $OUT_DIR/lh_GM_dilated.nii -sub $OUT_DIR/lh_GM.nii -mul $OUT_DIR/lh_WM.nii -add $OUT_DIR/lh_GM.nii $MOVING_IMG

    antsApplyTransforms --default-value 0 \
        --input $OUT_DIR/lh_GM_inflated.nii \
        --interpolation NearestNeighbor \
        --output $OUT_DIR/lh_GM_inflated_ds.nii \
        --reference-image $OUT_DIR/meanFunc.nii \
        --transform identity
  fi

fi

# ###########################################################################################################################
# ## make occipital label from volumetric (manual) mask
# SUBJECTS_DIR=$FS_DIR

# # lh
# mri_vol2surf --mov $OCC_MASK --regheader $SUBJECT_ID --hemi lh --out $OUT_DIR/lh.occMask.mgh --interp nearest --projfrac 0.5 --cortex
# mri_vol2label --i $OUT_DIR/lh.occMask.mgh --id 1 --l $OUT_DIR/lh.occMask.label --surf $SUBJECT_ID lh

# # rh
# mri_vol2surf --mov $OCC_MASK --regheader $SUBJECT_ID --hemi rh --out $OUT_DIR/rh.occMask.mgh --interp nearest --projfrac 0.5 --cortex
# mri_vol2label --i $OUT_DIR/rh.occMask.mgh --id 1 --l $OUT_DIR/rh.occMask.label --surf $SUBJECT_ID rh

