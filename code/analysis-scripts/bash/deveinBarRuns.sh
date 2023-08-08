SUBJECT_ID="sub-0$1"
PROJECT_ID=project-00-7t-pipeline-dev
# do layerification in anat space or func space?
SPACE=anat
#anat
#func

# how many layers?
NLAYERS=6

# which deveining method to use?
DEVEINING_METHOD=linear
#linear
#cbv
#leakage

# define directories
PROJECT_DIR=/home/mayajas/scratch/${PROJECT_ID}
LAYNII_DIR=/home/mayajas/Documents/programs/laynii
RIM_SCRIPT_DIR=/home/mayajas/Documents/${PROJECT_ID}/code/analysis-scripts/python
OUT_DIR=${PROJECT_DIR}/output/deveining/$SUBJECT_ID

if [ ! -d "${PROJECT_DIR}/output/deveining" ]; then
  mkdir ${PROJECT_DIR}/output/deveining
  echo "Making directory ${PROJECT_DIR}/output/deveining..."
else
  echo "Outputs will be saved to ${PROJECT_DIR}/output/deveining..."
fi

COREG_DIR=${PROJECT_DIR}/manualcorr/coreg/_subject_id_$SUBJECT_ID/struct2func
if [ "$SUBJECT_ID" = "sub-04" ]; then
    FS_DIR=${PROJECT_DIR}/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUBJECT_ID/autorecon_pial_rerun
else
    FS_DIR=${PROJECT_DIR}/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial_rerun
fi

if [ "$SPACE" = "func" ]; then
    LAYER_DIR=$OUT_DIR/layerification_func
else
    LAYER_DIR=$OUT_DIR/layerification
fi
DEVEIN_DIR=$OUT_DIR/deveined_${DEVEINING_METHOD}


# define func images
FUNC_IMG=${PROJECT_DIR}/output/func/meanFunc/_subject_id_${SUBJECT_ID}/merged_func_mcf.nii_mean_reg.nii
BAR1_IMG=${PROJECT_DIR}/output/func/sliceTimeCorr/_subject_id_${SUBJECT_ID}/_sess_id_task-bar_run-01_sess_nr_0_sess_nvol_124/atask-bar_run-01_roi_warp4D.nii
BAR2_IMG=${PROJECT_DIR}/output/func/sliceTimeCorr/_subject_id_${SUBJECT_ID}/_sess_id_task-bar_run-02_sess_nr_1_sess_nvol_124/atask-bar_run-02_roi_warp4D.nii

# make layer output directory
if [ ! -d "$OUT_DIR" ]; then
  mkdir $OUT_DIR
  mkdir $LAYER_DIR
  mkdir $DEVEIN_DIR
  echo "Making directory ${OUT_DIR}..."
else
  echo "Deveining outputs will be saved to ${OUT_DIR}..."
fi

# convert ribbon mgz to nii
mri_convert $FS_DIR/$SUBJECT_ID/mri/ribbon.mgz $LAYER_DIR/ribbon.nii


if [ "$SPACE" = "anat" ]; then
    # make laynii required rim file in anat space
    python $RIM_SCRIPT_DIR/makeLayniiRim.py $LAYER_DIR/ribbon.nii

    # laynii layerification
    $LAYNII_DIR/LN2_LAYERS -rim $LAYER_DIR/ribbon_rim.nii -nr_layers ${NLAYERS} -equivol

    # coreg to func space
    # equidistant image:
    INPUT=$LAYER_DIR/ribbon_rim_layers_equidist.nii
    OUTPUT=$LAYER_DIR/func_ribbon_rim_layers_equidist.nii   
    antsApplyTransforms --default-value 0 \
    -d 3 \
    --input $INPUT \
    --interpolation NearestNeighbor \
    --output $OUTPUT \
    --reference-image $FUNC_IMG \
    --transform $COREG_DIR/regT12UNI_Composite.h5 \
    --transform $COREG_DIR/regUNI2FUNC_Composite.h5

    # equivol image:
    INPUT=$LAYER_DIR/ribbon_rim_layers_equivol.nii
    OUTPUT=$LAYER_DIR/func_ribbon_rim_layers_equivol.nii
    antsApplyTransforms --default-value 0 \
    -d 3 \
    --input $INPUT \
    --interpolation NearestNeighbor \
    --output $OUTPUT \
    --reference-image $FUNC_IMG \
    --transform $COREG_DIR/regT12UNI_Composite.h5 \
    --transform $COREG_DIR/regUNI2FUNC_Composite.h5

else
    # coreg ribbon file to func space
    INPUT=$LAYER_DIR/ribbon.nii
    OUTPUT=$LAYER_DIR/func_ribbon.nii
        
    antsApplyTransforms --default-value 0 \
    -d 3 \
    --input $INPUT \
    --interpolation NearestNeighbor \
    --output $OUTPUT \
    --reference-image $FUNC_IMG \
    --transform $COREG_DIR/regT12UNI_Composite.h5 \
    --transform $COREG_DIR/regUNI2FUNC_Composite.h5

    # make laynii required rim file in anat space
    python $RIM_SCRIPT_DIR/makeLayniiRim.py $LAYER_DIR/func_ribbon.nii

    # laynii layerification
    $LAYNII_DIR/LN2_LAYERS -rim $LAYER_DIR/func_ribbon_rim.nii -nr_layers ${NLAYERS} -equivol
    

fi

# copy bar runs to deveining directory
if [ ! -f "${DEVEIN_DIR}/bar1.nii" ]; then
    cp $BAR1_IMG $DEVEIN_DIR/bar1.nii
    cp $BAR2_IMG $DEVEIN_DIR/bar2.nii
fi

# run deveining
if [ "$DEVEINING_METHOD" = "linear" ]; then
    $LAYNII_DIR/LN2_DEVEIN -input $DEVEIN_DIR/bar1.nii -layer_file $LAYER_DIR/func_ribbon_rim_layers_equidist.nii -linear
    $LAYNII_DIR/LN2_DEVEIN -input $DEVEIN_DIR/bar2.nii -layer_file $LAYER_DIR/func_ribbon_rim_layers_equidist.nii -linear
# elif [ "$DEVEINING_METHOD" = "cbv" ]; then
fi