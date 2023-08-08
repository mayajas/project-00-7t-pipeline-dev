SUBJECT_ID="sub-0$1"
PRF_MODEL=prfpy_vol_fit_hrf_True
#_LinearDeveining_4layers

# do layerification in anat space or func space?
SPACE=func
#anat
#func
NLAYERS=6

# define directories
LAYNII_DIR=/home/mayajas/Documents/programs/laynii
RIM_SCRIPT_DIR=/home/mayajas/Documents/project-00-7t-pipeline-dev/code/analysis-scripts/python
PRFPY_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/${PRF_MODEL}/$SUBJECT_ID
COREG_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/manualcorr/coreg/_subject_id_$SUBJECT_ID/struct2func
if [ "$SUBJECT_ID" = "sub-04" ]; then
    FS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUBJECT_ID/autorecon_pial
else
    FS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial
fi

if [ "$SPACE" = "func" ]; then
    LAYER_DIR=$PRFPY_DIR/layerification_func
else
    LAYER_DIR=$PRFPY_DIR/layerification
fi


# define anatomical and func images
FUNC_IMG=${PRFPY_DIR}/meanFunc.nii

# make layer output directory
if [ ! -d "$LAYER_DIR" ]; then
  mkdir $LAYER_DIR
  echo "Making directory ${LAYER_DIR}..."
else
  echo "Layerification outputs will be saved to ${LAYER_DIR}..."
fi

# convert ribbon mgz to nii
mri_convert $FS_DIR/$SUBJECT_ID/mri/ribbon.mgz $LAYER_DIR/ribbon.nii


if [ "$SPACE" = "anat" ]; then
    # make laynii required rim file in anat space
    python $RIM_SCRIPT_DIR/makeLayniiRim.py $LAYER_DIR/ribbon.nii

    # laynii layerification
    $LAYNII_DIR/LN2_LAYERS -rim $LAYER_DIR/ribbon_rim.nii -nr_layers ${NLAYERS} -equivol
    # $LAYNII_DIR/LN_GROW_LAYERS -rim $LAYER_DIR/ribbon_rim.nii -N ${NLAYERS} -thin

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
    #$LAYNII_DIR/LN_GROW_LAYERS -rim $LAYER_DIR/func_ribbon_rim.nii -N ${NLAYERS} -thin
    
    

fi