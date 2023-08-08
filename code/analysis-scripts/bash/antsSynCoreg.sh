SUBJECT_ID="sub-0$1"
PROJECT_ID=project-00-7t-pipeline-dev

WHICHCOREG=struct2func

echo "SyN coregistration for subject $SUBJECT_ID"

PROJECT_DIR=/home/mayajas/scratch/$PROJECT_ID
if [ "$SUBJECT_ID" = "sub-04" ]; then
    FS_DIR=$PROJECT_DIR/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUBJECT_ID/autorecon_pial_rerun
else
    FS_DIR=$PROJECT_DIR/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial_rerun
fi

T1_IMG=$PROJECT_DIR/output/func/anat/_subject_id_${SUBJECT_ID}/T1_out.nii
UNI_IMG=$PROJECT_DIR/output/anat/_subject_id_${SUBJECT_ID}/UNI_corrected.nii
MEANFUNC_IMG=$PROJECT_DIR/output/func/meanFunc/_subject_id_${SUBJECT_ID}/merged_func_mcf.nii_mean_reg.nii

OUTDIR=$PROJECT_DIR/manualcorr/coreg/_subject_id_${SUBJECT_ID}/${WHICHCOREG}_syn

# make prfpy output directory
if [ ! -d "$OUTDIR" ]; then
  mkdir $OUTDIR
  echo "Making directory ${OUTDIR}..."
else
  echo "Outputs will be saved to ${OUTDIR}..."
fi

# copy images to outdir
cd $OUTDIR
if [ ! -f "$OUTDIR/T1.nii" ]; then
    cp ${T1_IMG} $OUTDIR/T1.nii
    echo "Copied T1.nii to output directory"
else
    echo "Output directory already contains T1.nii"
fi
if [ ! -f "$OUTDIR/UNI.nii" ]; then
    cp ${UNI_IMG} $OUTDIR/UNI.nii
    echo "Copied UNI.nii to output directory"
else
    echo "Output directory already contains UNI.nii"
fi
if [ ! -f "$OUTDIR/meanFunc.nii" ]; then
    cp ${MEANFUNC_IMG} $OUTDIR/meanFunc.nii
    echo "Copied meanFunc.nii to output directory"
else
    echo "Output directory already contains meanFunc.nii"
fi

# bias correct mean functional
N4BiasFieldCorrection -d 3 -v 1 -s 4 -b [ 180 ] -c [ 50x50x50x50, 0.0 ] \
  -i $OUTDIR/meanFunc.nii -o [ corrected_meanFunc.nii, meanFunc_BiasField.nii ]

if [ "$WHICHCOREG" = "struct2func" ]; then
    # coreg T1 to UNI corrected
    PREFIX=T12UNI
    ITKSNAP_TRANSFORM=$OUTDIR/coreg_itksnap_${PREFIX}.txt
    REF_IMG=$OUTDIR/UNI.nii
    MOVING_IMG=$OUTDIR/T1.nii
    antsRegistration --dimensionality 3 --float 1 \
            --initial-moving-transform [ $ITKSNAP_TRANSFORM, 0 ] \
            --initialize-transforms-per-stage 0 --interpolation Linear --output [ reg${PREFIX}_, reg${PREFIX}_Warped.nii.gz, reg${PREFIX}_InverseWarped.nii.gz ] \
            --transform Rigid[ 0.1 ] \
            --metric MI[ $REF_IMG, $MOVING_IMG, 1, 32, Regular, 0.25 ] \
            --convergence [ 1000x500x250x100, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
            --transform Affine[ 0.1 ] \
            --metric MI[ $REF_IMG, $MOVING_IMG, 1, 32, Regular, 0.25 ] \
            --convergence [ 1000x500x250x100, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
            --winsorize-image-intensities [ 0.005, 0.995 ]  --write-composite-transform 1 --collapse-output-transforms 0


    # coreg UNI to mean func
    PREFIX=UNI2FUNC
    ITKSNAP_TRANSFORM=$OUTDIR/coreg_itksnap_${PREFIX}.txt
    REF_IMG=$OUTDIR/corrected_meanFunc.nii
    MOVING_IMG=$OUTDIR/UNI.nii
    MANUAL_MASK=$PROJECT_DIR/manualcorr/coreg/_subject_id_${SUBJECT_ID}/fixedImageMask_bigger.nii
    antsRegistration --dimensionality 3 --float 1 \
            --initial-moving-transform [ $ITKSNAP_TRANSFORM, 0 ] \
            --initialize-transforms-per-stage 0 --interpolation Linear --output [ reg${PREFIX}_, reg${PREFIX}_Warped.nii.gz, reg${PREFIX}_InverseWarped.nii.gz ] \
            --transform Rigid[ 0.1 ] \
            --metric MI[ $REF_IMG, $MOVING_IMG, 1, 32, Regular, 0.25 ] \
            --convergence [ 1000x500x250x100, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
            --masks [ NULL, ${MANUAL_MASK} ] \
            --transform Affine[ 0.1 ] \
            --metric MI[ $REF_IMG, $MOVING_IMG, 1, 32, Regular, 0.25 ] \
            --convergence [ 1000x500x250x100, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
            --masks [ NULL, ${MANUAL_MASK} ] \
            --winsorize-image-intensities [ 0.005, 0.995 ]  --write-composite-transform 1 --collapse-output-transforms 0 \
            --transform SyN[ 0.1, 3.0, 0.0 ] \
            --metric CC[ $REF_IMG, $MOVING_IMG, 1, 4, None ] \
            --convergence [ 100x70x50x20, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
            --masks [ NULL, ${MANUAL_MASK} ]

    # apply both transforms to T1 (keep T1 resolution)
    MOVING_IMG=$OUTDIR/T1.nii
    REF_IMG=$OUTDIR/T1.nii
    OUTPUT_IMG=$OUTDIR/funcReg_T1.nii
    antsApplyTransforms --default-value 0 \
        -d 3 \
        --input $MOVING_IMG \
        --interpolation Linear \
        --output $OUTPUT_IMG \
        --reference-image $REF_IMG \
        --transform $OUTDIR/regT12UNI_Composite.h5 \
        --transform $OUTDIR/regUNI2FUNC_Composite.h5
else
    # coreg meanFunc to UNI corrected
    PREFIX=FUNC2UNI
    ITKSNAP_TRANSFORM=$OUTDIR/coreg_itksnap_${PREFIX}.txt
    REF_IMG=$OUTDIR/UNI.nii
    MOVING_IMG=$OUTDIR/corrected_meanFunc.nii
    MANUAL_MASK=$PROJECT_DIR/manualcorr/coreg/_subject_id_${SUBJECT_ID}/fixedImageMask_bigger.nii
    antsRegistration --dimensionality 3 --float 1 \
        --initial-moving-transform [ $ITKSNAP_TRANSFORM, 0 ] \
        --initialize-transforms-per-stage 0 --interpolation Linear --output [ reg${PREFIX}_, reg${PREFIX}_Warped.nii.gz, reg${PREFIX}_InverseWarped.nii.gz ] \
        --transform Rigid[ 0.1 ] \
        --metric MI[ $REF_IMG, $MOVING_IMG, 1, 32, Regular, 0.25 ] \
        --convergence [ 1000x500x250x100, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
        --masks [ ${MANUAL_MASK}, NULL ] \
        --transform Affine[ 0.1 ] \
        --metric MI[ $REF_IMG, $MOVING_IMG, 1, 32, Regular, 0.25 ] \
        --convergence [ 1000x500x250x100, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
        --masks [ ${MANUAL_MASK}, NULL ] \
        --winsorize-image-intensities [ 0.005, 0.995 ]  --write-composite-transform 1 --collapse-output-transforms 0 \
        --transform SyN[ 0.1, 3.0, 0.0 ] \
        --metric CC[ $REF_IMG, $MOVING_IMG, 1, 4, None ] \
        --convergence [ 100x70x50x20, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
        --masks [ ${MANUAL_MASK}, NULL ]

    # coreg UNI to T1
    PREFIX=UNI2T1
    ITKSNAP_TRANSFORM=$OUTDIR/coreg_itksnap_${PREFIX}.txt
    REF_IMG=$OUTDIR/T1.nii
    MOVING_IMG=$OUTDIR/UNI.nii
    MANUAL_MASK=$PROJECT_DIR/manualcorr/coreg/_subject_id_${SUBJECT_ID}/fixedImageMask_bigger.nii
    antsRegistration --dimensionality 3 --float 1 \
            --initial-moving-transform [ $ITKSNAP_TRANSFORM, 0 ] \
            --initialize-transforms-per-stage 0 --interpolation Linear --output [ reg${PREFIX}_, reg${PREFIX}_Warped.nii.gz, reg${PREFIX}_InverseWarped.nii.gz ] \
            --transform Rigid[ 0.1 ] \
            --metric MI[ $REF_IMG, $MOVING_IMG, 1, 32, Regular, 0.25 ] \
            --convergence [ 1000x500x250x100, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
            --transform Affine[ 0.1 ] \
            --metric MI[ $REF_IMG, $MOVING_IMG, 1, 32, Regular, 0.25 ] \
            --convergence [ 1000x500x250x100, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
            --winsorize-image-intensities [ 0.005, 0.995 ]  --write-composite-transform 1 --collapse-output-transforms 0 

    # apply both transforms to meanFunc (keep meanFunc resolution)
    MOVING_IMG=$OUTDIR/meanFunc.nii
    REF_IMG=$OUTDIR/meanFunc.nii
    OUTPUT_IMG=$OUTDIR/reg_meanFunc.nii
    antsApplyTransforms --default-value 0 \
        -d 3 \
        --input $MOVING_IMG \
        --interpolation Linear \
        --output $OUTPUT_IMG \
        --reference-image $REF_IMG \
        --transform $OUTDIR/regFUNC2UNI_Composite.h5 \
        --transform $OUTDIR/regUNI2T1_Composite.h5

fi