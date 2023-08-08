SUBJECT_ID="sub-0$1"
PROJECT_ID=project-00-7t-pipeline-dev

echo "Preparing fieldmap for subject $SUBJECT_ID"

PROJECT_DIR=/home/mayajas/scratch/$PROJECT_ID

FUNC_IMG=$PROJECT_DIR/raw/${SUBJECT_ID}/func/task-bar_run-01.nii
PHASEDIFF_IMG=$PROJECT_DIR/raw/${SUBJECT_ID}/fmap/phasediff.nii
MAGNITUDE_IMG=$PROJECT_DIR/raw/${SUBJECT_ID}/fmap/magnitude1.nii
n_dummy=4 

OUTDIR=$PROJECT_DIR/output/fmap/_subject_id_${SUBJECT_ID}

# make output directory
if [ ! -d "$OUTDIR" ]; then
  mkdir $OUTDIR
  echo "Making directory ${OUTDIR}..."
else
  echo "Outputs will be saved to ${OUTDIR}..."
fi

# copy images to outdir
cd $OUTDIR
if [ ! -f "$OUTDIR/phasediff.nii" ]; then
    cp ${PHASEDIFF_IMG} $OUTDIR/phasediff.nii
    echo "Copied phasediff.nii to output directory"
else
    echo "Output directory already contains phasediff.nii"
fi
if [ ! -f "$OUTDIR/magnitude1.nii" ]; then
    cp ${MAGNITUDE_IMG} $OUTDIR/magnitude1.nii
    echo "Copied magnitude1.nii to output directory"
else
    echo "Output directory already contains magnitude1.nii"
fi
if [ ! -f "$OUTDIR/ref_func.nii" ]; then
    echo ${FUNC_IMG} 
    fslroi ${FUNC_IMG} $OUTDIR/ref_func.nii $n_dummy 1
    echo "Copied ref_func.nii to output directory"
else
    echo "Output directory already contains ref_func.nii"
fi

######################################################################
## prepare fieldmap
cd $OUTDIR

# first extract brain from magnitude image
mri_synthstrip -i magnitude1.nii -o magnitude1_brain.nii -b 0

# prepare fmap
deltaTE=2.46
fsl_prepare_fieldmap SIEMENS phasediff.nii magnitude1_brain.nii fmap_rads $deltaTE

#####################################################################
# reslice magnitude and phasediff images to ref_func

# bias correct mean functional
N4BiasFieldCorrection -d 3 -v 1 -s 4 -b [ 180 ] -c [ 50x50x50x50, 0.0 ] \
-i $OUTDIR/ref_func.nii -o [ corrected_ref_func.nii, ref_func_BiasField.nii ]

# coreg fieldmap to mean func
PREFIX=MAGNITUDE2FUNC
#ITKSNAP_TRANSFORM=$OUTDIR/coreg_itksnap_${PREFIX}.txt
REF_IMG=$OUTDIR/corrected_ref_func.nii
MOVING_IMG=$OUTDIR/magnitude1.nii

antsRegistration --dimensionality 3 \
--float 1 --initialize-transforms-per-stage 0 \
--interpolation Linear --output [ reg${PREFIX}_, reg${PREFIX}_Warped.nii.gz, reg${PREFIX}_InverseWarped.nii.gz ] \
--transform Rigid[ 0.1 ] --metric MI[ $REF_IMG, $MOVING_IMG, 1, 32, Regular, 0.25 ] \
--convergence [ 1000x500x250x100, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
--transform Affine[ 0.1 ] --metric MI[ $REF_IMG, $MOVING_IMG, 1, 32, Regular, 0.25 ] \
--convergence [ 1000x500x250x100, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
--winsorize-image-intensities [ 0.005, 0.995 ]  --write-composite-transform 1 --collapse-output-transforms 0
#    --initial-moving-transform [ $ITKSNAP_TRANSFORM, 0 ] \

# apply the transform to magnitude image (use ref_func resolution)
MOVING_IMG=$OUTDIR/magnitude1.nii
OUTPUT_IMG=$OUTDIR/funcReg_magnitude1.nii
antsApplyTransforms --default-value 0 -d 3 --input $MOVING_IMG \
--interpolation Linear \
--output $OUTPUT_IMG \
--reference-image $REF_IMG \
--transform $OUTDIR/regMAGNITUDE2FUNC_Composite.h5

# apply the transform to fieldmap (use ref_func resolution)
MOVING_IMG=$OUTDIR/fmap_rads.nii
OUTPUT_IMG=$OUTDIR/funcReg_fmap_rads.nii
antsApplyTransforms --default-value 0 -d 3 --input $MOVING_IMG \
--interpolation Linear \
--output $OUTPUT_IMG \
--reference-image $REF_IMG \
--transform $OUTDIR/regMAGNITUDE2FUNC_Composite.h5

######################################################################
## apply fieldmap correction to mean functional

fugue -i $OUTDIR/ref_func.nii --loadfmap=funcReg_fmap_rads.nii --dwell=0.00499241 --unwarp=$OUTDIR/unwarped_ref_func.nii -v --unwarpdir=z-