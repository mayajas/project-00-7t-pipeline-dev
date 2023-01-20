#!/bin/bash

#SBATCH --job-name=antsReg3
#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=500MB
#SBATCH --cpus-per-task=20
#SBATCH --time=0:30:00
#SBATCH --qos=prio

# activate py36 environment
# include for this reason: https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
eval "$(conda shell.bash hook)"
conda activate py36

# load modules
module purge
module add ANTs/2.3.1-foss-2018b-Python-3.6.6

# get environment variable of num threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$SLURM_CPUS_PER_TASK

# unset Java max heap space option
unset _JAVA_OPTIONS

# define dirs
SUB_ID=sub-04
REF_IMG=/scratch/mayaaj90/project-00-7t-pipeline-dev/output/anat/_subject_id_$SUB_ID/UNI_corrected.nii
#REF_IMG=/scratch/mayaaj90/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_$SUB_ID/maskedUNI/UNI_corrected_masked.nii
MEANFUNC_IMG=/scratch/mayaaj90/project-00-7t-pipeline-dev/derivatives/wf_laminar_fMRI_func_pRF/_subject_id_$SUB_ID/meanFunc/merged_func_mcf.nii_mean_reg.nii

MASK_DIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/manualcorr/coreg/_subject_id_$SUB_ID
MASK=midoccMask.nii
#MASK=funcSpaceOccipitalMask.nii
OUTDIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/manualcorr/coreg/_subject_id_$SUB_ID
ITKSNAP_TRANSFORM=/scratch/mayaaj90/project-00-7t-pipeline-dev/manualcorr/coreg/_subject_id_$SUB_ID/coreg_itksnap_func2struct.txt

# make output dir
if [[ ! -d $OUTDIR ]]
then
	mkdir -p $OUTDIR
fi
cd $OUTDIR

# copy necessary files
cp $REF_IMG $OUTDIR
cp $MEANFUNC_IMG $OUTDIR


# run registration
if [[ ! -f "$OUTDIR/registered_1Warp.nii.gz" ]]
then
	echo "Ants registration..."
	antsRegistration --collapse-output-transforms 1 --dimensionality 3 --float 1 \
	--initial-moving-transform [ $ITKSNAP_TRANSFORM, 0 ] \
	--initialize-transforms-per-stage 0 --interpolation Linear --output [ registered_, registered_Warped.nii.gz, registered_InverseWarped.nii.gz ] \
	--transform Rigid[ 0.1 ] \
	--metric MI[ $REF_IMG, $MEANFUNC_IMG, 1, 32, Regular, 0.25 ] \
	--convergence [ 1000x500x250x100, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
	--masks [ $MASK_DIR/$MASK, NULL ] \
	--transform Affine[ 0.1 ] \
	--metric MI[ $REF_IMG, $MEANFUNC_IMG, 1, 32, Regular, 0.25 ] \
	--convergence [ 1000x500x250x100, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
	--masks [ $MASK_DIR/$MASK, NULL ] \
	--transform SyN[ 0.1, 3.0, 0.0 ] \
	--metric CC[ $REF_IMG, $MEANFUNC_IMG, 1, 4, None ] \
	--convergence [ 100x70x50x20, 1e-06, 10 ] --smoothing-sigmas 3.0x2.0x1.0x0.0 --shrink-factors 8x4x2x1 --use-histogram-matching 0 \
	--masks [ $MASK_DIR/$MASK, NULL ] -v --winsorize-image-intensities [ 0.005, 0.995 ]  --write-composite-transform 0
fi

# apply registration
if [[ ! -f "$OUTDIR/reg_meanFunc.nii" ]]
then
	antsApplyTransforms --default-value 0 --float 0 --input $MEANFUNC_IMG --interpolation BSpline[ 5 ] --output $OUTDIR/reg_meanFunc_nonlinear.nii \
	--reference-image $REF_IMG --transform $OUTDIR/registered_0GenericAffine.mat --transform $OUTDIR/registered_1Warp.nii.gz

	antsApplyTransforms --default-value 0 --float 0 --input $MEANFUNC_IMG --interpolation BSpline[ 5 ] --output $OUTDIR/reg_meanFunc_affine.nii \
        --reference-image $REF_IMG --transform $OUTDIR/registered_0GenericAffine.mat
fi
