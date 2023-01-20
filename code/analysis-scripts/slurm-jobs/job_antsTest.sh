#!/bin/bash

#SBATCH --job-name=ants
#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=500MB
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
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

# run the thing
SUB_ID=sub-01
PRFPY_DIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/output/prfpy/$SUB_ID
REF_IMG=/scratch/mayaaj90/project-00-7t-pipeline-dev/output/anat/_subject_id_$SUB_ID/UNI_corrected.nii
TRANSFORM_DIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/derivatives/wf_laminar_fMRI_func_pRF/_subject_id_sub-01/coreg

cd $PRFPY_DIR

PARAM=l_ecc
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/$PARAM.nii --interpolation BSpline[ 5 ] --output $PRFPY_DIR/reg_$PARAM.nii --reference-image $REF_IMG --transform $TRANSFORM_DIR/registered_0GenericAffine.mat --transform $TRANSFORM_DIR/registered_1Warp.nii.gz

PARAM=l_pol
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/$PARAM.nii --interpolation BSpline[ 5 ] --output $PRFPY_DIR/reg_$PARAM.nii --reference-image $REF_IMG --transform $TRANSFORM_DIR/regist$

PARAM=l_rsq
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/$PARAM.nii --interpolation BSpline[ 5 ] --output $PRFPY_DIR/reg_$PARAM.nii --reference-image $REF_IMG --transform $TRANSFORM_DIR/regist$

PARAM=l_sigma
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/$PARAM.nii --interpolation BSpline[ 5 ] --output $PRFPY_DIR/reg_$PARAM.nii --reference-image $REF_IMG --transform $TRANSFORM_DIR/regist$

PARAM=l_x
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/$PARAM.nii --interpolation BSpline[ 5 ] --output $PRFPY_DIR/reg_$PARAM.nii --reference-image $REF_IMG --transform $TRANSFORM_DIR/regist$

PARAM=l_y
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/$PARAM.nii --interpolation BSpline[ 5 ] --output $PRFPY_DIR/reg_$PARAM.nii --reference-image $REF_IMG --transform $TRANSFORM_DIR/regist$


