#!/bin/bash

#SBATCH --job-name=prfpy-params
#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=500M
#SBATCH --cpus-per-task=1
#SBATCH --time=0:15:00
#SBATCH --qos=prio

# activate py36 environment
# include for this reason: https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
eval "$(conda shell.bash hook)"
conda activate py36

# load modules
module purge
module load FreeSurfer/dev-centos7_x86_64 ; source $FREESURFER_HOME/SetUpFreeSurfer.sh
module add ANTs/2.3.1-foss-2018b-Python-3.6.6


# get environment variable of num threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$SLURM_CPUS_PER_TASK

# unset Java max heap space option
unset _JAVA_OPTIONS

# array of combinations for sbatch array
declare -a combinations
index=0
for sub in 0 1 2 3
do
    for hem in 0 1
    do
      	combinations[$index]="$sub $hem"
        index=$((index + 1))
    done
done

parameters=(${combinations[${SLURM_ARRAY_TASK_ID}]})

sub=${parameters[0]}
hem=${parameters[1]}

# vars and directories
SUBJECTS=(sub-01 sub-02 sub-03 sub-04)
HEMISPHERES=(lh rh)

SUB_ID=${SUBJECTS[sub]}
HEM_ID=${HEMISPHERES[hem]}

PROJ_FRAC=0
INTERP_ANTS=NearestNeighbor
#BSpline[5]
INTERP_VOL2SURF=trilin
#trilin
#nearest

PRFPY_DIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/output/prfpy/$SUB_ID
REF_IMG=/scratch/mayaaj90/project-00-7t-pipeline-dev/output/anat/_subject_id_$SUB_ID/UNI_corrected.nii
TRANSFORM_DIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/manualcorr/coreg/_subject_id_$SUB_ID

if [ "$SUB_ID" = "sub-04" ]; then
    SUBJECTS_DIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUB_ID/autorecon_pial_rerun
else
    SUBJECTS_DIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_$SUB_ID/autorecon_pial_rerun
fi

MEAN_FUNC=/scratch/mayaaj90/project-00-7t-pipeline-dev/output/func/meanFunc/_subject_id_$SUB_ID/merged_func_mcf.nii_mean_reg.nii

# change directory to prfpy outputs
cd $PRFPY_DIR

# apply transform and surface project mean functional
antsApplyTransforms --default-value 0 --float 0 --input $MEAN_FUNC --interpolation $INTERP_ANTS --output $PRFPY_DIR/reg_meanFunc.nii \
--reference-image $REF_IMG --transform $TRANSFORM_DIR/registered_0GenericAffine.mat

mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.meanFunc.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_meanFunc.nii


# apply transform and surface project bar data
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/bar.nii --interpolation $INTERP_ANTS --output $PRFPY_DIR/reg_bar.nii \
--reference-image $REF_IMG --transform $TRANSFORM_DIR/registered_0GenericAffine.mat --input-image-type 3 --verbose

mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.bar.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_bar.nii

fslreorient2std $PRFPY_DIR/reg_bar.nii $PRFPY_DIR/reg_bar_fsl.nii 

# apply transform and surface project pRF parameters
PARAM=ecc
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/${HEM_ID}_$PARAM.nii --interpolation $INTERP_ANTS --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
--reference-image $REF_IMG --transform $TRANSFORM_DIR/registered_0GenericAffine.mat 

mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii



PARAM=pol
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/${HEM_ID}_$PARAM.nii --interpolation $INTERP_ANTS --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
--reference-image $REF_IMG --transform $TRANSFORM_DIR/registered_0GenericAffine.mat 

mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii



PARAM=sigma
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/${HEM_ID}_$PARAM.nii --interpolation $INTERP_ANTS --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
--reference-image $REF_IMG --transform $TRANSFORM_DIR/registered_0GenericAffine.mat 

mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii


PARAM=rsq
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/${HEM_ID}_$PARAM.nii --interpolation $INTERP_ANTS --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
--reference-image $REF_IMG --transform $TRANSFORM_DIR/registered_0GenericAffine.mat 

mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii



PARAM=x
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/${HEM_ID}_$PARAM.nii --interpolation $INTERP_ANTS --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
--reference-image $REF_IMG --transform $TRANSFORM_DIR/registered_0GenericAffine.mat 

mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii



PARAM=y
antsApplyTransforms --default-value 0 --float 0 --input $PRFPY_DIR/${HEM_ID}_$PARAM.nii --interpolation $INTERP_ANTS --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
--reference-image $REF_IMG --transform $TRANSFORM_DIR/registered_0GenericAffine.mat 

mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii
