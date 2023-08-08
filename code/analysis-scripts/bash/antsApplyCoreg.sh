SUB_ID="sub-0$1"
HEM_ID="$2h"

# set which model was fit
PRF_MODEL=Iso2DGaussianModel_LinearDeveining
#Iso2DGaussianModel
#DoG_Iso2DGaussianModel

# set desired surface projection params
PROJ_FRAC=0.5
FWHM=0.5

INTERP_ANTS=NearestNeighbor
#BSpline[5]
#NearestNeighbor

INTERP_VOL2SURF=nearest
#trilin
#nearest

PRFPY_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/prfpy_${PRF_MODEL}/$SUB_ID
REF_IMG=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/func/meanFunc/_subject_id_$SUB_ID/merged_func_mcf.nii_mean_reg.nii
#REF_IMG=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/anat/_subject_id_$SUB_ID/UNI_corrected.nii
TRANSFORM_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/manualcorr/coreg/_subject_id_$SUB_ID/func2struct

if [ "$SUB_ID" = "sub-04" ]; then
    SUBJECTS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUB_ID/autorecon_pial_rerun
else
    SUBJECTS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_$SUB_ID/autorecon_pial_rerun
fi

# change directory to prfpy outputs
cd $PRFPY_DIR

# apply transform and surface project mean functional
MEAN_FUNC=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/func/meanFunc/_subject_id_$SUB_ID/merged_func_mcf.nii_mean_reg.nii
antsApplyTransforms --default-value 0 -d 3 \
	    --input $MEAN_FUNC \
        --interpolation $INTERP_ANTS \
        --output $PRFPY_DIR/reg_meanFunc.nii \
        --reference-image $REF_IMG \
        --transform ${TRANSFORM_DIR}/regFUNC2UNI_Composite.h5 \
        --transform ${TRANSFORM_DIR}/regUNI2T1_Composite.h5
mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.meanFunc.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_meanFunc.nii --fwhm $FWHM
fslreorient2std $PRFPY_DIR/reg_meanFunc.nii $PRFPY_DIR/reg_meanFunc_fsl.nii

# # apply transform and surface project bar data
# antsApplyTransforms --default-value 0 \
#         --input-image-type 3 \
#         --input $PRFPY_DIR/${HEM_ID}_bar.nii \
#         --interpolation $INTERP_ANTS \
#         --output $PRFPY_DIR/reg_${HEM_ID}_bar.nii \
#         --reference-image $REF_IMG \
#         --transform ${TRANSFORM_DIR}/regFUNC2UNI_Composite.h5 \
#         --transform ${TRANSFORM_DIR}/regUNI2T1_Composite.h5
# mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.bar.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_bar.nii --fwhm $FWHM
# fslreorient2std $PRFPY_DIR/reg_${HEM_ID}_bar.nii $PRFPY_DIR/reg_${HEM_ID}_bar_fsl.nii

# apply transform and surface project pRF parameters
PARAM=ecc
antsApplyTransforms --default-value 0 \
        -d 3 \
	--input $PRFPY_DIR/${HEM_ID}_$PARAM.nii \
        --interpolation $INTERP_ANTS \
        --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
        --reference-image $REF_IMG \
        --transform ${TRANSFORM_DIR}/regFUNC2UNI_Composite.h5 \
        --transform ${TRANSFORM_DIR}/regUNI2T1_Composite.h5
mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii --fwhm $FWHM

PARAM=pol
antsApplyTransforms --default-value 0 \
        -d 3 \
	--input $PRFPY_DIR/${HEM_ID}_$PARAM.nii \
        --interpolation $INTERP_ANTS \
        --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
        --reference-image $REF_IMG \
        --transform ${TRANSFORM_DIR}/regFUNC2UNI_Composite.h5 \
        --transform ${TRANSFORM_DIR}/regUNI2T1_Composite.h5
mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii --fwhm $FWHM

PARAM=rsq
antsApplyTransforms --default-value 0 \
        -d 3 \
	--input $PRFPY_DIR/${HEM_ID}_$PARAM.nii \
        --interpolation $INTERP_ANTS \
        --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
        --reference-image $REF_IMG \
        --transform ${TRANSFORM_DIR}/regFUNC2UNI_Composite.h5 \
        --transform ${TRANSFORM_DIR}/regUNI2T1_Composite.h5
mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii --fwhm $FWHM

# pRF, sRF
if [ "$PRF_MODEL" = "Iso2DGaussianModel" ] || [ "$PRF_MODEL" = "Iso2DGaussianModel_LinearDeveining" ]; then
    PARAM=sigma
    antsApplyTransforms --default-value 0 \
            -d 3 \
            --input $PRFPY_DIR/${HEM_ID}_$PARAM.nii \
            --interpolation $INTERP_ANTS \
            --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
            --reference-image $REF_IMG \
            --transform ${TRANSFORM_DIR}/regFUNC2UNI_Composite.h5 \
            --transform ${TRANSFORM_DIR}/regUNI2T1_Composite.h5
    mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii --fwhm $FWHM
else
    PARAM=prf_size
    antsApplyTransforms --default-value 0 \
            -d 3 \
            --input $PRFPY_DIR/${HEM_ID}_$PARAM.nii \
            --interpolation $INTERP_ANTS \
            --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
            --reference-image $REF_IMG \
            --transform ${TRANSFORM_DIR}/regFUNC2UNI_Composite.h5 \
            --transform ${TRANSFORM_DIR}/regUNI2T1_Composite.h5
    mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii --fwhm $FWHM

    PARAM=srf_size
    antsApplyTransforms --default-value 0 \
            -d 3 \
            --input $PRFPY_DIR/${HEM_ID}_$PARAM.nii \
            --interpolation $INTERP_ANTS \
            --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
            --reference-image $REF_IMG \
            --transform ${TRANSFORM_DIR}/regFUNC2UNI_Composite.h5 \
            --transform ${TRANSFORM_DIR}/regUNI2T1_Composite.h5
    mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii --fwhm $FWHM


    PARAM=prf_amp
    antsApplyTransforms --default-value 0 \
            -d 3 \
            --input $PRFPY_DIR/${HEM_ID}_$PARAM.nii \
            --interpolation $INTERP_ANTS \
            --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
            --reference-image $REF_IMG \
            --transform ${TRANSFORM_DIR}/regFUNC2UNI_Composite.h5 \
            --transform ${TRANSFORM_DIR}/regUNI2T1_Composite.h5
    mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii --fwhm $FWHM

    PARAM=srf_amp
    antsApplyTransforms --default-value 0 \
            -d 3 \
            --input $PRFPY_DIR/${HEM_ID}_$PARAM.nii \
            --interpolation $INTERP_ANTS \
            --output $PRFPY_DIR/reg_${HEM_ID}_$PARAM.nii \
            --reference-image $REF_IMG \
            --transform ${TRANSFORM_DIR}/regFUNC2UNI_Composite.h5 \
            --transform ${TRANSFORM_DIR}/regUNI2T1_Composite.h5
    mri_vol2surf --hemi ${HEM_ID} --interp $INTERP_VOL2SURF --o ${HEM_ID}.$PARAM.mgh --out_type mgh --regheader $SUB_ID --projfrac $PROJ_FRAC --mov reg_${HEM_ID}_$PARAM.nii --fwhm $FWHM
fi