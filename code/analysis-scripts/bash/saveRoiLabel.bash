SUBJECT_ID="sub-0$1"
HEMI="$2h"

PRF_MODEL=prfpy_vol_fit_hrf_False

echo "Converting manually-defined ROI labels to niftis in functional space for subject $SUBJECT_ID"

# define directories
PROJECT_ID=project-00-7t-pipeline-dev
PROJECT_DIR=/home/mayajas/scratch/$PROJECT_ID
SUB_DIR=$PROJECT_DIR/output/${PRF_MODEL}/$SUBJECT_ID

if [ "$SUBJECT_ID" = "sub-04" ]; then
    SUBJECTS_DIR=$PROJECT_DIR/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUBJECT_ID/autorecon_pial
else
    SUBJECTS_DIR=$PROJECT_DIR/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial
fi
COREG_DIR=$PROJECT_DIR/manualcorr/coreg/_subject_id_$SUBJECT_ID/struct2func

# define anatomical and func images
# ANAT_IMG=$PROJECT_DIR/output/${PRF_MODEL}/$SUBJECT_ID/UNI_corrected.nii
ANAT_IMG=$PROJECT_DIR/output/func/anat/_subject_id_${SUBJECT_ID}/T1_out.nii
FUNC_IMG=$PROJECT_DIR/output/${PRF_MODEL}/$SUBJECT_ID/meanFunc.nii


declare -a ROIs=("V1" "V2" "V3" "V4" "V1d" "V1v" "V2d" "V2v" "V3d" "V3v")

for ROI in ${ROIs[@]}; do
    ## ROI label 
    # convert manually-defined label to nifti (anat space)
    mri_label2vol --label $SUB_DIR/${HEMI}.$ROI.label --temp $ANAT_IMG \
    --o $SUB_DIR/${HEMI}.$ROI.nii --subject $SUBJECT_ID --hemi ${HEMI} --identity --proj frac -0.5 1.5 0.01

    # apply inverse transform to func space
    antsApplyTransforms --default-value 0 \
        -d 3 \
        --input $SUB_DIR/${HEMI}.$ROI.nii \
        --interpolation NearestNeighbor \
        --output $SUB_DIR/func_${HEMI}.$ROI.nii \
        --reference-image $FUNC_IMG \
        --transform $COREG_DIR/regT12UNI_Composite.h5 \
        --transform $COREG_DIR/regUNI2FUNC_Composite.h5
done

## 2deg iso-ecc label 
# convert manually-defined label to nifti (anat space)
mri_label2vol --label $SUB_DIR/${HEMI}.ecc2.label --temp $ANAT_IMG \
--o $SUB_DIR/${HEMI}.ecc2.nii --subject $SUBJECT_ID --hemi ${HEMI} --identity --proj frac -0.5 1.5 0.01

# apply inverse transform to func space
antsApplyTransforms --default-value 0 \
    -d 3 \
    --input $SUB_DIR/${HEMI}.ecc2.nii \
    --interpolation NearestNeighbor \
    --output $SUB_DIR/func_${HEMI}.ecc2.nii \
    --reference-image $FUNC_IMG \
    --transform $COREG_DIR/regT12UNI_Composite.h5 \
    --transform $COREG_DIR/regUNI2FUNC_Composite.h5