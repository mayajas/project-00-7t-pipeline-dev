SUBJECT_ID="sub-0$1"
HEMI="$2h"

echo "Converting manually-defined ROI labels to niftis in functional space for subject $SUBJECT_ID"

# define directories
if [ "$SUBJECT_ID" = "sub-04" ]; then
    SUBJECTS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_$SUBJECT_ID/autorecon_pial_rerun
else
    SUBJECTS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial_rerun
fi
ROI_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/prfpy/$SUBJECT_ID/ROIs
COREG_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/manualcorr/coreg/_subject_id_$SUBJECT_ID

# define anatomical and func images
ANAT_IMG=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/prfpy/$SUBJECT_ID/UNI_corrected.nii
FUNC_IMG=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/prfpy/$SUBJECT_ID/meanFunc.nii

# convert manually-defined label to nifti (anat space)
mri_label2vol --label $ROI_DIR/${HEMI}_V1.label --temp $ANAT_IMG \
--o $ROI_DIR/${HEMI}_V1.nii --subject $SUBJECT_ID --hemi ${HEMI} --fill-ribbon --regheader $ANAT_IMG

#--proj frac -0.5 1.5 0.01 

# apply inverse transform to func space
${ANTSPATH}antsApplyTransforms \
  -d 3 \
  -i $ROI_DIR/${HEMI}_V1.nii \
  -r $FUNC_IMG \
  -t [ $COREG_DIR/registered_0GenericAffine.mat, 1 ] \
  -o $ROI_DIR/func_${HEMI}_V1.nii \
  --interpolation NearestNeighbor
