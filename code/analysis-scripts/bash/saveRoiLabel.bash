SUBJECT_ID=sub-01
SUBJECTS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial_rerun

# convert manually-defined label to nifti (anat space)
mri_label2vol --label lh_V1.label --temp $SUBJECTS_DIR/$SUBJECT_ID/mri/orig.mgz --o lh_V1.nii --proj frac 0 0 0 --subject $SUBJECT_ID --hemi lh --identity

# apply inverse transform to func space
REF_IMG=meanFunc.nii
COREG_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_laminar_fMRI_func_pRF_realign2mean/_subject_id_$SUBJECT_ID/coreg

antsApplyTransforms --default-value 0 --float 0 --input lh_V1.nii --interpolation NearestNeighbor \
--output func_lh_V1.nii \
--reference-image $REF_IMG \
--transform [ $COREG_DIR/registered_0GenericAffine.mat, 1 ] \
--transform $COREG_DIR/registered_1InverseWarp.nii.gz \
--verbose





SUBJECT_ID=sub-01
SUBJECTS_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_$SUBJECT_ID/autorecon_pial_rerun

# apply inverse transform to func space
REF_IMG=meanFunc.nii
COREG_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/derivatives/wf_laminar_fMRI_func_pRF_realign2mean/_subject_id_$SUBJECT_ID/coreg

antsApplyTransforms --default-value 0 --float 0 --input UNI_corrected.nii --interpolation BSpline[5] \
--output func_UNI_corrected.nii \
--reference-image $REF_IMG \
--transform [ $COREG_DIR/registered_0GenericAffine.mat, 1 ] \
--transform $COREG_DIR/registered_1InverseWarp.nii.gz \
--verbose