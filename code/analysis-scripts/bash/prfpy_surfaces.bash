SUB_ID=sub-01




cd /home/mayajas/scratch/project-00-7t-pipeline-dev/output/prfpy/$SUB_ID
SUBJECTS_DIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_$SUB_ID/autorecon_pial_rerun

VAR=sigma
mri_vol2surf --hemi lh --interp nearest --o lh.$VAR.mgh --out_type mgh --regheader $SUB_ID --projfrac 0.5 --mov reg_lh_$VAR.nii
VAR=pol
mri_vol2surf --hemi lh --interp nearest --o lh.$VAR.mgh --out_type mgh --regheader $SUB_ID --projfrac 0.5 --mov reg_lh_$VAR.nii
VAR=rsq
mri_vol2surf --hemi lh --interp nearest --o lh.$VAR.mgh --out_type mgh --regheader $SUB_ID --projfrac 0.5 --mov reg_lh_$VAR.nii
VAR=ecc
mri_vol2surf --hemi lh --interp nearest --o lh.$VAR.mgh --out_type mgh --regheader $SUB_ID --projfrac 0.5 --mov reg_lh_$VAR.nii