#!/bin/bash

# e.g., ./surface_project_pRF_runs.sh sub-01 /usr/local/freesurfer/subjects /home/mayajas/Documents/project-00-7t-pipeline-dev/data/pRF/data /home/mayajas/Documents/project-00-7t-pipeline-dev/data/output/func/coreg/_subject_id_sub-01/registered_0GenericAffine.mat /home/mayajas/Documents/project-00-7t-pipeline-dev/data/output/func/sliceTimeCorr/_subject_id_sub-01/_sess_id_task-bar_run-01_sess_nr_0_sess_nvol_124/atask-bar_run-01_roi_warp4D.nii /home/mayajas/Documents/project-00-7t-pipeline-dev/data/output/func/sliceTimeCorr/_subject_id_sub-01/_sess_id_task-bar_run-02_sess_nr_1_sess_nvol_124/atask-bar_run-02_roi_warp4D.nii /home/mayajas/scratch/project-00-7t-pipeline-dev/manualcorr/sub-01/occipital.nii

SUBJECT=$1
FSDIR=$2
OUTDIR=$3
REG_NIIGZ=$4
MEAN_FUNC=$5
PRF_BAR1=$6
PRF_BAR2=$7
PRF_ECC1=$8
PRF_ECC2=$9
PRF_POL1=${10}
PRF_POL2=${11}
REST=${12}
OCC=${13}

export SUBJECTS_DIR=$FSDIR
ANAT=$FSDIR/$SUBJECT/mri/T1.mgz
cd $OUTDIR/$SUBJECT

echo "Occipital file: ${OCC}"
echo "Bar1 file: ${PRF_BAR1}"
echo "Bar2 file: ${PRF_BAR2}"

# Convert ANTS/ITK nonlinear transform to .m3z (FreeSurfer):
# convert non-linear deformation field warp file formats.
mri_warp_convert --initk $REG_NIIGZ --outm3z $OUTDIR/$SUBJECT/registered_1Warp.m3z --insrcgeom $MEAN_FUNC
#mri_warp_convert --inm3z registered_1Warp.m3sz --outitk out.nii.gz

# Convert ANTS/ITK transform to LTA (FreeSurfer):
# First convert the ANTS binary mat file to ITK text file format and then to lta (adding src and trg geometry info, from images that were used to create the transform in ANTS):
$ANTSPATH/ConvertTransformFile 3 $REG_MAT $OUTDIR/$SUBJECT/registered_0GenericAffine.txt
lta_convert --initk $OUTDIR/$SUBJECT/registered_0GenericAffine.txt --outlta $OUTDIR/$SUBJECT/registered_0GenericAffine.lta --src $MEAN_FUNC --trg $ANAT --subject $SUBJECT

# Coregister functional runs to structural space (volume)
# mean functional
mri_vol2vol --mov $MEAN_FUNC --lta $OUTDIR/$SUBJECT/registered_0GenericAffine.lta --o $OUTDIR/$SUBJECT/reg_meanFunc_intermediary.nii --trilin --fstarg --no-resample
# the above works, i.e., gives the same result as 
# antsApplyTransforms --interpolation BSpline[5] -d 3 -i $MEAN_FUNC -r $ANAT -t $REG_MAT -o registered_applied_lin.nii
mri_vol2vol --mov $OUTDIR/$SUBJECT/reg_meanFunc_intermediary.nii --m3z $OUTDIR/$SUBJECT/registered_1Warp.m3z --noDefM3zPath --fstarg --o $OUTDIR/$SUBJECT/reg_meanFunc.nii --trilin --no-resample
# this on the other hand...does not work

## FS
mri_vol2vol --mov $MEAN_FUNC --m3z $OUTDIR/$SUBJECT/registered_1Warp.m3z --lta $OUTDIR/$SUBJECT/registered_0GenericAffine.lta --noDefM3zPath --o $OUTDIR/$SUBJECT/reg_meanFunc_fs_nonlin_lin.nii --trilin --fstarg --no-resample
mri_vol2vol --mov $MEAN_FUNC --lta $OUTDIR/$SUBJECT/registered_0GenericAffine.lta --m3z $OUTDIR/$SUBJECT/registered_1Warp.m3z --noDefM3zPath --o $OUTDIR/$SUBJECT/reg_meanFunc_fs_lin_nonlin.nii --trilin --fstarg --no-resample

## ANTS
antsApplyTransforms --interpolation BSpline[5] -d 3 -i $MEAN_FUNC -r $ANAT -t $REG_NIIGZ -t $REG_MAT -o reg_meanFunc_ants_nonlin_lin.nii
antsApplyTransforms --interpolation BSpline[5] -d 3 -i $MEAN_FUNC -r $ANAT -t $REG_MAT -t $REG_NIIGZ -o reg_meanFunc_ants_lin_nonlin.nii

## FS with intermediary step
mri_vol2vol --mov $MEAN_FUNC --lta $OUTDIR/$SUBJECT/registered_0GenericAffine.lta --o $OUTDIR/$SUBJECT/reg_meanFunc_fs_intermediary.nii --trilin --targ $ANAT --no-resample
mri_vol2vol --mov $OUTDIR/$SUBJECT/reg_meanFunc_fs_intermediary.nii --m3z $OUTDIR/$SUBJECT/registered_1Warp.m3z --noDefM3zPath --o $OUTDIR/$SUBJECT/reg_meanFunc_fs_composite_lin_nonline.nii --trilin --targ $ANAT --no-resample

# bar sess 1
mri_vol2vol --mov $PRF_BAR1 --lta $OUTDIR/$SUBJECT/registered_0GenericAffine.lta --o $OUTDIR/$SUBJECT/reg_std-prf_bar1.nii --trilin --fstarg --no-resample
# bar sess 2
mri_vol2vol --mov $PRF_BAR2 --lta $OUTDIR/$SUBJECT/registered_0GenericAffine.lta --o $OUTDIR/$SUBJECT/reg_std-prf_bar2.nii --trilin --fstarg --no-resample
# occipital mask
mri_vol2vol --mov $OCC --lta $OUTDIR/$SUBJECT/registered_0GenericAffine.lta --o $OUTDIR/$SUBJECT/reg_occ.nii --trilin --fstarg --no-resample


# Coregistration and surface sampling across depths
depthTxt=("m1.5" "m1.0" "m0.5" "0.0" "0.5" "1.0" "1.5" "2.0" "2.5")
depth=(-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0 2.5)

for i in $(seq 0 1 $((${#depth[@]}-1)))
do
	echo "Depth: ${depth[$i]}"
	for hemis in lh rh
	do
		method=trilinear
		# bar sess 1
		mri_vol2surf --mov $PRF_BAR1 --srcreg $OUTDIR/$SUBJECT/registered_0GenericAffine.lta --hemi ${hemis} --cortex \
		--out $OUTDIR/$SUBJECT/${hemis}_bar_sess1_depth${depthTxt[$i]}.mgh --interp ${method} --projdist ${depth[$i]}
		# bar sess 2
		mri_vol2surf --mov $PRF_BAR2 --srcreg $OUTDIR/$SUBJECT/registered_0GenericAffine.lta --hemi ${hemis} --cortex \
		--out $OUTDIR/$SUBJECT/${hemis}_bar_sess2_depth${depthTxt[$i]}.mgh --interp ${method} --projdist ${depth[$i]}
		
		method=nearest
		# occipital label
		mri_vol2surf --mov $OCC --srcreg $OUTDIR/$SUBJECT/registered_0GenericAffine.lta --hemi ${hemis} --cortex \
		--out $OUTDIR/$SUBJECT/${hemis}_occ_depth${depthTxt[$i]}.mgh --interp ${method} --projdist ${depth[$i]}
		mri_vol2label --i $OUTDIR/$SUBJECT/${hemis}_occ_depth${depthTxt[$i]}.mgh --id 1  --surf ${SUBJECT} ${hemis}  \
		--l $OUTDIR/$SUBJECT/${hemis}_occ_depth${depthTxt[$i]}.label
	done
done
