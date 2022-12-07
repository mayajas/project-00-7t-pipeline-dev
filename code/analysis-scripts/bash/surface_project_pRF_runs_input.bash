#!/bin/bash

SUBJECT=sub-01
which_coreg=midocc

project=project-00-7t-pipeline-dev

if [ "$SUBJECT" = "sub-04" ]; then
    wf_name=wf_advanced_skullstrip_sub-04
else
    wf_name=wf_advanced_skullstrip
fi

FSDIR=/home/mayajas/scratch/${project}/derivatives/$wf_name/_subject_id_${SUBJECT}/autorecon_pial
OUTDIR=/home/mayajas/scratch/${project}/pRF/data
REG_MAT=/home/mayajas/scratch/${project}/output/func/coreg_$which_coreg/_subject_id_${SUBJECT}/registered_0GenericAffine.mat
REG_NIIGZ=/home/mayajas/scratch/${project}/output/func/coreg_$which_coreg/_subject_id_${SUBJECT}/registered_1Warp.nii.gz
MEAN_FUNC=/home/mayajas/scratch/${project}/output/func/meanFunc/_subject_id_${SUBJECT}/merged_func_mcf.nii_mean_reg.nii
PRF_BAR1=/home/mayajas/scratch/${project}/output/func/sliceTimeCorr/_subject_id_${SUBJECT}/_sess_id_task-bar_run-01_sess_nr_0_sess_nvol_124/atask-bar_run-01_roi_warp4D.nii
PRF_BAR2=/home/mayajas/scratch/${project}/output/func/sliceTimeCorr/_subject_id_${SUBJECT}/_sess_id_task-bar_run-02_sess_nr_1_sess_nvol_124/atask-bar_run-02_roi_warp4D.nii
PRF_ECC1=/home/mayajas/scratch/${project}/output/func/sliceTimeCorr/_subject_id_${SUBJECT}/_sess_id_task-ecc_run-01_sess_nr_4_sess_nvol_131/atask-ecc_run-01_roi_warp4D.nii
PRF_ECC2=/home/mayajas/scratch/${project}/output/func/sliceTimeCorr/_subject_id_${SUBJECT}/_sess_id_task-ecc_run-02_sess_nr_4_sess_nvol_131/atask-ecc_run-02_roi_warp4D.nii
PRF_POL1=/home/mayajas/scratch/${project}/output/func/sliceTimeCorr/_subject_id_${SUBJECT}/_sess_id_task-pol_run-01_sess_nr_2_sess_nvol_254/atask-pol_run-01_roi_warp4D.nii
PRF_POL2=/home/mayajas/scratch/${project}/output/func/sliceTimeCorr/_subject_id_${SUBJECT}/_sess_id_task-pol_run-02_sess_nr_2_sess_nvol_254/atask-pol_run-02_roi_warp4D.nii
REST=/home/mayajas/scratch/${project}/output/func/sliceTimeCorr/_subject_id_${SUBJECT}/_sess_id_rest_run-01_sess_nr_6_sess_nvol_50/arest_run-01_roi_warp4D.nii
OCC=/home/mayajas/scratch/${project}/manualcorr/func/_subject_id_${SUBJECT}/occipitalMask.nii

if [[ ! -d $OUTDIR/${SUBJECT} ]]
then
	mkdir -p $OUTDIR/${SUBJECT}
fi
./surface_project_pRF_runs.bash $SUBJECT $FSDIR $OUTDIR $REG_MAT $MEAN_FUNC $PRF_BAR1 $PRF_BAR2 $PRF_ECC1 $PRF_ECC2 $PRF_POL1 $PRF_POL2 $REST $OCC

#PRF_FSDIR=/home/mayajas/scratch/${project}/pRF/data_FS
#if [[ ! -d $PRF_FSDIR/${SUBJECT} ]]
#then
#	mkdir -p $PRF_FSDIR/${SUBJECT}
#	cp -r $FSDIR/${SUBJECT} $PRF_FSDIR
#fi
