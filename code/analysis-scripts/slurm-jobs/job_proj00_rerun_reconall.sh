#!/bin/bash

#SBATCH --job-name=fsmaned
#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=5
#SBATCH --time=5:00:00
#SBATCH --qos=standard

# set subject id
SUBJECTS=(sub-01 sub-02 sub-03 sub-04)

SUBJECT_ID=${SUBJECTS[${SLURM_ARRAY_TASK_ID}]}
#sub-01

# white matter or pial edits?
manual_edits=pial

echo "Running $SUBJECT_ID $manual_edits surface edits."

# activate py36 environment
# include for this reason: https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
eval "$(conda shell.bash hook)"
conda activate py36

# load modules
module purge
module load FreeSurfer/dev-centos7_x86_64 ; source $FREESURFER_HOME/SetUpFreeSurfer.sh

# set subject directory and output directory for outputing new surfaces
if [ "$SUBJECT_ID" = "sub-04" ]; then
    ORIG_FS_DIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_${SUBJECT_ID}/autorecon_pial
    SUBJECTS_DIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip_sub-04/_subject_id_${SUBJECT_ID}/autorecon_pial_rerun
else
    ORIG_FS_DIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_${SUBJECT_ID}/autorecon_pial
    SUBJECTS_DIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/derivatives/wf_advanced_skullstrip/_subject_id_${SUBJECT_ID}/autorecon_pial_rerun
fi

OUTDIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/output/freesurfer_reregenerated/_subject_id_$SUBJECT_ID/
MANUALEDDIR=/scratch/mayaaj90/project-00-7t-pipeline-dev/manualcorr/freesurfer_externalbrainmask/_subject_id_${SUBJECT_ID}

# copy FS folder if doesn't already exist
if [ ! -d "$SUBJECTS_DIR" ]; then
    echo "Copying previous FS directory."
    cp -r $ORIG_FS_DIR $SUBJECTS_DIR/
else
    echo "New FS directory already exists."
fi

# run recon-all
if [ "$manual_edits" = "pial" ]; then
    echo "Running pial surface reconstruction."
    cp $MANUALEDDIR/brain.finalsurfs.mgz $SUBJECTS_DIR/$SUBJECT_ID/mri/ --verbose
    echo "Copied manually-edited brain.finalsurfs.mgz file"
    recon-all -autorecon-pial -no-isrunning -hires -parallel -nomaskbfs -subjid $SUBJECT_ID -sd $SUBJECTS_DIR
    # -nomotioncor -notalairach -nonuintensitycor -nonormalization -noskullstrip -nogcareg -nocanorm -nocareg -nocalabel -nonormalization2 -nomaskbfs -nosegmentation -nofill -notessellate -nosmooth1 -noinflate1 -noqsphere -nofix -nowhite -nosmooth2 -noinflate2 -nocurvHK -nocurvstats -nosphere -nosurfreg -nojacobian_white -noavgcurv -nocortparc -noparcstats -nocortparc2 -noparcstats2 -nocortparc3 -noparcstats3 -nopctsurfcon -nocortribbon -nohyporelabel -noaparc2aseg -noapas2aseg -nosegstats -nowmparc -nobalabels
elif [ "$manual_edits" = "white" ]; then
    echo "Running WM surface reconstruction."
    cp $MANUALEDDIR/wm.mgz $SUBJECTS_DIR/$SUBJECT_ID/mri/ --verbose
    echo "Copied manually-edited wm.mgz file"
    recon-all -autorecon2-wm -autorecon3 -no-isrunning -hires -parallel -subjid $SUBJECT_ID -sd $SUBJECTS_DIR
fi

# copy surface outfiles
if [ ! -d "$OUTDIR" ]; then
    echo "Making outdir: $OUTDIR."
    mkdir $OUTDIR
fi
if [ "$manual_edits" = "pial" ] || [ "$manual_edits" = "white" ]; then
    cp $SUBJECTS_DIR/$SUBJECT_ID/surf/rh.pial $OUTDIR
    cp $SUBJECTS_DIR/$SUBJECT_ID/surf/rh.white $OUTDIR
    cp $SUBJECTS_DIR/$SUBJECT_ID/surf/lh.pial $OUTDIR
    cp $SUBJECTS_DIR/$SUBJECT_ID/surf/lh.white $OUTDIR
fi
