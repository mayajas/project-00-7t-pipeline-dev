#!/bin/bash

#SBATCH --job-name=prfpy
#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --mail-type=none
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=20
#SBATCH --time=8:00:00
#SBATCH --qos=prio

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

SUB_ID=${parameters[0]}
HEM_ID=${parameters[1]}
#
#SUB_ID=0
#HEM_ID=1
# 0 - lh, 1 - rh

#Print out parameters
echo sub $SUB_ID
echo hem $HEM_ID

# activate py36 environment
# include for this reason: https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
eval "$(conda shell.bash hook)"
conda activate py38prf

# load modules
module purge
module load MATLAB/2021a
module load FreeSurfer/dev-centos7_x86_64 ; source $FREESURFER_HOME/SetUpFreeSurfer.sh
module load FSL/5.0.11-foss-2018b-Python-3.6.6
module add ANTs/2.3.1-foss-2018b-Python-3.6.6
module add AFNI/18.3.00-foss-2018b-Python-3.6.6

# get environment variable of num threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# unset Java max heap space option
unset _JAVA_OPTIONS

# run the thing
cd /home/mayaaj90/projects/project-00-7t-pipeline-dev/code/analysis-scripts/python
#ipython pRF-mapping.py $SUB_ID $HEM_ID > log_prf_sub${SUB_ID}_hem${HEM_ID}.txt
ipython pRF_mapping_DoG_Iso2DGaussian.py $SUB_ID $HEM_ID > log_prf_sub${SUB_ID}_hem${HEM_ID}.txt
