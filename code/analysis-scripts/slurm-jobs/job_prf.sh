#!/bin/bash

#SBATCH --job-name=laminar_fMRI_pipeline
#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --mail-type=none
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --qos=standard

# print tunneling instructions jupyter-log
# activate py36 environment
# include for this reason: https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
eval "$(conda shell.bash hook)"
conda activate py36_prf

# get environment variable of num threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$SLURM_CPUS_PER_TASK
echo "CPUs per task=$OMP_NUM_THREADS"

# Run pyprf analysis
pyprf -config /scratch/mayaaj90/7T-ret/pRF_config/config_project-00-7t-pipeline-dev.csv
