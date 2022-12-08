#!/bin/bash

#SBATCH --job-name=p00-anat
#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=192G
#SBATCH --cpus-per-task=2
#SBATCH --time=40:00:00
#SBATCH --qos=standard

# activate py36 environment
# include for this reason: https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
eval "$(conda shell.bash hook)"
conda activate py36

# load modules
module purge
module load MATLAB/2021a
module load FreeSurfer/dev-centos7_x86_64 ; source $FREESURFER_HOME/SetUpFreeSurfer.sh
module load FSL/5.0.11-foss-2018b-Python-3.6.6
module add ANTs/2.3.1-foss-2018b-Python-3.6.6
module add AFNI/18.3.00-foss-2018b-Python-3.6.6

# set environment variable of num threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "Available threads: $OMP_NUM_THREADS"

unset _JAVA_OPTIONS
ipython project-00-7t-pipeline-dev-advanced-skullstrip.py > log_p00_anat.txt
