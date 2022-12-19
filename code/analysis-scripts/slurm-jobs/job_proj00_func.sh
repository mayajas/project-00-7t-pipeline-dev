#!/bin/bash

#SBATCH --job-name=p00-func
#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --qos=prio

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

# get environment variable of num threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# unset Java max heap space option
unset _JAVA_OPTIONS

# run the thing
cd /home/mayaaj90/projects/project-00-7t-pipeline-dev/code/analysis-scripts/python
ipython project-00-7t-pipeline-dev-functional-pRF.py > log_p00_func.txt
