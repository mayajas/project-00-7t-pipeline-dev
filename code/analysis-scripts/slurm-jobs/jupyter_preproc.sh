#!/bin/bash

#SBATCH --job-name=jn
#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --mail-type=none
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=200M
#SBATCH --cpus-per-task=1
#SBATCH --time=0:05:00
#SBATCH --qos=prio

PORT=8888
NODE=$(hostname -s)

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

echo  "Open ssh connection"
echo "NODE=$NODE"
ssh -N -f -R $PORT:localhost:$PORT login

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "CPUs per task=$OMP_NUM_THREADS"
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$SLURM_CPUS_PER_TASK

# this variable needs to be unset (see: https://github.com/jupyter/notebook/issues/1318)
unset XDG_RUNTIME_DIR

# change directory to project dir
cd ~/projects/project-00-7t-pipeline-dev/code/

# run jupyter notebook
jupyter notebook --no-browser --port $PORT --ip=$(NODE)
