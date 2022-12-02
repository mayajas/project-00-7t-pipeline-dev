#!/bin/bash

#SBATCH --job-name=laminar_fMRI_pipeline
#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --mail-type=none
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --qos=standard

XDG_RUNTIME_DIR=""
node=$(hostname -s)
user=$(whoami)
cluster="curta"
port=8888

# print tunneling instructions jupyter-log
echo -e "
Command to create ssh tunnel:
ssh -N -f -L ${port}:${node}:${port} ${user}@${cluster}.zedat.fu-berlin.de
Use a Browser on your local machine to go to:
localhost:${port}  (prefix w/ https:// if using password)
"

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
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$SLURM_CPUS_PER_TASK
echo "CPUs per task=$OMP_NUM_THREADS"
#export available_mem=$SLURM_MEM_PER_CPU
#echo "Available memory=$available_mem"

# Run Jupyter Notebook
jupyter notebook --no-browser --port=${port} --ip=${node} --NotebookApp.max_buffer_size=1000000000
