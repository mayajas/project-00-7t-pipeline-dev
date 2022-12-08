#!/bin/bash

#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --job-name=p00d02
#SBATCH --mail-type=end
#SBATCH --mem=1GB
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --qos=standard

project="project-00-7t-pipeline-dev"

### Set up runtime environment
module add MATLAB/2021a

# wait a bit so it doesn't crash
sleep 50

cd /home/mayaaj90/scripts/project-00-7t-pipeline-dev/

### Start job
matlab -nosplash -noFigureWindows -r "getPRFsizes_curta(project)"
echo set to run
### Output core and memory efficiency

seff $SLURM_JOBID
