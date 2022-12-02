#!/bin/bash

#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --job-name=p00d02
#SBATCH --mail-type=end
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=30
#SBATCH --time=2-00:00:00
#SBATCH --qos=standard

sub=1
depth=2

#Extract parameters
echo sub $sub
echo depth $depth
echo using $SLURM_CPUS_PER_TASK cpus

### Set up runtime environment
module add MATLAB/2021a

# wait a bit so it doesn't crash
sleep 50

cd /home/mayaaj90/scripts/project-00-7t-pipeline-dev/

### Start job
matlab -nosplash -noFigureWindows -r "pRFmapping_curta(${sub},${depth},$SLURM_CPUS_PER_TASK)" > serial.out #this worked
echo set to run
### Output core and memory efficiency

seff $SLURM_JOBID
