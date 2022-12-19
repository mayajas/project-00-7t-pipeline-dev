#!/bin/bash

#SBATCH --mail-user=mayaaj90@zedat.fu-berlin.de
#SBATCH --job-name=p0prf
#SBATCH --mem=10G
#SBATCH --cpus-per-task=20
#SBATCH --time=1-00:00:00
#SBATCH --qos=standard
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

declare -a combinations
index=0
for sub in 1 2 3 4
do
    for depth in 1 2 3 4 5 6 7 8 9
    do
      	combinations[$index]="$sub $depth"
        index=$((index + 1))
    done
done

parameters=(${combinations[${SLURM_ARRAY_TASK_ID}]})

sub=${parameters[0]}
depth=${parameters[1]}

#Print out parameters
echo sub $sub
echo depth $depth
echo using $SLURM_CPUS_PER_TASK cpus

### Set up runtime environment
module add MATLAB/2021a

# wait a bit so it doesn't crash
sleep 50

# change directory to the one containing the scripts
cd /home/mayaaj90/projects/project-00-7t-pipeline-dev/code/analysis-scripts/matlab

### Start job
echo set to run
matlab -nosplash -noFigureWindows -r "pRFmapping_curta($sub,$depth,$SLURM_CPUS_PER_TASK)" > serial_sub${sub}_depth${depth}.out #this worked
