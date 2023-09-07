#!/bin/bash

#SBATCH -o PATH_TO_LOGS/heat_%a.slurm.out
#SBATCH -e PATH_TO_LOGS/heat_%a.slurm.err
#SBATCH -n N_TASKS
#SBATCH --gres=gpu:N_GPUS
#SBATCH --cpus-per-gpu=N_CPUS
#SBATCH --mem MEMORY

export MEZEHOME=PATH_TO_MEZE

LIG_NUMBER=$SLURM_ARRAY_TASK_ID

python $MEZEHOME/equilibrate.py $LIG_NUMBER INPUTFILE
