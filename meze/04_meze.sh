#!/bin/bash

#SBATCH -o PATH_TO_LOGS/meze_%a.slurm.out
#SBATCH -e PATH_TO_LOGS/meze_%a.slurm.err
#SBATCH -n N_TASKS
#SBATCH --gres=gpu:N_GPUS
#SBATCH --cpus-per-gpu=N_CPUS

eval "$(conda shell.bash hook)"
conda activate obss


export MEZEHOME=PATH_TO_MEZE

start=`date +%s`

id=$SLURM_ARRAY_TASK_ID

source $MEZEHOME/parse.sh
transformation=${transformations_array[$id]}

python $MEZEHOME/meze.py PROTOCOLFILE $transformation

end=`date +%s`
runtime=$((end - start))

echo "Finished in $runtime seconds, or $((runtime/60)) minutes"
