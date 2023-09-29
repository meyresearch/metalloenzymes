#!/bin/bash

#SBATCH -o PATH_TO_LOGS/afe_%a.slurm.out
#SBATCH -e PATH_TO_LOGS/afe_%a.slurm.err
#SBATCH -n N_TASKS
#SBATCH --gres=gpu:N_GPUS
#SBATCH --cpus-per-gpu=N_CPUS
#SBATCH --mem MEMORY

export MEZEHOME=PATH_TO_MEZE

start=`date +%s`

python $MEZEHOME/meze.py PROTOCOLFILE

end=`date +%s`
runtime=$((end - start))

echo "Finished in $runtime seconds, or $((runtime/60)) minutes"