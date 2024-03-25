#!/bin/bash

#SBATCH -o PATH_TO_LOGS/heat_%a.slurm.out
#SBATCH -e PATH_TO_LOGS/heat_%a.slurm.err
#SBATCH -n N_TASKS
#SBATCH --gres=gpu:N_GPUS
#SBATCH --cpus-per-gpu=N_CPUS

export ligands_dat_file=LIGANDS_DATA_FILE
export transformations_file=TRANSFORMATIONS_DATA_FILE

source $MEZEHOME/parse.sh

export MEZEHOME=PATH_TO_MEZE

lig_i=$SLURM_ARRAY_TASK_ID
ligand=${ligand_array[$lig_i]}

start=`date +%s`

python $MEZEHOME/equilibrate.py "$ligand" PROTOCOLFILE

end=`date +%s`
runtime=$((end - start))

echo "Finished in $runtime seconds, or $((runtime/60)) minutes"
