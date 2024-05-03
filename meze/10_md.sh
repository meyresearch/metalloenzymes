#!/bin/bash
#SBATCH -o PATH_TO_LOGS/md_%a.slurm.out
#SBATCH -e PATH_TO_LOGS/md_%a.slurm.err
#SBATCH -n N_TASKS
#SBATCH --gres=gpu:N_GPUS
#SBATCH --cpus-per-gpu=N_CPUS

export ligands_dat_file=LIGANDS_DATA_FILE
export transformations_file=TRANSFORMATIONS_DATA_FILE
export MEZEHOME=PATH_TO_MEZE

source $MEZEHOME/parse.sh

lig_i=$SLURM_ARRAY_TASK_ID
ligand=${ligand_array[$lig_i]}

start=`date +%s`

python $MEZEHOME/md.py $ligand --input-pdb-file INPUTPROTEINFILE --ligand-charge -2 --group-name vim2

end=`date +%s`
runtime=$((end - start))

echo "Finished in $runtime seconds, or $((runtime/60)) minutes"
