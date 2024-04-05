#!/bin/bash

#SBATCH -o PATH_TO_LOGS/add_water_%a.slurm.out
#SBATCH -e PATH_TO_LOGS/add_water_%a.slurm.err

export ligands_dat_file=LIGANDS_DATA_FILE
export transformations_file=TRANSFORMATIONS_DATA_FILE

source $MEZEHOME/parse.sh 

export MEZEHOME=PATH_TO_MEZE

lig_i=$SLURM_ARRAY_TASK_ID
ligand=${ligand_array[$lig_i]}

start=`date +%s`

python $MEZEHOME/solvate.py "$ligand" PROTOCOLFILE

end=`date +%s`
runtime=$((end - start))

echo "Finished in $runtime seconds, or $((runtime/60)) minutes"
