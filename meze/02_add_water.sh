#!/bin/bash

#SBATCH -o PATH_TO_LOGS/add_water_%a.slurm.out
#SBATCH -e PATH_TO_LOGS/add_water_%a.slurm.err



export MEZEHOME=PATH_TO_MEZE

LIG_NUMBER=$SLURM_ARRAY_TASK_ID

start=`date +%s`

python $MEZEHOME/solvate.py "ligand_$LIG_NUMBER" PROTOCOLFILE

end=`date +%s`
runtime=$((end - start))

echo "Finished in $runtime seconds, or $((runtime/60)) minutes"
