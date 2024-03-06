#!/bin/bash

#SBATCH -o PATH_TO_LOGS/extra_meze_%a.slurm.out
#SBATCH -e PATH_TO_LOGS/extra_meze_%a.slurm.err

extra_transformations=EXTRA_TRANSFORMATIONS

export MEZEHOME=PATH_TO_MEZE

start=`date +%s`

id=$SLURM_ARRAY_TASK_ID

source $MEZEHOME/parse_extra_edges.sh
transformation=${transformations_array[$id]}

python $MEZEHOME/meze.py PROTOCOLFILE $transformation --extra-transformations $extra_transformations

end=`date +%s`
runtime=$((end - start))

echo "Finished in $runtime seconds, or $((runtime/60)) minutes"
