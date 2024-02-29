#!/bin/bash

###TODO IN SETUP.PY:
# setup $MEZEHOME
###
###TODO How to activate correct environment? 

N_LIGANDS=NUMBER_OF_LIGANDS
ENGINE=ENGINE
AFE_INPUT_DIR=AFE_PATH

export ligands_dat_file=LIGANDS_DATA_FILE
export transformations_file=TRANSFORMATIONS_DATA_FILE

source $MEZEHOME/parse.sh 

meze_job_id=$(sbatch --parsable --array=0-$((${#transformations_array[@]}-1)) $AFE_INPUT_DIR/04_meze.sh)
echo "Preparing AFE with slurm job $prepafe_job_id"

for i in "${!transformations_array[@]}"
do
    afe_job_id=$(sbatch --dependency=afterok:${meze_job_id} --parsable --job-name=${transformations_array[i]} --array=0-$((${n_windows_array[i]}-1)) $AFE_INPUT_DIR/05_run_$ENGINE.sh ${transformations_array[i]} "${lambdas_array[i]}")
    echo "Submitted AFE slurm job $afe_job_id"
done
