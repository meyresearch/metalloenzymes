#!/bin/bash

###TODO IN SETUP.PY:
# setup $MEZEHOME
###
###TODO How to activate correct environment? 


ENGINE=ENGINE
AFE_INPUT_DIR=AFE_PATH

export transformations_file=TRANSFORMATIONS_DATA_FILE

dos2unix "$transformations_file"

transformations_array=()
lambdas_array=()
n_windows_array=()

while IFS="," read -a line
do
    transformation=${line[1]}~${line[3]}
    lambda=${line[7]}
    n_windows=${line[6]}
    transformations_array+=("$transformation")
    lambdas_array+=("$lambda")
    n_windows_array+=("$n_windows")
done < <(tail -n +2 "$transformations_file")

meze_job_id=$(sbatch --parsable --array=0-$((${#transformations_array[@]}-1)) $AFE_INPUT_DIR/04_meze.sh)
echo "Preparing AFE with slurm job $prepafe_job_id"

for i in "${!transformations_array[@]}"
do
    afe_job_id=$(sbatch --dependency=afterok:${meze_job_id} --parsable --job-name=${transformations_array[i]} --array=0-$((${n_windows_array[i]}-1)) $AFE_INPUT_DIR/05_run_$ENGINE.sh ${transformations_array[i]} "${lambdas_array[i]}")
    echo "Submitted AFE slurm job $afe_job_id"
done
