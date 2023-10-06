#!/bin/bash

###TODO IN SETUP.PY:
# setup $MEZEHOME
###

N_LIGANDS=NUMBER_OF_LIGANDS
QMMM_INPUT_DIR=AFE_PATH

# 01_setup_quantum.sh

solvation_job_id=$(sbatch --parsable --array=1-$N_LIGANDS $QMMM_INPUT_DIR/02_add_water.sh)
echo "Adding water with slurm job $solvation_job_id"

# 03_prepare.sh


# minimise 
# heat

# production qm/mm
    

heating_job_id=$(sbatch --dependency=afterok:${solvation_job_id} --parsable --array=1-$N_LIGANDS $QMMM_INPUT_DIR/03_heat_meze.sh)
echo "Heating meze with slurm job $heating_job_id"

meze_job_id=$(sbatch --dependency=afterok:${heating_job_id} --parsable --array=1-$N_LIGANDS $QMMM_INPUT_DIR/04_meze.sh)
echo "Preparing AFE with slurm job $prepafe_job_id"

qmmm_job_id=$(sbatch --dependency=afterok:${meze_job_id} --parsable --job-name=1-$N_LIGANDS --array=0-$((${n_windows_array[i]}-1)) $QMMM_INPUT_DIR/run_QMMM.sh ${transformations_array[i]} "${lambdas_array[i]}")
echo "Submitted AFE slurm job $qmmm_job_id"
