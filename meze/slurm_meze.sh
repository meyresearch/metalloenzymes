#! /bin/bash

### TO BE DONE IN SETUP.PY:

AFE_INPUT_DIR="/home/jguven/projects/alchemistry/kpc2/partially_protonated_ligand/afe/"
# setup $MEZEHOME
###

### Needs to be figured out:
N_LIGANDS=16
ENGINE="SOMD"
# AND HOW INPUT FILE NAMES FOR LIGANDS.DAT AND MEZE_NETWORK.DAT ARE TAKEN

INPUT_FILE=$1
LIGAND_CHARGE=$2

python $MEZEHOME/prepare.py --input-pdb-file $INPUT_FILE --ligand-charge $LIGAND_CHARGE

export ligands_dat_file=$3
export transformations_file=$4

source $MEZEHOME/parse.sh 

solvation_job_id=$(sbatch --parsable --array=0-$(($N_LIGANDS-1)) $AFE_INPUT_DIR/slurm_add_water.sh)
echo "Adding water with slurm job $solvation_job_id"

heating_job_id=$(sbatch --dependency=afterok:${solvation_job_id} --parsable --array=1-$N_LIGANDS $AFE_INPUT_DIR/slurm_heat_meze.sh)
echo "Heating meze with slurm job $heating_job_id"

prepafe_job_id=$(sbatch --dependency=afterok:${heating_job_id} --parsable $AFE_INPUT_DIR/slurm_prepafe.sh)
echo "Preparing AFE with slurm job $prepafe_job_id"

for i in "${!transformations_array[@]}"
do
    afe_job_id=$(sbatch --parsable --job-name=${transformations_array[i]} --array=0-$((${n_windows_array[i]}-1)) $AFE_INPUT_DIR/run_$ENGINE.sh ${transformations_array[i]} "${lambdas_array[i]}")
    echo "Submitted AFE slurm job $afe_job_id"
done