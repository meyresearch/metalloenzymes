#! /bin/bash

### TO BE DONE IN SETUP.PY:
# N_LIGANDS=NUMBER_OF_LIGANDS
N_LIGANDS=16
AFE_INPUT_DIR="/home/jguven/projects/alchemistry/kpc2/partially_protonated_ligand/afe/"
# setup $MEZEHOME
###
INPUT_FILE=$1
LIGAND_CHARGE=$2

python $MEZEHOME/prepare.py --input-pdb-file $INPUT_FILE --ligand-charge $LIGAND_CHARGE

solvation_job_id=$(sbatch --parsable --array=0-$(($N_LIGANDS-1)) $AFE_INPUT_DIR/slurm_add_water.sh)
echo "Adding water with slurm job $solvation_job_id"

heating_job_id=$(sbatch --dependency=afterok:${solvation_job_id} --parsable --array=1-$N_LIGANDS $AFE_INPUT_DIR/slurm_heat_meze.sh)
echo "Heating meze with slurm job $heating_job_id"


# FEP prep 
# jidfep=$(sbatch --dependency=afterany:${jidlig} --parsable --array=0-$((${#trans_array[@]}-1)) run_fepprep_slurm.sh)
# echo "FEP prep jobid is $jidfep"


### run_fepprep_slurm.sh
#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --job-name=fepprep
#SBATCH -o ../slurm_logs/fepprep_%A_%a.out
#SBATCH -e ../slurm_logs/fepprep_%A_%a.err

# date
# echo "Folder for these runs is : $MAINDIRECTORY"
# echo "Network file is : $net_file"

# trans=${trans_array[$SLURM_ARRAY_TASK_ID]}

# python $scripts_dir/fepprep.py $trans

###

### fepprep.py

# network.afe_prep()

###

# production run & analysis