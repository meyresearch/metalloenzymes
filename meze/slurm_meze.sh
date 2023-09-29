#! /bin/bash

INPUT_FILE=$1
echo $INPUT_FILE

python $MEZEHOME/prepare.py --input-pdb-file $INPUT_FILE

# export MAINDIRECTORY <- working_directory
# export scripts_dir <- set to $MEZEHOME in setup.py
# export protein_file <- same as my --input-file but I do it with PDB file and as a command line argument
# export ligands_folder <- -lwd / default value os.getcwd() + /inputs/ligands/

# run network.prepare_network() here: -> it outputs some files that can be used in later python scripts

# ligand prep
# jidlig=$(sbatch --parsable --array=1-$N_LIGANDS run_ligprep_slurm.sh)

### run_ligprep_slurm.sh

# #!/bin/bash
# #SBATCH -n 1
# #SBATCH --gres=gpu:1
# #SBATCH --job-name=ligprep
# #SBATCH -o ../slurm_logs/ligprep_%A_%a.out
# #SBATCH -e ../slurm_logs/ligprep_%A_%a.err

# # sourcing
# source $BSS
# source $amber
# source $gromacs
# source extract_execution_model_bash.sh

# date
# echo "Folder for these runs is : $MAINDIRECTORY"
# echo "Ligands file is : $lig_file"

# lig=${lig_array[$SLURM_ARRAY_TASK_ID]}

# echo "prep for $lig..."
# python $scripts_dir/ligprep.py $lig

###

### ligprep.py

# take as user input file?
# minim_steps = 250
# runtime_short_nvt = 5 # ps
# runtime_nvt = 50 # ps 
# runtime_npt = 200 # ps

# main_dir = os.environ["MAINDIRECTORY"]

# amber_home = os.environ["AMBERHOME"]
# pmemd_path = amber_home + "/bin/pmemd.cuda" 

# protein_file = os.environ["protein_file"]
# ligands_folder = os.environ["ligands_folder"]

# ### Open 'ligands.dat', find sys.argv[1] 'th entry
# print (f"{sys.argv[0]} {sys.argv[1]}")
# lig_name = sys.argv[1]

# solvation and equilibration here
# e.g. could do solvation.py & equilibration.py 
# saves solvated & equilibrated files in inputs/ligands, inputs/protein/ and equilibration/ <- can be used as input for next step
###

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