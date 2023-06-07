import glob
import functions 
import shutil
import subprocess


path_to_ligands = "/home/jguven/projects/metalloenzymes/all_ligands/vim2/"
ligand_files = sorted(glob.glob(path_to_ligands+"*.pdb"))
ligand_names = [filename.split("/")[-1].replace(".pdb", "") for filename in ligand_files]

nb_directory = "/home/jguven/projects/metalloenzymes/nonbonded_model_vim2/"
prep_directory = nb_directory + "system_preparation/"
protein_file = prep_directory + "vim2.pdb"

for i in range(len(ligand_names)):
    ligand_directory = nb_directory + ligand_names[i] 
    functions.create_dirs(ligand_directory)

    for j in [1, 2]:
        step_directory = ligand_directory + f"/step_4n{j}/"
        wat_pdb_file = prep_directory + f"WAT.pdb"
        wat_mol2_file = prep_directory + f"WAT.mol2"
        functions.create_dirs(step_directory)
        shutil.copy(ligand_files[i], step_directory)
        shutil.copy(protein_file, step_directory)
        shutil.copy(wat_pdb_file, step_directory)
        shutil.copy(wat_mol2_file, step_directory)
        shutil.copy(nb_directory+f"step_4n{j}/vim2.in", step_directory)

prepare_ligands = subprocess.call("./prepare_ligands.sh")
prepare_protein = subprocess.call("./prepare_protein.sh")
run_mcbp_step_4 = subprocess.call("./nonbonded_mcpb.sh")



