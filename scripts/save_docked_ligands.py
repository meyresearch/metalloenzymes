# Welcome to Flare Python Interpreter.
# Documentation can be found at Python > Documentation.
# This default python script can be edited at:
# '/Users/jasminguven/Library/Application Support/Cresset BMD/Flare/python/default-scripts/InterpreterDefaultScript.py'
from cresset import flare
import numpy as np
import os

project = flare.main_window().project

title = str(project.proteins[0].title)

protein_name = ""
if "6DD0" in title:
    protein_name = "vim_2"
elif "6D14" in title:
    protein_name = "kpc_2"

print(protein_name)
ligand_numbers = np.arange(1,17,1)

for i in range(len(project.ligands)):
    ligand = project.ligands[i]

    ligand_properties = ligand.properties
    ligand_title = ligand_properties["Title"]

    if "_D" in str(ligand_title):
        filename = str(ligand_properties["Filename"])
        ligand_name_no_extension = filename.split(".")[0]
        ligand_number = ligand_name_no_extension.split("_")[-1]
        print(ligand_number)
        if ligand_number != "":
            ligand.write_file(f"/Users/jasminguven/projects/metalloenzymes/dft/ligand_{ligand_number}/ligand_{ligand_number}.pdb", "pdb")
    elif str(ligand_title) == "A P9T 403":
        ligand.write_file(f"/Users/jasminguven/projects/metalloenzymes/dft/ligand_8/ligand_8.pdb", "pdb")
                