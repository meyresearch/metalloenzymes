# Welcome to Flare Python Interpreter.
# Documentation can be found at Python > Documentation.
# This default python script can be edited at:
# '/Users/jasminguven/Library/Application Support/Cresset BMD/Flare/python/default-scripts/InterpreterDefaultScript.py'
from cresset import flare
import numpy as np
import os

project = flare.main_window().project

for i in range(len(project.ligands)):
    ligand = project.ligands[i]
    ligand_properties = ligand.properties
    ligand_name = ligand_properties["Title"]
    print(ligand_name)

    ligand.write_file(f"/home/jguven/projects/metalloenzymes/all_ligands/vim2/{ligand_name}.mol2", "mol2")
    ligand.write_file(f"/home/jguven/projects/metalloenzymes/all_ligands/vim2/{ligand_name}.pdb", "pdb")   
    ligand.write_file(f"/home/jguven/projects/metalloenzymes/all_ligands/vim2/{ligand_name}.sdf", "sdf")