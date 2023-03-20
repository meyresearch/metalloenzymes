# Welcome to Flare Python Interpreter.
# Documentation can be found at Python > Documentation.
# This default python script can be edited at:
# '/Users/jasminguven/Library/Application Support/Cresset BMD/Flare/python/default-scripts/InterpreterDefaultScript.py'
from cresset import flare
import numpy as np
import os

project = flare.main_window().project

ligand_numbers = np.arange(1,17,1)

for i in range(len(project.ligands)):
    ligand = project.ligands[i]
    ligand_properties = ligand.properties
    filename = str(ligand_properties["Title"])
    print(filename)
    home_dir = os.path.expanduser("~")
    ligand.write_file(f"{home_dir}/projects/metalloenzymes/all_ligands/kpc2/{filename}.sdf", "sdf")
