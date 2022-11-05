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
    ligand.write_file(f"/Users/jasminguven/projects/metalloenzymes/kpc2/inputs/ligands/docked_{ligand_numbers[i]}.mol2",
                       "mol2")
