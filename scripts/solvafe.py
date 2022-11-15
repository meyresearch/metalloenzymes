import warnings
warnings.filterwarnings("ignore")
import BioSimSpace as bss 
import pandas as pd
import MDAnalysis as mda


ligand_directory = "../inputs/ligands/"
protein_directory = "../inputs/protein/"
xtal_filename = ligand_directory + "ligand_8_xtal"
parm_filename = protein_directory + "PT1.mol2"


parm_u = mda.Universe(protein_directory+"PT1.mol2")
xtal_u = mda.Universe(ligand_directory+"ligand_8_xtal.mol2")
new_coords = xtal_u.atoms.positions
new_universe = parm_u.load_new(ligand_directory+"ligand_8_xtal.mol2", format="MOL2")

new_universe.select_atoms("all").write(ligand_directory+"fixed_ligand_8.mol2")


ligand_1 = bss.IO.readMolecules(ligand_directory+"fixed_ligand_8.mol2")[0]

