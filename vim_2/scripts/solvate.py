import warnings
warnings.filterwarnings("ignore")
import BioSimSpace as bss 
import pandas as pd
import MDAnalysis as mda


def get_parmed_charge(parmed_file: str) -> float:
    """
    Get net charge of molecule from PARMED output file
    """
    with open(parmed_file, "r") as file:
        lines = file.readlines()
    mask_line = [i for i in range(len(lines)) if "mask" in lines[i]][0]
    info = [line.split() for line in lines[mask_line:] if "mask" not in line][1:-2]
    column_headers = [element for element in info[0] if element != "GB" and element != "LJ"]
    df = pd.DataFrame(info[1:], columns=column_headers)
    return df["Charge"].astype(float).sum()

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

# add hydroxide ions to the ligand to make the charge integer and then parameterise using gaff2!!!