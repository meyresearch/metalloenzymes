import warnings
warnings.filterwarnings("ignore")
import BioSimSpace as bss 
import pandas as pd


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


def get_mol2_info(mol2_file: str) -> tuple:
    """
    Get 3 lists containing the molecule, atom and remaining info from a mol2 file
    """
    remaining_info = None
    with open(mol2_file, "r") as ifile:
        mol2_lines = ifile.readlines()
    get_index = lambda kword, lines: [i for i in range(len(lines)) if kword in lines[i]][0]
    atom_info_index = get_index("ATOM", mol2_lines) + 1
    molecule_info = mol2_lines[:atom_info_index]
    try:
        bond_info_index = get_index("BOND", mol2_lines)
        remaining_info = mol2_lines[bond_info_index:]
        atom_info = [element.split() for element in mol2_lines[atom_info_index:bond_info_index]]
        # atom_info = mol2_lines[atom_info_index:bond_info_index]
    except IndexError:
        print(f"No bond info found for file {mol2_file}")
        substructure_index = get_index("SUBSTRUCTURE", mol2_lines)
        atom_info = [element.split() for element in mol2_lines[atom_info_index:substructure_index]]
        # atom_info = mol2_lines[atom_info_index:substructure_index]
        remaining_info = mol2_lines[substructure_index:]
    return molecule_info, atom_info, remaining_info


ligand_directory = "../inputs/ligands/"
protein_directory = "../inputs/protein/"
xtal_filename = ligand_directory + "ligand_8_xtal"
parm_filename = protein_directory + "PT1.mol2"

ligand_1 = bss.IO.readMolecules(xtal_filename+".sdf")[0]
bss.IO.saveMolecules(xtal_filename,ligand_1, "mol2")
xtal_molecule, xtal_atom, xtal_remaining = get_mol2_info(xtal_filename+".mol2")
parm_molecule, parm_atom, parm_remaining = get_mol2_info(parm_filename)






x, y, z = 2, 3, 4

# for i in range(len(parm_atom)):
#     parm_atom[i][x] = xtal_atom[i][x]
#     parm_atom[i][y] = xtal_atom[i][y]
#     parm_atom[i][z] = xtal_atom[i][z]
    
# with open(ligand_directory+"test_8.mol2", "w") as ofile:
#     ofile.writelines(parm_molecule)
#     ofile.writelines(parm_atom)
#     ofile.writelines(parm_remaining)
