import BioSimSpace as bss
import argparse
import functions
from argparse import RawTextHelpFormatter


def create_network(ligands, ligand_names, ligand_path):
    """
    
    """
    transformations, lomap_scores = bss.Align.generateNetwork(ligands, plot_network=True, names=ligand_names, work_dir=ligand_path)
    network_dict = {}
    named_transformations = [(ligand_names[transformation[0]], ligand_names[transformation[1]]) for transformation in transformations]
    
    for transformation, score in zip(named_transformations, lomap_scores):
        network_dict[transformation] = score

    with open(ligand_path+f"/meze_lomap_network.csv", "w") as lomap_out:
        for key, value in network_dict.items():
            lomap_out.write(f"{key}: {value}\n")
    
    return network_dict


def main():

    parser = argparse.ArgumentParser(description="MEZE: MEtalloenZymE FF-builder for alchemistry\nCreate and edit a LOMAP network.",
                                     formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("-l",
                        "--ligand-path",
                        dest="ligand_path",
                        help="path to ligand files")
    arguments = parser.parse_args()

    ligand_path = functions.path_exists(arguments.ligand_path)
    ligand_files = functions.get_ligand_files(ligand_path)
    ligands = [bss.IO.readMolecules(file)[0] for file in ligand_files]
    ligand_names = [functions.get_filenames(filepath) for filepath in ligand_files]

    network_dict = create_network(ligands, ligand_names, ligand_path)



if __name__ == "__main__":
    main()