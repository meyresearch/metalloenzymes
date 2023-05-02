import functions
import argparse
import ast
import os

def clean_arguments(arguments):
    """
    Check arguments and clean them

    Parameters:
    -----------
    arguments: Namespace
        command-line arguments
    Return:
    -------
    cleaned Namespace of arguments
    
    """
    if arguments.ligand_path is None:
        arguments.ligand_path = os.getcwd()
    elif not os.path.isdir(arguments.ligand_path):
        raise argparse.ArgumentTypeError(f"{arguments.ligand_path} does not exist")
    return arguments


def main():

    parser = argparse.ArgumentParser(description="MEZE: MEtalloenZymE FF-builder for alchemistry")
   
    parser.add_argument("-i",
                        "--input-pdb",
                        dest="protein",
                        required=True, 
                        type=functions.file_exists,
                        help="input pdb file for the metalloenzyme/protein")
    
    parser.add_argument("-l",
                        "--ligand-path",
                        dest="ligand_path",
                        help="path to ligand files")
    
    arguments = parser.parse_args()
    # arguments = clean_arguments(arguments)
    ligand_path = functions.path_exists(arguments.ligand_path)

    functions.prepare_system(ligand_path, arguments.protein)

if __name__ == "__main__":
    main()