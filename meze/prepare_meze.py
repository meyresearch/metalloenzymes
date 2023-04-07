import functions
import argparse


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
                        required=True, 
                        type=functions.file_exists,
                        help="path to ligand files")
    
    arguments = parser.parse_args()


if __name__ == "__main__":
    main()