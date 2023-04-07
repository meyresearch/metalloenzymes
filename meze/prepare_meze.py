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

    arguments = parser.parse_args()


if __name__ == "__main__":
    main()