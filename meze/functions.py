""" Useful functions """

import argparse
import os


def file_exists(file):
    """ 
    Check that given file exists
    
    Parameters:
    -----------
    file: str
        full path to file
    Return:
    -------
    file: str
        existing filepath
    """
    if not os.path.isfile(file):
        raise argparse.ArgumentTypeError(f"{file} does not exist")
    return file


def parse_cli_arguments():
    """
    Function to handle all argparse CLI arguments
    Return: 
    -------
    parser.parse_args(): Namespace
        parsed CLI arguments
    """
    parser = argparse.ArgumentParser(description="MEtalloenZymE FF-builder for alchemistry")


    parser.add_argument("-i",
                        "--input-pdb",
                        dest="protein",
                        required=True, 
                        type=file_exists,
                        help="input pdb file for the metalloenzyme/protein")

    return parser.parse_args()

