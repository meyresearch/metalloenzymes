""" Useful functions """

import argparse
import os
import glob
import BioSimSpace as bss
import pathlib
import csv
import numpy as np
import network
from definitions import ROOT_DIRECTORY


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


def path_exists(path):
    """ 
    Check that given path exists

    Parameters:
    -----------
    path: str
        full path 

    Return:
    -------
    path: str
        existing path
    """
    if path is None:
        path = os.getcwd()
    elif not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"{path} does not exist")
    return path


def get_absolute_path(path):
    """
    Take a path, check it exists and convert it to absolute 

    Parameters:
    -----------
    path: str
        absolute/relative path 

    Return:
    -------
    str
        absolute path
    """
    correct_path = path_exists(path)
    clean_path = correct_path.replace("..", "")
    return ROOT_DIRECTORY + clean_path


def read_files(path):
    """
    Read and sort files in a given path

    Parameters:
    -----------
    path: str
        full path to directory containing files

    Return:
    -------
    list
        sorted list of files in path  
    """
    return sorted(glob.glob(path))


def get_filenames(path):
    """
    Read path to file and remove extension.

    Parameters:
    -----------
    path: str
        full path to directory containing files

    Return:
    -------
    str
        file name without extension
    """
    file_name = path.split("/")[-1]
    extension = file_name.split(".")[-1]
    return file_name.replace("." + extension, "")


def prepare_system(protein, ligands, afe): 
    """
    Prepare ligands and protein for AFE calculations.

    Parameters:
    -----------
    working_directory: str
        project working directory
    protein: 
        input pdb file for metalloenzyme/protein

    Return:
    -------
    """
    # ligand prep
    ligand_molecules = ligands.get_molecules()
    ligand_names = ligands.get_names()
    network_dict = network.create_network(ligand_molecules, ligand_names, ligands.path)

    # protein prep
    protein_water_complex = protein.create_complex()
    tleap_command = protein.write_tleap_input(protein_water_complex)
    protein.tleap(tleap_command)

    # create afe production directory and add .dat files to it
    afe_directory = afe.create_directory()
    ligands_dat_file = ligands.create_dat_file(afe_directory)
    forward_network_dat_file, backward_network_dat_file = network.create_network_files(afe_directory, network_dict, afe.engine)
    protocol_dat_file = afe.create_dat_file(Protein=protein, Ligands=ligands)



