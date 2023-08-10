""" Useful functions """

import argparse
import os
import glob
import BioSimSpace as bss
import pathlib
import csv
import numpy as np
import Network 
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


def prepare(Protein, Network, AFE): 
    """
    Prepare ligands and protein for AFE calculations.

    Parameters:
    -----------
    Protein: Protein
        Protein class object
    Ligands: Ligands
        Ligands class object
    AFE: AlchemicalFreeEnergy
        AlchemicalFreeEnergy class object

    Return:
    -------
    """
    protein_water_complex_file = Protein.create_complex()
    tleap_command = Protein.write_tleap_input(protein_water_complex_file)
    Protein.tleap(tleap_command)

    afe_directory = AFE.create_directory()
    Network.create_ligand_dat_file(afe_directory)
    Network.create_network_files(AFE.engine)
    AFE.create_dat_file(Protein=Protein, Network=Network)


# def solvate(Protein, Ligands):

#     for i in range(Ligands.n_ligands()):
#         ligand_number = Ligands.names[i].split("_")[-1]
#         print(f"Solvating ligand {ligand_number}")
#         ligand_parameters = bss.Parameters.gaff2(Ligands.molecules[i], net_charge=Ligands.charge).getMolecule()


