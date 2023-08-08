""" Useful functions """

import argparse
import os
import glob
import BioSimSpace as bss
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


def get_ligand_files(ligand_path):
    """
    Read in all .sdf and .mol2 files in the provided path 
    Parameters:
    -----------
    ligand_path: str
        path to directory containing mol2 and/or sdf files
    Return:
    -------
    list
        list of ligand filenames
    """
    # Adapted from dbmol.py:
    ligand_files = read_files(f"{ligand_path}/*.sdf")
    ligand_files += read_files(f"{ligand_path}/*.mol2")
    
    if len(ligand_files) < 2:
        raise IOError(f"Path {ligand_path} must contain at least two sdf or mol2 files.")
    
    return ligand_files


def prepare_system(ligand_path, protein): # should probably separate to ligand and protein
    """
    
    """
    ligand_files = get_ligand_files(ligand_path)
    ligands = [bss.IO.readMolecules(file)[0] for file in ligand_files]
    ligand_names = [get_filenames(filepath) for filepath in ligand_files]

    network_dict = network.create_network(ligands, ligand_names, ligand_path)
    network = network.edit_network(network_dict)
    
    # create afe production directory and add .dat files to it

    # create protocol and add to above directory



