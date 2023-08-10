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
    network_dictionary = Network.create_dictionary()
    #TODO Edit dictionary?
    protein_water_complex_file = Protein.create_complex()
    prepared_protein_file = Protein.tleap(protein_water_complex_file)

    afe_directory = AFE.create_directory()
    ligands_datfile = Network.create_ligand_dat_file(afe_directory)
    network_forward, network_backward = Network.create_network_files(AFE.engine)
    protocol_file = AFE.create_dat_file(Protein=Protein, Network=Network)


def solvate(Protein, Network, AFE):
    """
    Solvate unbound and bound systems.

    Parameters:
    -----------
    Protein: 
        Protein class object
    Network: 
        Network class object
    AFE: 
        AlchemicalFreeEnergy class object

    Return:
    -------
    """
    ligands = Network.ligands
    for i in range(Network.n_ligands):
        ligand_number = Network.names[i].split("_")[-1]
        print(f"Solvating unbound ligand {ligand_number}")
        ligand_parameters = ligands[i].parameterise(Network.forcefield, Network.charge)
        unbound_box, unbound_box_angles = AFE.create_box(ligand_parameters)
        solvated_ligand = bss.Solvent.solvate(model=Protein.water_model, 
                                              molecule=ligand_parameters, 
                                              box=unbound_box,
                                              angles=unbound_box_angles)
        print(f"Solvating bound ligand {ligand_number}")        
        system_parameters = ligand_parameters + Protein.get_prepared_protein()
        bound_box, bound_box_angles = AFE.create_box(system_parameters)
        solvated_system = bss.Solvent.solvate(model=Protein.water_model,
                                              molecule=system_parameters,
                                              box=bound_box,
                                              angles=bound_box_angles)
        ligand_savename = Network.path + "ligand_" + ligand_number + "_solvated"
        system_savename = Protein.path + "system_" + ligand_number + "_solvated"
        bss.IO.saveMolecules(ligand_savename, solvated_ligand, ["PRM7", "RST7"])
        bss.IO.saveMolecules(system_savename, solvated_system, ["PRM7", "RST7"])


def equilibrate():
    pass
