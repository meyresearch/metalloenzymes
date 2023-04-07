""" Useful functions """

import argparse
import os
import glob
import BioSimSpace as bss


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


# def path_exists(path):
#     """ 
#     Check that given path exists
    
#     Parameters:
#     -----------
#     path: str
#         full path 
#     Return:
#     -------
#     path: str
#         existing path
#     """
#     if path is None:
#         path = os.getcwd()
#     elif not os.path.isdir(path):
#         raise argparse.ArgumentTypeError(f"{path} does not exist")
#     return path


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


def get_names(path):
    """
    
    """
    file_name = path.split("/")[-1]
    extension = file_name.split(".")[-1]
    print(extension)
    return file_name.replace("." + extension, "")


def prepare_system(ligand_path, protein):
    """
    
    """
    print(ligand_path, protein)
    # Adapted from dmbol.py:
    ligand_files = read_files(f"{ligand_path}/*.sdf")
    ligand_files += read_files(f"{ligand_path}/*.mol2")
    if len(ligand_files) < 2:
        raise IOError(f"Path {ligand_path} must contain at least two sdf or mol2 files.")
    ligands = [bss.IO.readMolecules(file)[0] for file in ligand_files]

    ligand_names = [get_names(filepath) for filepath in ligand_files]
    print(ligand_names)


    # read ligand files 
    # ligand_files = sorted(glob.glob(f"{ligand_path}*.sdf"))
    # ligands = [bss.IO.readMolecules(filepath)[0] for filepath in ligand_files]
    # ligand_names = [filepath.split("/")[-1].replace(".sdf","") for filepath in ligand_files]
    # transformations, lomap_scores = bss.Align.generateNetwork(ligands, plot_network=True, names=ligand_names, work_dir=ligand_path)

    # perturbation_network_dict = {}
    # transformations_named = [(ligand_names[transf[0]], ligand_names[transf[1]]) for transf in transformations]
    # for transformation, score in zip(transformations_named, lomap_scores):
    #     perturbation_network_dict[transformation] = score

    # with open(ligand_path+f"/lomap_{system}.csv", "w") as lomap_out:
    #     for key, value in perturbation_network_dict.items():
    #         lomap_out.write(f"{key}: {value}\n")

