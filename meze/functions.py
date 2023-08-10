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


def write_tleap_setup_file(group_name, protein_water_complex, working_directory, forcefield="ff14SB", water_model="tip3p"):
    """
    Write tleap input file for setting up protein with given parameters

    Parameters:
    -----------
    group_name: str
        group name for system, e.g. vim2/kpc2/ndm1/etc
    protein_water_complex: str
        pdb file with water and protein combined
    working_directory: str
        project working directory
    forcefield: str
        protein force field, default is ff14sb
    water_model: str
        name of water model, default is tip3p

    Return:
    -------
    tleap_command: str
        tleap run command
    """
    tleap_input_file = working_directory + "/tleap.in"
    tleap_output_file = working_directory + "/tleap.out"
    save_file = working_directory + "/" + group_name

    with open(tleap_input_file, "w") as tleap_in:
        tleap_in.write(f"source leaprc.protein.{forcefield}\n")
        tleap_in.write(f"source leaprc.water.{water_model}\n")
        tleap_in.write(f"complex = loadpdb {protein_water_complex}\n")
        tleap_in.write(f"saveamberparm complex {save_file}_tleap.prm7 {save_file}_tleap.rst7\n")
        tleap_in.write("quit")

    return f"tleap -s -f {tleap_input_file} > {tleap_output_file}"


def run_tleap(command):
    """
    Run tleap 

    Parameters:
    -----------
    command: str
        tleap command for the current system

    Return:
    -------
    """
    os.system(command)


def add_xtal_waters(protein_file, water_file, working_directory):
    """
    Combine protein file with crystallographic waters

    Parameters:
    -----------
    protein: str
        protein pdb file
    water: str
        xtal water pdb file

    Return:
    -------
    pdb_file: str
        pdb file with water and protein combined
    """
    protein = bss.IO.readMolecules(protein_file)
    water = bss.IO.readMolecules(water_file)
    output = working_directory + "/protein_water_complex"
    complex = protein + water
    bss.IO.saveMolecules(output, complex, fileformat="pdb")
    return output + ".pdb"


def get_water_file(working_directory):
    """
    Get water.pdb file from current working directory
    #TODO maybe get it as input 

    Parameters:
    -----------
    Return:
    -------
    water_pdb: str
        crystal waters in pdb file
    """
    return working_directory + "/water.pdb"


def create_afe_directory(working_directory):
    """
    Create afe directory in the project working directory
    This is the directory where network, protocol and ligand datafiles are saved

    Parameters:
    -----------
    working_directory: str
        project working directory

    Return:
    -------
    afe_directory: str
    """
    try:
        afe_directory = working_directory + "/afe/"
        pathlib.Path(afe_directory).mkdir(parents=False, exist_ok=False)
    except FileNotFoundError as e:
        print(f"Could not create afe directory. Pathlib raised error: {e}")
    except FileExistsError as e:
        print(f"Could not create afe directory. Pathlib raised error: {e}")
    return afe_directory


def create_ligands_file(working_directory, names):
    """
    Create ligands.dat file in afe directory from ligand names

    Parameters:
    -----------
    working_directory: str 
        afe directory
    names: list
        ligand names

    Return:
    -------
    ligands_dat: str
        ligands datafile
    """
    ligands_dat = working_directory + "ligands.dat"
    with open(ligands_dat, "w") as ligands_file:
        writer = csv.writer(ligands_file)
        for ligand in names:
            writer.writerow([ligand])
    return ligands_dat


def set_n_windows(lomap_score, threshold=0.4, n_normal=11, n_difficult=17):
    """
    Set the number of lambda windows for given transformation depending on lomap score

    Parameters:
    -----------
    lomap_score: float
        lomap score for transformation
    threshold: float
        threshold value of lomap score to define difficult transformations
    n_normal: int
        number of lambda windows for a "normal" transformation
    n_difficult: int
        number of lambda windows for difficult transfromations

    Return:
    -------
    n_windows: int
        number of lambda windows for given 
    """
    if lomap_score == None or lomap_score < float(threshold):
        n_windows = n_difficult
    else:
        n_windows = n_normal
    return n_windows


def create_lambda_list_bash(n_windows):
    """
    Create a bash-readable list of evenly spaced lambda values between 0 and n_windows

    Parameters:
    -----------
    n_windows: 
        number of lambda windows

    Return:
    -------
    bash_list: str
        a bash-readable list of lambda values
    """
    lambda_list_numpy = list(np.linspace(0, 1, int(n_windows)))
    lambda_list = [format(item, ".4f") for item in lambda_list_numpy]
    return " ".join(lambda_list)



def create_network_files(working_directory, network_dict, engine):
    """
    _summary_

    Parameters:
    -----------
    working_directory: str 
        project working directory
    network: dict
        lomap network
    engine: str
        MD engine
    
    Return:
    -------
    networks: tuple
        forward transformation network, backward transformation network
    """
    forward = working_directory + "network_fwd.dat"
    backward = working_directory + "network_bwd.dat"

    with open(forward, "w") as network_file:
        for transformation, lomap_score in network_dict.items():
            n_windows = set_n_windows(lomap_score)
            lambda_array_bash = create_lambda_list_bash(n_windows)
            network_file.write(f"{transformation[0]}, {transformation[1]}, {n_windows}, {lambda_array_bash}, {engine}\n")
    with open(backward, "w") as network_file:
        for transformation, lomap_score in network_dict.items():
            n_windows = set_n_windows(lomap_score)
            lambda_array_bash = create_lambda_list_bash(n_windows)
            network_file.write(f"{transformation[1]}, {transformation[0]}, {n_windows}, {lambda_array_bash}, {engine}\n")
    
    return forward, backward


def create_protocol_dat(working_directory, ligand_ff, forcefield, water_model, engine, box_edges, box_shape, sampling_time):
    """
    Create a protocol.dat file for afe runs

    Parameters:
    -----------
    ligand_ff: str
        ligand force field, default gaff2
    forcefield: str
        protein force field, default ff14sb
    water_model: str
        water model, default tip3p
    engine: str
        default SOMD
    box_edges: str
        box edges in angstroms, default 20 A
    box_shape: str
        default orthorhombic
    
    Return:
    -------
    protocol_dat: str
        protocol.dat file
    """
    protocol = [f"ligand forcefield = {ligand_ff}", 
            f"protein forcefield = {forcefield}", 
            f"solvent = {water_model}", 
            f"box edges = {box_edges}*angstrom", 
            f"box shape = {box_shape}", 
            f"protocol = default",
            f"sampling = {sampling_time}*ns",
            f"engine = {engine}"]
    protocol_file = working_directory + "/protocol.dat"

    with open(protocol_file, "w") as protocol_file:
        writer = csv.writer(protocol_file)
        for protocol_line in protocol:
            writer.writerow([protocol_line])
    return protocol_file


def prepare_system(group_name, 
                   working_directory, 
                   protein, 
                   forcefield, 
                   water_model, 
                   engine, 
                   ligand_ff,
                   box_edges,
                   box_shape,
                   sampling_time,
                   ): 
    """
    Prepare ligands and protein for AFE calculations
    #TODO should probably separate to ligand and protein

    Parameters:
    -----------
    group_name: str
        group name for system, e.g. vim2/kpc2/ndm1/etc
    working_directory: str
        project working directory
    protein: 
        input pdb file for metalloenzyme/protein
    forcefield: str
        protein force field
    water_model: 
        water model

    Return:
    -------
    """
    ligand_path = path_exists(working_directory + "/inputs/ligands/") 
    ligand_files = get_ligand_files(ligand_path)
    ligands = [bss.IO.readMolecules(file)[0] for file in ligand_files]
    ligand_names = [get_filenames(filepath) for filepath in ligand_files]

    network_dict = network.create_network(ligands, ligand_names, ligand_path)
    # network_dict = network.edit_network(ligand_path, network_dict)
    
    # protein prep
    protein_path = path_exists(working_directory + "/inputs/protein/")
    water_file = get_water_file(protein_path)
    protein_water_complex = add_xtal_waters(protein, water_file, protein_path)
    tleap_command = write_tleap_setup_file(group_name, protein_water_complex, protein_path, forcefield, water_model)
    run_tleap(tleap_command)

    # create afe production directory and add .dat files to it
    afe_directory = create_afe_directory(working_directory)
    ligands_dat_file = create_ligands_file(afe_directory, ligand_names)
    forward_network_dat_file, backward_network_dat_file = create_network_files(afe_directory, network_dict, engine)
    protocol_dat_file = create_protocol_dat(afe_directory, ligand_ff, forcefield, water_model, engine, box_edges, box_shape, sampling_time)



