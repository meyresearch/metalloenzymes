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
    return os.path.abspath(file)


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
    return path + "/"


def check_int(input):
    """
    Check if given input is an integer

    Parameters:
    -----------
    input: any
        user input

    Return:
    -------
    value: int
        user input converted to int
    """
    try:
        value = int(input)
    except ValueError:
        print("Error: number of minimisation steps and/or runtime should be integer.")
    return value


def check_float(input):
    """
    Check if given input is a float

    Parameters:
    -----------
    input: any
        user input

    Return:
    -------
    value: float
        user input converted to float
    """
    try:
        value = float(input)
    except ValueError:
        print("Error: number of minimisation steps and/or runtime should be a number.")
    return value


def check_positive(number):
    """
    Check that given number is positive.

    Parameters:
    -----------
    number: int/float
        user inputted number

    Return:
    -------
    value: int/float
        positive integer value
    """
    try:
        if number <= 0:
            raise argparse.ArgumentTypeError(f"{number} is an invalid positive integer value")
    except argparse.ArgumentTypeError as message:
        print(message)
    return number


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


def mkdir(directory):
    """
    Create directory

    Parameters:
    -----------
    directory: str
        full path to directory to be created

    Return:
    -------
    directory: str
        newly created directory
    """
    path_exists = os.path.exists(directory)
    if not path_exists:
        os.makedirs(directory)
    return directory


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


def convert_to_units(value, units):
    """
    Convert given value to given units

    Parameters:
    -----------
    value: any
        measure of time, temperature, length, etc.
    units: constant
        one of bss.Units defined in definitions.py

    Return:
    -------
    One of bss.Units
    """
    return value * units


def input_to_dict(file):
    """
    Convert an input file to dictionary

    Parameters:
    -----------
    file: str
        full path to input file

    Return:
    -------
    dict:
        input options as a dictionart
    """
    with open(file, "r") as file:
        lines = file.readlines()
    clean_lines = [line.strip().split("=") for line in lines if "=" in line]
    
    dictionary = {}
    for _key, _value in clean_lines:
        key = _key.strip()
        value = _value.strip()
        if value.isdigit():
            value = int(value)
            dictionary[key] = value
        else:
            try:
                value = float(value)
                dictionary[key] = value
            except ValueError:
                dictionary[key] = value
    return dictionary


def write_slurm_script(template_file, path, log_dir, protocol_file, extra_options=None, extra_lines=None):
    """
    Write a slurm script for running different stages of MEZE
    Parameters:
    -----------
    template_file: str
        name of the template run script
    path: str
        path where slurm script will be saved
    log_dir: str
        directory for outputting slurm logs
    protocol_file: str
        full path to the protocol file
    extra_options: None or dict
        a dictionary of options to add to the slurm script; "REPLACE": value. Overrides defaults.
    extra_lines: None or dict
        a dictionary of options to append to the default options

    Return:
    -------
    file: str
        slurm script 
    """
    output = path + "/" + template_file
    meze = os.environ["MEZEHOME"]
    template = meze + "/" + template_file
    with open(template, "r") as file:
        lines = file.readlines()
    
    if extra_options:
        options = extra_options
    elif not extra_options:
        options = {"PATH_TO_LOGS": log_dir,
                   "N_TASKS": str(1),
                   "N_GPUS": str(1), 
                   "N_CPUS": str(10),
                   "MEMORY": str(4069),
                   "PATH_TO_MEZE": meze,
                   "PROTOCOLFILE": protocol_file}    
    if extra_lines:
        for key, value in extra_lines.items():
            options[key] = str(value)

    with open(output, "w") as file:
        for line in lines:
            for key, value in options.items():
                line = line.replace(key, value)
            file.write(line)
    os.system(f"chmod +x {output}")
    return output


def separate(string, sep="~"):
    """
    Split given string by given separator

    Parameters:
    -----------
    string: str

    Return:
    -------
    list: 
        separated string as a list
    """
    return string.split(sep)
