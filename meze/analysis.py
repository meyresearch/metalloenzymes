import BioSimSpace as bss
import glob
import numpy as np
import pandas as pd
from definitions import BOLTZMANN_CONSTANT, AVOGADROS_NUMBER
import matplotlib.pyplot as plt
import seaborn as sbn
import sklearn.metrics
import scipy
import argparse
import functions
import os
import time
import multiprocessing.pool
import tqdm



def inhibition_to_ddg(ki_a, ki_b, temperature=300.0) -> float:
    """
    Convert experimental inhibition constant (K_i) values to relative binding free energy

    Parameters:
    -----------
    ki_a: float
        experimental K_i of ligand a
    ki_b: float
        experimental K_i of ligand b 
    temperature: float
        temperature in kelvin

    Return:
    -------
    float
        relative binding free energy
    """

    ic50_a = 2 * ki_a
    ic50_b = 2 * ki_b

    return (BOLTZMANN_CONSTANT * AVOGADROS_NUMBER * temperature / 4184) * np.log(ic50_b / ic50_a)


def get_experimental_error(error_a, ki_a, error_b, ki_b, temperature=300):
    """
    Propagate experimental error to obtain error on experimental relative binding free energy

    Parameters:
    -----------
    error_a: float
        experimentall error in K_i for ligand a
    ki_a: float
        experimental K_i of ligand a
    error_b: float
        experimental error in K_i for ligand b
    ki_b: float
        experimental K_i of ligand b

    Return:
    -------
    : _type_
        _description_
    """
    fraction = ki_b / ki_a
    fraction_error = fraction * np.sqrt((error_b / ki_b) ** 2 + (error_a / ki_a) ** 2)
    return (BOLTZMANN_CONSTANT * temperature * fraction_error / fraction) * AVOGADROS_NUMBER / 4184


def bootstrap_statistics(experimental: np.array, calculated: np.array, n_samples = 10000, alpha_level = 0.05):
    """
    _summary_

    Parameters:
    -----------
    experimental: np.array
        _description_
    calculated: np.array
        _description_
    n_samples: int
        number of bootstrapping samples
    alpha_level: float
        confidence level
    Return:
    -------
    results: dict
        dictionary of Pearson's r, MUE and Spearman's r coefficient with the real value, 
        bootsrapped mean and bootstrapped lower and upper bounds
    """
    n_data_samples = len(experimental)
    statistics_dict = {"r": [],
                       "mue": [],
                       "rho": []}

    for i in range(n_samples):
        if i==0:
            experimental_samples = experimental
            calculated_samples = calculated
        else:
            bootstrap_sample = np.random.choice(range(n_data_samples), size = n_samples)
            experimental_samples = [experimental[i] for i in bootstrap_sample]
            calculated_samples = [calculated[i] for i in bootstrap_sample]

        pearson, _ = scipy.stats.pearsonr(experimental_samples, calculated_samples)
        mue = sklearn.metrics.mean_absolute_error(experimental_samples, calculated_samples)
        spearman = scipy.stats.spearmanr(experimental_samples, calculated_samples)
        statistics_dict["r"].append(pearson)
        statistics_dict["mue"].append(mue)
        statistics_dict["rho"].append(spearman)

    results = {"r": {},
               "mue": {},
               "rho": {}}
    
    lower_fraction = alpha_level/2.0
    upper_fraction = 1 - lower_fraction

    for statistic in statistics_dict.keys():
        results[statistic]["real"] = statistics_dict[statistic][0]
        statistics_dict[statistic] = sorted(statistics_dict[statistic])
        results[statistic]["mean"] = np.mean(statistics_dict[statistic])
        results[statistic]["lower"] = statistics_dict[statistic][int(n_samples * lower_fraction)]
        results[statistic]["upper"] = statistics_dict[statistic][int(n_samples * upper_fraction)]
    return results


def create_correct_header(minimisation_simfile, unbound_simfile, bound_simfile, template_header):

    with open(minimisation_simfile, "r") as file:
        minimisation_lines = file.readlines()

    with open(minimisation_simfile, "r") as file:
        for j, line in enumerate(file):
            if "lambda" in line:
                start = j
            if "#" not in line:
                end = j
                break

    lambda_lines = minimisation_lines[start:end]
    correct_lambda_header = template_header + lambda_lines

    write_header(unbound_simfile, correct_lambda_header)
    write_header(bound_simfile, correct_lambda_header)


def write_header(simfile, correct_header):
    """
    Take the new, correct header and append it to the start of the simfile

    Parameters:
    -----------
    simfile: str
        full path to the simfile in the specific lambda window
    correct_header: list
        list of the correct header lines for the specific lambda window

    Return:
    -------
    """
    with open(simfile, "r+") as file:
        data = file.readlines()
        file.seek(0, 0)
        header_and_data = correct_header + data
        file.writelines(header_and_data)


def read_protocol(file):
    """
    Read the protocol.dat file in to a dictionary 

    Parameters:
    -----------
    file: str
        full path to the protocol.dat file written by meze.py

    Return:
    -------
    protocol: dict
        protocol file as a dictionary
    """
    protocol = {}
    with open(file, "r") as f:
        for line in f: 
            key, value = line.rstrip().split(" = ")
            protocol[key] = value
    return protocol


def fix_simfile(protocol): # ONLY FOR SOMD
    """
    Get the headers from minimisation simfile for each lambda file, which contain the correct headers required by MBAR.

    Parameters:
    -----------
    protocol: dict
        the protocol file as a dictionary

    Return:
    -------
    """
    # The minimisation simfile for each lambda file contains the correct headers required by MBAR
    #TODO: How do I set the path to this file?
    template = "/home/jguven/projects/metalloenzymes/meze/simfile_header.txt"
    with open(template, "r") as file:
        template_header = file.readlines()

    outputs = functions.path_exists(protocol["outputs"])

    engine = protocol["engine"]
    repeat_paths = functions.read_files(outputs + "/" + engine + "_*/*/")

    for path in repeat_paths:

        unbound_directory = path + "unbound/"
        bound_directory = path + "bound/"

        unbound_minimisation_simfiles = functions.read_files(unbound_directory + "/minimisation/lambda_*/simfile.dat")
        unbound_simfiles = functions.read_files(unbound_directory + "lambda_*/simfile.dat")
        bound_simfiles = functions.read_files(bound_directory + "lambda_*/simfile.dat")

        arguments = [(unbound_minimisation_simfiles[i], unbound_simfiles[i], bound_simfiles[i], template_header) for i in range(len(unbound_minimisation_simfiles))]

        # with multiprocessing.pool.Pool() as pool:
        #     pool.starmap(create_correct_header, arguments)
            
        for i in range(len(unbound_minimisation_simfiles)):

            with open(unbound_minimisation_simfiles[i], "r") as file:
                minimisation_simfile = file.readlines()
            
            with open(unbound_minimisation_simfiles[i], "r") as file:
                for j, line in enumerate(file):
                    if "lambda" in line:
                        start = j
                    if "#" not in line:
                        end = j
                        break

            lambda_lines = minimisation_simfile[start:end]
            correct_lambda_header = template_header + lambda_lines

            write_header(unbound_simfiles[i], correct_lambda_header)
            write_header(bound_simfiles[i], correct_lambda_header)
            

def get_results(protocol):
    
    engine = protocol["engine"]
    outputs = protocol["outputs"]

    results_file = outputs + engine + "_results.txt"
    print(results_file)

    transformation_paths = functions.read_files(outputs + engine + "_*/*/")
    
    



def write_results(path):

    engine = path.split("/")[-3].split("_")[0]

    
    

    









        

def main():
    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")

    parser.add_argument("-pf",
                        "--protocol-file",
                        dest="protocol",
                        default=os.getcwd() + "/afe/protocol.dat",
                        help="input pdb file for the metalloenzyme/protein")

    arguments = parser.parse_args()
    protocol_file = functions.file_exists(arguments.protocol)
    protocol = read_protocol(protocol_file)
    # s = time.time()
    # # fix_simfile(protocol)
    # print(f"fix_simfile took {time.time() - s} seconds")
    get_results(protocol)

if __name__ == "__main__":
    main()

