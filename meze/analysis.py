import BioSimSpace as bss
import glob
import numpy as np
import pandas as pd
from definitions import BOLTZMANN_CONSTANT, AVOGADROS_NUMBER
import matplotlib.pyplot as plt
import seaborn as sbn
import sklearn.metrics
import scipy


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


def fix_simfile():
    
    #TODO: How do I set the path to this file?
    template = "simfile_header.txt"
    with open(template, "r") as file:
        template_header = file.readlines()

    # How do I get the path to outputs?
    # How do I get the number of repeats?
        # Could I write a file in meze.py that outputs all that information and I just have to read it in?
    

def main():
    pass

if __name__ == "__main__":
    main()

