import BioSimSpace as bss
import numpy as np
import pandas as pd
from definitions import BOLTZMANN_CONSTANT, AVOGADROS_NUMBER, COLOURS
import matplotlib.pyplot as plt
import seaborn as sbn
import sklearn.metrics
import scipy
import argparse
import functions
import os
import time
import tqdm
import seaborn as sns




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


def propagate_experimental_error(error_a, ki_a, error_b, ki_b, temperature=300):
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
    float
        propagated error in kcal / mol
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
    #TODO: How do I set the path to this file?
    template = "/home/jguven/projects/metalloenzymes/meze/simfile_header.txt"
    with open(template, "r") as file:
        template_header = file.readlines()

    outputs = functions.path_exists(protocol["outputs"])

    engine = protocol["engine"]
    repeat_paths = functions.read_files(outputs + "/" + engine + "_*/*/")
    print("\n")
    for path in tqdm.tqdm(repeat_paths, desc="Fixing headers"): 

        unbound_directory = path + "unbound/"
        bound_directory = path + "bound/"

        unbound_minimisation_simfiles = functions.read_files(unbound_directory + "/minimisation/lambda_*/simfile.dat")
        unbound_simfiles = functions.read_files(unbound_directory + "lambda_*/simfile.dat")
        bound_simfiles = functions.read_files(bound_directory + "lambda_*/simfile.dat")

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
    """
    Using the protocol, get results from each repeat.
    Save relative binding free energy values and errors into two dataframes.

    Parameters:
    -----------
    protocol: dic
        protocol file as a dictionary

    Return:
    -------
    values, errors: tuple
        tuple containing the values and errors dataframes
    """
    engine = protocol["engine"]
    outputs = protocol["outputs"]
    repeat_paths = functions.read_files(outputs + engine + "_*/")
    get_transformations = functions.read_files(repeat_paths[0] + "*/")
    transformations = [path.split("/")[-2] for path in get_transformations] # the transformations are the same for each repeat
    values_dictionary = {"transformations": transformations}
    errors_dictionary = {"transformations": transformations}
    print("\n")
    for i in tqdm.tqdm(range(len(repeat_paths)), desc="Getting results"):
        transformation_paths = functions.read_files(repeat_paths[i] + "*/")
        repeat = "repeat_"+ str(i + 1)
        values, errors = [], []
        for transformation in transformation_paths:
            unbound = transformation + "unbound/"
            bound = transformation + "bound/"
            unbound_pmf, unbound_matrix = bss.FreeEnergy.Relative.analyse(unbound)
            bound_pmf, bound_matrix = bss.FreeEnergy.Relative.analyse(bound)
            relative_binding_free_energy, error = bss.FreeEnergy.Relative.difference(bound_pmf, unbound_pmf)
            values.append(relative_binding_free_energy._value)
            errors.append(error._value)
        values_dictionary[repeat] = values
        errors_dictionary[repeat] = errors

    values_dataframe = pd.DataFrame.from_dict(values_dictionary)
    values_dataframe.to_csv(outputs + "/" + engine + "_values.csv", index=False)
    errors_dataframe = pd.DataFrame.from_dict(errors_dictionary)
    errors_dataframe.to_csv(outputs + "/" + engine + "_errors.csv", index=False)
    return values_dataframe, errors_dataframe


def average(values):
    """
    Average the binding free energy values accross repeats

    Parameters:
    -----------
    values: pd.DataFrame
        original dataframe of free energy differences

    Return:
    -------
    df: pd.DataFrame
        averaged binding free energies
    """
    df = pd.DataFrame()
    df["average"] = values.mean(axis=1, numeric_only=True)
    return df

def add_in_quadrature(errors):
    """
    Take the errors in each repeat and propagate by adding them in quadrature 

    Parameters:
    -----------
    errors: pd.DataFrame
        original errors dataframe

    Return:
    -------
    df: pd.DataFrame
        propagated errors 
    """
    df = pd.DataFrame()
    squares = pd.DataFrame()
    repeats = list(errors.drop(columns=["transformations"]).columns)
    square_names = []
    for repeat in repeats:
        column = repeat.split("_")[-1]
        squares[f"square_{column}"] = errors[repeat] ** 2
        square_names.append(f"square_{column}")
    squares["sum"] = squares[square_names].sum(axis=1)
    df["error_in_quadrature"] = (1/len(repeats)) * np.sqrt(squares["sum"])
    return df

def standard_deviation(values):
    """
    Take the standard deviation of the binding free energies accross repeats

    Parameters:
    -----------
    values: pd.DataFrame 
        original dataframe of free energy differences 

    Return:
    -------

    """
    df = pd.DataFrame()
    df["std"] = values.std(axis=1, numeric_only=True)
    return df


def main():
    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")

    parser.add_argument("-pf",
                        "--protocol-file",
                        dest="protocol",
                        default=os.getcwd() + "/afe/protocol.dat",
                        help="input pdb file for the metalloenzyme/protein")

    parser.add_argument("-ef",
                        "--experimental-file",
                        dest="experimental_file",
                        default=os.getcwd() + "/afe/experimental_K_i.csv",
                        help="file containing experimental inhibition constants and errors")
    
    arguments = parser.parse_args()
    protocol_file = functions.file_exists(arguments.protocol)
    protocol = read_protocol(protocol_file)
    # s = time.time()
    # fix_simfile(protocol)
    # print(f"fix_simfile took {time.time() - s} seconds")
    # values, errors = get_results(protocol)

    ### FOR DEBUGGING: ###
    path = "/home/jguven/projects/alchemistry/kpc2/partially_protonated_ligand/test_outputs/"
    values = pd.read_csv(path + "SOMD_values.csv")
    errors = pd.read_csv(path + "SOMD_errors.csv")
    #######

    calculated_means = average(values)
    propagated_errors = add_in_quadrature(errors)
    standard_deviations = standard_deviation(values)
 
    results = pd.concat([values["transformations"], calculated_means, standard_deviations, propagated_errors], axis=1)
    results.to_csv(protocol["outputs"] + "/" + protocol["engine"] + "_results.csv", index=False)

    experimental_file = arguments.experimental_file
    experimental_free_energy, experimental_error = get_experimental_data(experimental_file, values["transformations"].tolist())

    plot_bar(results, experimental_free_energy, experimental_error)
    plot_correlation(results, experimental_free_energy, experimental_error)
    pearson_r = scipy.stats.pearsonr(experimental_free_energy, results["average"].to_numpy())
    spearman = scipy.stats.spearmanr(experimental_free_energy, results["average"].to_numpy())
    mue = sklearn.metrics.mean_absolute_error(experimental_free_energy, results["average"].to_numpy())
    # stats = bootstrap_statistics(experimental_free_energy, results["average"].to_numpy())
    # print(stats)

def plot_repeat_correlation():

    plt.figure(figsize=(10, 10))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
    plt.scatter(experimental_free_energy, values["repeat_1"].to_numpy(), s=50, label="somd 1")
    plt.scatter(experimental_free_energy, values["repeat_2"].to_numpy(), s=50, label="somd 2")
    plt.scatter(experimental_free_energy, values["repeat_3"].to_numpy(), s=50, label="somd 3")
    # plt.scatter(3, 2.5, s=0)

    (_, caps, _) = plt.errorbar(experimental_free_energy,
                            results["average"].to_numpy(),
                            color="#D0006F",
                            yerr=results["std"].to_numpy(),
                            capsize=3,
                            linestyle="",
                            zorder=-1)

    plt.plot([-4.5, 4.5], [-4.5, 4.5], color="#0099AB", linestyle=":", zorder=-1)
    plt.xlabel("$\Delta \Delta$ G$_\mathrm{EXP}$ (kcal mol⁻¹)")
    plt.ylabel("$\Delta \Delta$ G$_\mathrm{AFE}$ (kcal mol⁻¹)")
    plt.vlines(0, -3.5, 3.5, color = "silver", linestyle="--", zorder=-1)
    plt.hlines(0, -3.5, 3.5, color = "silver", linestyle="--", zorder=-1)
    plt.xlim(-3.5, 3.5)
    plt.ylim(-3.5, 3.5)
    plt.tight_layout()
    labels = [transformation.strip().replace("_", "").replace("ligand", "").replace("~", " to ") for transformation in values["transformations"].tolist()]
    for i in range(len(labels)):
        plt.annotate(labels[i], (experimental_free_energy[i], results["average"].to_numpy()[i]))
    # plt.savefig("corr.png", dpi=1200, transparent=True)

    plt.legend()
    plt.savefig("test-CORR-all.png")



def plot_bar(outputs, afe_df, exp_free_energy, exp_error):
    """
    Plot the bar plot of experimental binding free energies vs calculated free energies

    Parameters:
    -----------
    afe_df: pd.DataFrame
        calculated free energy values and errors 
    exp_free_energy: np.array
        experimental free energy values (converted from inhibition constants)
    exp_error: np.array
        propagated error on the experimental free energy values 

    Return:
    -------
    """
    means = afe_df["average"].to_numpy() 
    std = afe_df["std"].to_numpy()
    transformations = afe_df["transformations"].to_list()
    labels = [transformation.strip().replace("_", "").replace("ligand", "").replace("~", " to ") for transformation in transformations]
    n_x_labels = np.arange(len(means))
    bar_width = 0.35
    calculated_x = n_x_labels - (bar_width / 2)
    experimental_x = n_x_labels + (bar_width / 2)

    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)

    ax.bar(x=calculated_x,
           height=means,
           yerr=std,
           width=bar_width,
           label="AFE",
           color=COLOURS["PINK"],
           linewidth=0)
    
    (_, caps, _) = ax.errorbar(x=calculated_x,
                               y=means,
                               yerr=std,
                               capsize=3,
                               linestyle="",
                               color="black")
    
    ax.bar(x=experimental_x,
           height=exp_free_energy,
           width=bar_width,
           yerr=exp_error,
           label="EXP",
           color=COLOURS["BLUE"],
           linewidth=0)
    
    (_, caps, _) = ax.errorbar(x=experimental_x,
                               y=exp_free_energy,
                               yerr=exp_error,
                               capsize=3,
                               linestyle="",
                               color="black")
        
    ax.axhline(0, 0, 1, c="black")
    max_calculated = max(np.absolute(means)) +1
    max_experimental = max(np.absolute(exp_free_energy)) + 1
    max_y = max(max_calculated, max_experimental)
    ax.set_ylim(-max_y, max_y)
    ax.set_xticks(calculated_x, labels, rotation = 85, ha="center")
    ax.legend()
    ax.set_xlabel("Transformation")
    ax.set_ylabel("$\Delta \Delta$ G (kcal mol \u207B \u00B9)")
    fig.tight_layout()
    fig.savefig(f"{outputs}/meze_AFE_barplot.png")


def plot_correlation(afe_df, exp_free_energy, exp_error):
    
    means = afe_df["average"].to_numpy() 
    std = afe_df["std"].to_numpy()

    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)

    ax.scatter(exp_free_energy, 
               means, 
               s=50, 
               color=COLOURS["PINK"])

    ax.errorbar(exp_free_energy,
                means,
                color=COLOURS["PINK"],
                yerr=std,
                xerr=exp_error,
                capsize=3,
                linestyle="",
                zorder=-1)
    max_calculated = max(np.absolute(means)) +1
    max_experimental = max(np.absolute(exp_free_energy)) + 1
    max_y = max(max_calculated, max_experimental)
    ax.plot([-max_y, max_y], [-max_y, max_y], color=COLOURS["BLUE"], linestyle=":", zorder=-1)
    ax.vlines(0, -max_y, max_y, color="silver", linestyle="--", zorder=-1)
    ax.hlines(0, -max_y, max_y, color="silver", linestyle="--", zorder=-1)
    ax.set_xlim(-max_y, max_y)
    ax.set_ylim(-max_y, max_y)
    ax.set_xlabel("$\Delta \Delta$ G$_\mathrm{EXP}$ (kcal mol \u207B \u00B9)")
    ax.set_ylabel("$\Delta \Delta$ G$_\mathrm{AFE}$ (kcal mol \u207B \u00B9)")
    
    fig.tight_layout()

    fig.savefig("test-CORR.png")




def get_ligand_indices(transformations):
    """
    Take the list of transformation strings and extract the ligand numbers.
    The ligand numbers are used as indicies for comparing with experiments.

    Parameters:
    -----------
    transformations: list
        list of transformations between two ligands

    Return:
    -------
    first_indices, second_indices: tuple
        tuple of two lists containing the first and second ligand indices, respectively
    """
    first_indices, second_indices = [], []
    for transformation in transformations:
        clean_line = transformation.strip("\n")
        ligand_1 = clean_line.split("~")[0].split("_")[-1]
        ligand_2 = clean_line.split("~")[1].split("_")[-1]
        first_indices.append(int(ligand_1) - 1)
        second_indices.append(int(ligand_2) - 1)
    return first_indices, second_indices


def get_experimental_data(file, transformations): #SOMETHING IS WRONG HERE! ERRORS WRONG AND THE BARS ARE WRONG
    """
    Use the transformations to get the indices of ligands.
    Convert experimental inhibition constants to free energies 

    Parameters:
    -----------
    file: str
        full path to the experimental datafile
    transformations: list
        list of transformation strings

    Return:
    -------
    free_energies, errors: tuple
        ki converted to free energies and propagated errors in kcal / mol
    """
    df = pd.read_csv(file)
    columns = list(df.columns)
    ki, ki_error = df[columns[1]], df[columns[2]]
    first_indices, second_indices = get_ligand_indices(transformations)

    free_energies, errors = [], []
    for i in range(len(first_indices)):
        i_1, i_2 = first_indices[i], second_indices[i]
        free_energy = inhibition_to_ddg(ki_a=ki[i_1], ki_b=ki[i_2])
        error = propagate_experimental_error(error_a=ki_error[i_1], error_b=ki_error[i_2], ki_a=ki[i_1], ki_b=ki[i_2])
        free_energies.append(free_energy)
        errors.append(error)
    return free_energies, errors


if __name__ == "__main__":
    main()

