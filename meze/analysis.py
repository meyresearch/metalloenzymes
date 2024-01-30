import math
import shutil
import numpy as np
import pandas as pd
from definitions import BOLTZMANN_CONSTANT, AVOGADROS_NUMBER, COLOURS
import matplotlib.pyplot as plt
import sklearn.metrics
import scipy
import argparse
import functions
import os
import seaborn as sns
import subprocess as sp
import MDAnalysis as mda
import MDAnalysis.analysis.rms
import warnings
import logging
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)


def read_mbar(mbar_textfile):
    """
    Read in mbar.txt and return binding free energy and an error estimate from MBART

    Parameters:
    -----------
    mbar_textfile: 
        full path to mbar.txt

    Return:
    -------
    free_energy, error: tuple
        values for free energy and its associated error
    """
    with open(mbar_textfile, "r") as file:
        lines = file.readlines()
    result_line = 0
    for i, line in enumerate(lines):
        if "#MBAR free energy difference in kcal/mol" in line:
            result_line = i + 1
    free_energy, error = None, None
    try:
        if "#" in lines[result_line]:
            split_line = lines[result_line].split("#")
            warnings.warn(split_line[-1])
            try:
                free_energy = float(split_line[0].split(",")[0])
                error = float(split_line[0].split(",")[-1])
            except ValueError as error_message:
                print(f"Error: {error_message}")
    except IndexError as e:
        print(f"{e}: {mbar_textfile}")

    else:
        try:
            free_energy = float(lines[result_line].split(",")[0])
            error = float(lines[result_line].split(",")[-1])
        except ValueError as error_message:
            print(f"Error: {error_message}")  
   
    return free_energy, error 



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
    statistics_dict = {"pearson_r": [],
                       "mue": [],
                       "spearman_rho": []}

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
        spearman, _ = scipy.stats.spearmanr(experimental_samples, calculated_samples)
        statistics_dict["pearson_r"].append(pearson)
        statistics_dict["mue"].append(mue)
        statistics_dict["spearman_rho"].append(spearman)

    results = {"pearson_r": {},
               "mue": {},
               "spearman_rho": {}}
    
    lower_fraction = alpha_level/2.0
    upper_fraction = 1 - lower_fraction

    for statistic in statistics_dict.keys():
        results[statistic]["real"] = statistics_dict[statistic][0]
        statistics_dict[statistic] = sorted(statistics_dict[statistic])
        results[statistic]["mean_value"] = np.mean(statistics_dict[statistic])
        results[statistic]["lower_bound"] = statistics_dict[statistic][int(n_samples * lower_fraction)]
        results[statistic]["upper_bound"] = statistics_dict[statistic][int(n_samples * upper_fraction)]
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
    data = []
    with open(simfile, "r") as file:

        for line in file:
            if "#" not in line and not line.isspace():
                data.append(line)

        # Check last value of datafile
        # If it's not the last frame (#TODO remove hard-coding this), delete the line
        # This is because there seem to be some random numbers at the end of the files and I'm not sure where they are coming from 
        
        for i, line in enumerate(data):
            split_line = line.split()
            if "#" not in line and len(split_line) < 16: 
                del data[i]

        file.seek(0, 0)
        header_and_data = correct_header + data
    with open(simfile, "w") as file:
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


def fix_simfile(protocol): 
    """
    Get the headers from minimisation simfile for each lambda file, which contain the correct headers required by MBAR.
    Note: Only for SOMD

    Parameters:
    -----------
    protocol: dict
        the protocol file as a dictionary

    Return:
    -------
    """
    template = os.environ["MEZEHOME"] + "simfile_header.txt"
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
        
        old_unbound_files = [filename.replace(".dat", "_original") for filename in unbound_simfiles]
        old_bound_files = [filename.replace(".dat", "_original") for filename in bound_simfiles]

        _ = [shutil.copy(unbound_simfiles[i], old_unbound_files[i]) for i in range(len(old_unbound_files))]
        _ = [shutil.copy(bound_simfiles[i], old_bound_files[i]) for i in range(len(old_bound_files))]

        for i in range(len(unbound_minimisation_simfiles)):
            with open(unbound_minimisation_simfiles[i], "r") as file:
                minimisation_simfile = file.readlines()
            if "#" in minimisation_simfile[0]:
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
            else: 
                continue
                

def get_results(protocol):
    """
    Using the protocol, get results from each repeat.
    Save relative binding free energy values and errors into two dataframes.

    Adapted from BioSimSpace.

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

    values_file = outputs + "/" + engine + "_values.csv"
    errors_file = outputs + "/" + engine + "_errors.csv"

    if os.path.isfile(values_file) and os.path.isfile(errors_file):
        values_dataframe = pd.read_csv(values_file).dropna()
        errors_dataframe = pd.read_csv(errors_file).dropna()
    else:

        repeat_paths = functions.read_files(outputs + engine + "_*/")
        get_transformations = functions.read_files(repeat_paths[0] + "*/")
        transformations = [path.split("/")[-2] for path in get_transformations if "minimisations" not in path] # the transformations are the same for each repeat
        values_dictionary = {"transformations": transformations}
        errors_dictionary = {"transformations": transformations}

        for i in range(len(repeat_paths)):
            transformation_paths = functions.read_files(repeat_paths[i] + "ligand_*/")
            repeat = "repeat_"+ str(i + 1)
            values, errors = [], []
            for transformation in transformation_paths:

                unbound = transformation + "unbound/"
                bound = transformation + "bound/"
                # Move minimisation directory out of transformation directory to avoid BioSimSpace error
                try:
                    unbound_minimisation_directory = repeat_paths[i] + "unbound_minimisations"
                    bound_minimisation_directory = repeat_paths[i] + "bound_minimisations"
                    functions.mkdir(unbound_minimisation_directory)
                    functions.mkdir(bound_minimisation_directory)
                    target_path = unbound_minimisation_directory + "/" + transformation.split("/")[-2]
                    source_path = unbound + "minimisation"
                    shutil.move(source_path, target_path)

                    target_path = bound_minimisation_directory + "/" + transformation.split("/")[-2]
                    source_path = bound + "minimisation"
                    shutil.move(source_path, target_path)
                except FileNotFoundError as error_message:
                    print(f"{error_message}: Files have already been moved")
                
                analyser_path = os.environ["BSS_HOME"] + "analyse_freenrg"

                if not os.path.isfile(f"{unbound}/mbar.txt"):

                    try:
                        unbound_command = f"{analyser_path} mbar -i {unbound}/lambda*/simfile.dat -o {unbound}/mbar.txt --overlap --subsampling"
                        sp.check_output(unbound_command, shell=True)
                    except sp.CalledProcessError as error_message:
                        print(error_message.output)
                        print("Trying again without subsampling.")
                        warnings.warn(f"Warning: Disabling subsampling may meen results are unreliable. Please check the unbound transformation {transformation.split('/')[-2]}")
                        unbound_command = f"{analyser_path} mbar -i {unbound}/lambda*/simfile.dat -o {unbound}/mbar.txt --overlap"

                    with open(unbound + "mbar.out", "w") as file:
                        sp.run(unbound_command, shell=True, stdout=file)

                if not os.path.isfile(f"{bound}/mbar.txt"):
                    try:
                        bound_command = f"{analyser_path} mbar -i {bound}/lambda*/simfile.dat -o {bound}/mbar.txt --overlap --subsampling"
                        sp.check_output(bound_command, shell=True)
                    except sp.CalledProcessError as error_message:
                        print(error_message.output)
                        print("Trying again without subsampling.")
                        warnings.warn(f"Warning: Disabling subsampling may meen results are unreliable. Please check the bound transformation {transformation.split('/')[-2]}")
                        bound_command = f"{analyser_path} mbar -i {bound}/lambda*/simfile.dat -o {bound}/mbar.txt --overlap"

                    with open(bound + "mbar.out", "w") as file:
                        sp.run(bound_command, shell=True, stdout=file)

                unbound_free_energy, unbound_error = read_mbar(unbound + "/mbar.txt")
                bound_free_energy, bound_error = read_mbar(bound + "/mbar.txt")

                relative_binding_free_energy = None
                error = None

                try:
                    relative_binding_free_energy = bound_free_energy - unbound_free_energy
                    error = math.sqrt(bound_error ** 2 + unbound_error ** 2)
                except TypeError as error_message:
                    print(f"{transformation.split('/')[-2]}: {error_message}")

                values.append(relative_binding_free_energy)
                errors.append(error)

            values_dictionary[repeat] = values
            errors_dictionary[repeat] = errors

        values_dataframe = pd.DataFrame.from_dict(values_dictionary).dropna()
        values_dataframe.to_csv(outputs + "/" + engine + "_values.csv", index=False)
        errors_dataframe = pd.DataFrame.from_dict(errors_dictionary).dropna()
        errors_dataframe.to_csv(outputs + "/" + engine + "_errors.csv", index=False)

    return values_dataframe, errors_dataframe


def combine_results(protocol, values_dataframe, errors_dataframe):
    """
    Take values and errors from individual runs and average them and propagate errors.

    Parameters:
    -----------
    protocol: dict
        protocol file as a dictionary
    values_dataframe: pd.DataFrame
        dataframe containing free energy estimates from individual runs
    errors_dataframe: 
        dataframe containing error estimates from individual runs

    Return:
    -------
    results: pd.DataFrame
        dataframe containing averaged results and propagated errors
    """
    calculated_means = average(values_dataframe)
    propagated_errors = add_in_quadrature(errors_dataframe)
    standard_deviations = standard_deviation(values_dataframe)
    results = pd.concat([values_dataframe["transformations"], calculated_means, standard_deviations, propagated_errors], axis=1).dropna()
    results.to_csv(protocol["outputs"] + "/" + protocol["engine"] + "_results.csv", index=False)
    return results


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


def plot_bar(outputs, afe_df, exp_free_energy, exp_error):
    """
    Plot the bar plot of experimental binding free energies vs calculated free energies

    Parameters:
    -----------
    outputs: str
        full path to outputs directory
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
    
    (_, _, _) = ax.errorbar(x=calculated_x,
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
    
    (_, _, _) = ax.errorbar(x=experimental_x,
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
    fig.savefig(f"{outputs}/meze_AFE_barplot.png", dpi=1000) 


def plot_correlation(outputs, results, experimental_free_energy, exp_error, region=True):
    """
    Plot the correlation plot of experimental binding free energies vs calculated free energies

    Parameters:
    -----------
    outputs: str
        full path to outputs directory
    results: pd.DataFrame
        averaged calculated free energy values and errors 
    experimenta_free_energy: np.array
        experimental free energy values (converted from inhibition constants)
    experimental_error: np.array
        propagated error on the experimental free energy values 

    Return:
    -------
    
    """
    means = results["average"].to_numpy() 
    std = results["std"].to_numpy()

    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)

    ax.scatter(experimental_free_energy, 
               means, 
               s=50, 
               color=COLOURS["PINK"])

    ax.errorbar(experimental_free_energy,
                means,
                color=COLOURS["PINK"],
                yerr=std,
                xerr=exp_error,
                capsize=3,
                linestyle="",
                zorder=-1)
    max_calculated = max(np.absolute(means)) +1
    max_experimental = max(np.absolute(experimental_free_energy)) + 1
    max_y = max(max_calculated, max_experimental)
    ax.plot([-max_y, max_y], [-max_y, max_y], color=COLOURS["BLUE"], linestyle=":", zorder=-1)
    ax.vlines(0, -max_y, max_y, color="silver", linestyle="--", zorder=-1)
    ax.hlines(0, -max_y, max_y, color="silver", linestyle="--", zorder=-1)

    if region:
        top = np.arange(-max_y+0.5, max_y+1.5)
        bottom = np.arange(-max_y-0.5, max_y+0.5)
        x = np.arange(-max_y, max_y+1)
        ax.fill_between(x, bottom, top, alpha=0.2, zorder=-1)

    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(-3.5, 3.5)
    ax.set_xlabel("$\Delta \Delta$ G$_\mathrm{EXP}$ (kcal mol \u207B \u00B9)")
    ax.set_ylabel("$\Delta \Delta$ G$_\mathrm{AFE}$ (kcal mol \u207B \u00B9)")
    fig.tight_layout()
    fig.savefig(f"{outputs}/meze_AFE_correlation.png", dpi=1000) 


def plot_individual_runs(outputs, experimental_free_energy, experimental_error, afe_values, results):
    """
    Plot correlation of individual AFE runs against experimental values

    Parameters:
    -----------
    outputs: str
        full path to outputs directory
    experimental_free_energy: np.array
        experimental free energy values (converted from inhibition constants)
    experimental_error: np.array
        propagated error on the experimental free energy values
    afe_values: pd.DataFrame
        calculated free energy values and errors
    results: 
        averaged free energy values and propagated errors

    Return:
    -------
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
    
    #TODO don't assume 3 repeats
    ax.scatter(experimental_free_energy, afe_values["repeat_1"].to_numpy(), s=50, label="SOMD 1")
    ax.scatter(experimental_free_energy, afe_values["repeat_2"].to_numpy(), s=50, label="SOMD 2")
    ax.scatter(experimental_free_energy, afe_values["repeat_3"].to_numpy(), s=50, label="SOMD 3")

    (_, _, _) = plt.errorbar(experimental_free_energy,
                                results["average"].to_numpy(),
                                color="#D0006F",
                                xerr=experimental_error,
                                yerr=results["std"].to_numpy(),
                                capsize=3,
                                linestyle="",
                                zorder=-1)

    ax.plot([-4.5, 4.5], [-4.5, 4.5], color="#0099AB", linestyle=":", zorder=-1)
    ax.set_xlabel("$\Delta \Delta$ G$_\mathrm{EXP}$ (kcal mol⁻¹)")
    ax.set_ylabel("$\Delta \Delta$ G$_\mathrm{AFE}$ (kcal mol⁻¹)")
    ax.vlines(0, -3.5, 3.5, color = "silver", linestyle="--", zorder=-1)
    ax.hlines(0, -3.5, 3.5, color = "silver", linestyle="--", zorder=-1)
    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(-3.5, 3.5)
 
    labels = [transformation.strip().replace("_", "").replace("ligand", "").replace("~", " to ") for transformation in afe_values["transformations"].tolist()]
    for i in range(len(labels)):
        ax.annotate(labels[i], (experimental_free_energy[i], results["average"].to_numpy()[i]))
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{outputs}/individual_correlation.png", dpi=1000)
    fig.show()


def output_statistics(experimental_free_energy, results):
    """
    Output statistics in a nice way and show bootstrapped statistics.

    Parameters:
    -----------
    experimental_free_energy: list
        experimental free energy values (converted from inhibition constants)
    results: pd.DataFrame
        averaged calculated free energy values and errors

    Return:
    -------
    statistics_dataframe: pd.DataFrame
        dataframe containing the real statics values and the bootstrapped mean, lower and upper bounds
    """
    x = np.array([[1,2,3,4],
              [2,3,np.nan,5],
              [np.nan,5,2,3]])
    np.argwhere(np.isnan(x))
    
    remove_indices = np.argwhere(np.isnan(experimental_free_energy)).flatten()
    
    experimental_values = np.delete(experimental_free_energy, remove_indices)
    calculated_values = np.delete(results["average"].to_numpy(), remove_indices)

    pearson_r = scipy.stats.pearsonr(experimental_values, calculated_values)
    spearman = scipy.stats.spearmanr(experimental_values, calculated_values)
    mue = sklearn.metrics.mean_absolute_error(experimental_values, calculated_values)
    print("\n")
    print("Bootstrapping statistics...\n")
    print("==============================================================")
    print("|                                                            |")
    print("|                        Statistics                          |")
    print("|                                                            |")
    print("==============================================================")
    stats = bootstrap_statistics(experimental_values, calculated_values)
    print("\n")
    print("--------------------------------------------------------------")
    print(f"Pearson R:                                               {pearson_r[0]:.3f}")
    print(f"                                          p value:   {pearson_r[1]:.3E}")
    print(f"Bootstrapped statistics:")
    print("\n")
    print(f"                      Mean:        {stats['pearson_r']['mean_value']:.3f}") 
    print(f"                   Lower Bound:    {stats['pearson_r']['lower_bound']:.3f}")
    print(f"                   Upper bound:    {stats['pearson_r']['upper_bound']:.3f}")
    print("\n")
    print("--------------------------------------------------------------")
    print(f"Spearman rho:                                            {spearman[0]:.3f}")
    print(f"                                          p value:   {spearman[1]:.3E}")
    print(f"Bootstrapped statistics:")
    print("\n")
    print(f"                      Mean:        {stats['spearman_rho']['mean_value']:.3f}") 
    print(f"                   Lower Bound:    {stats['spearman_rho']['lower_bound']:.3f}")
    print(f"                   Upper bound:    {stats['spearman_rho']['upper_bound']:.3f}")
    print("\n")
    print("--------------------------------------------------------------")
    print(f"Mean unsigned error:                                     {mue:.3f}")
    print("\n")
    print(f"Bootstrapped statistics:")
    print("\n")
    print(f"                      Mean:        {stats['mue']['mean_value']:.3f}") 
    print(f"                   Lower Bound:    {stats['mue']['lower_bound']:.3f}")
    print(f"                   Upper bound:    {stats['mue']['upper_bound']:.3f}")
    print("\n")
    print("--------------------------------------------------------------")
    
    statistics_dataframe = pd.DataFrame.from_dict(stats).transpose()
    statistics_dataframe.rename(columns={0: "statistic"})
    statistics_dataframe["p_values"] = [pearson_r[1], "", spearman[1]]
    return statistics_dataframe



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


def get_experimental_data(file, transformations): 
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


def save_statistics_to_file(outputs, statistics):
    """
    Take the statistics dataframe and save it to the outputs directory as a csv 

    Parameters:
    -----------
    outputs: str
        full path to the outputs directory  
    statistics: pd.DataFrame
        statistics dataframe 

    Return:
    -------
    """
    statistics.to_csv(f"{outputs}/meze_statistics.csv")
    

def compute_rmsd(directory, engine):
    """
    Compute root mean square deviation in a simulation

    Parameters:
    -----------
    directory: 
        _description_
    engine: 
        _description_

    Return:
    -------
    : _type_
        _description_
    """
    
    trajectory_file = functions.read_files(directory + "*.dcd")[0]
    
    with mda.lib.formats.libdcd.DCDFile(trajectory_file) as trajectory:
        frames = [frame for frame in trajectory]
    
    first_frame = frames[0].xyz
    universe = mda.Universe(directory + f"{engine.lower()}.prm7", trajectory_file, topology_format="PARM7")
    reference_universe = mda.Universe(directory + f"{engine.lower()}.prm7", first_frame, topology_format="PARM7")

    ligand = universe.select_atoms("resname LIG")
    reference = reference_universe.select_atoms("resname LIG")

    rmsd = mda.analysis.rms.RMSD(ligand, reference)
    rmsd.run()
    rmsd_result = rmsd.results.T

    time = rmsd_result[1]
    rmsd_values = rmsd_result[2]

    return time, rmsd_values


def plot_individual_rmsd(plots_directory, savename, stage, time, rmsd):
    """
    Plot individual root mean square deviation against time.
    Plots are saved in /outputs/plots/ 

    Parameters:
    -----------
    plots_directory: str
        full path to /outputs/plots/
    savename: str
        name of transformation
    stage: str
        unbound or bound
    time: np.array
        array of time
    rmsd: np.array
        rmsd values over time

    Return:
    -------
    """

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)

    ax.plot(time, rmsd, "k-")

    ax.set_xlabel("Time (ps)", fontsize=14)
    ax.set_ylabel("RMSD ($\AA$)", fontsize=14)

    fig.savefig(f"{plots_directory}/{stage}_{savename}.png", dpi=1000)


def analyse_rmsds():
    pass


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
    fix_simfile(protocol)

    values, errors = get_results(protocol)[0].dropna(), get_results(protocol)[1].dropna()
    results = combine_results(protocol, values, errors)

    experimental_file = arguments.experimental_file
    experimental_free_energy, experimental_error = get_experimental_data(experimental_file, results["transformations"].tolist())

    plot_bar(protocol["outputs"], results, experimental_free_energy, experimental_error)
    plot_correlation(protocol["outputs"], results, experimental_free_energy, experimental_error)
    plot_individual_runs(protocol["outputs"], experimental_free_energy, experimental_error, values, results)
    statistics = output_statistics(experimental_free_energy, results)
    save_statistics_to_file(protocol["outputs"], statistics)

    

if __name__ == "__main__":
    main()

