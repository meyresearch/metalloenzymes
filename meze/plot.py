import warnings
import matplotlib
import numpy as np
import scipy.stats 
import sklearn.metrics
import pandas as pd
from definitions import BOLTZMANN_CONSTANT, AVOGADROS_NUMBER, COLOURS
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import meze
import functions


def read_results(protocol):
    """
    Read in dataframes of the raw results from each repeat.

    Parameters:
    -----------
    protocol: dict
        protocol file as a dictionary

    Return:
    -------
    free_energies, errors: tuple(np.array, np.array)
        free energies and errors from all repeats: [[repeat_1, repeat_2, repeat_3]...]
    """
    outputs = functions.path_exists(protocol["outputs"])
    engine = protocol["engine"]
    n_repeats = functions.check_int(protocol["repeats"])

    free_energies = []
    errors = []

    for i in range(1, n_repeats + 1):
        results_file = f"{outputs}/{engine}_{i}_raw.csv"
        dataframe = pd.read_csv(results_file, index_col=False).sort_values(by=["transformation"])
        transformations = dataframe["transformation"].tolist()
        free_energy = dataframe["free-energy"].to_numpy()
        error = dataframe["error"].to_numpy()
        free_energies.append(free_energy)
        errors.append(error)
    
    return transformations, np.array(free_energies).T, np.array(errors).T


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
    calculated_values = np.delete(results["average_ddg"].to_numpy(), remove_indices)

    pearson_r = scipy.stats.pearsonr(experimental_values, calculated_values)
    spearman = scipy.stats.spearmanr(experimental_values, calculated_values)

    mue = sklearn.metrics.mean_absolute_error(experimental_values, calculated_values)
    rmse = sklearn.metrics.root_mean_squared_error()
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
    # print(f"Pearson R:                                               {pearson_r[0]:.3f}")
    # print(f"                                          p value:   {pearson_r[1]:.3E}")
    # print(f"Bootstrapped statistics:")
    # print("\n")
    # print(f"                      Mean:        {stats['pearson_r']['mean_value']:.3f}") 
    # print(f"                   Lower Bound:    {stats['pearson_r']['lower_bound']:.3f}")
    # print(f"                   Upper bound:    {stats['pearson_r']['upper_bound']:.3f}")
    # print("\n")
    print("--------------------------------------------------------------")
    # print(f"Spearman rho:                                            {spearman[0]:.3f}")
    # print(f"                                          p value:   {spearman[1]:.3E}")
    # print(f"Bootstrapped statistics:")
    # print("\n")
    # print(f"                      Mean:        {stats['spearman_rho']['mean_value']:.3f}") 
    # print(f"                   Lower Bound:    {stats['spearman_rho']['lower_bound']:.3f}")
    # print(f"                   Upper bound:    {stats['spearman_rho']['upper_bound']:.3f}")
    # print("\n")
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
    # statistics_dataframe["p_values"] = [pearson_r[1], "", spearman[1]]
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


def combine_results(protocol):
    """
    Take values and errors from individual runs and average them and propagate errors.

    Parameters:
    -----------
    protocol: dict
        protocol file as a dictionary

    Return:
    -------
    results: pd.DataFrame
        dataframe containing averaged results and propagated errors
    """

    transformations, free_energies, errors = read_results(protocol)
    
    nan_indices = functions.check_nan(free_energies)
    calculated_means = functions.average(free_energies) 
    propagated_errors = functions.add_in_quadrature(errors) 
    standard_deviations = functions.standard_deviation(free_energies)

    results_dictionary = {"transformation": transformations, 
                          "average_ddg": calculated_means,
                          "standard_deviation": standard_deviations,
                          "propagated_error": propagated_errors}
    results = pd.DataFrame.from_dict(results_dictionary)
    results_file = protocol["outputs"] + "/" + protocol["engine"] + "_results.csv"
    results.to_csv(results_file, index=False)

    if len(nan_indices) > 0: 
        write_results_warning(nan_indices, results_file, transformations)
        warnings.warn(f"Some transformations contained NaNs. Please check {results_file} for details.")

    return results


def write_results_warning(nan_indices, results_file, transformations):
    """
    Write a warning message to the results dataframe csv file about transformations containgin NaNs.

    Parameters:
    -----------
    nan_indices: array
        array of indices of values that are NaNs
    results_file: str
        name and full path to the results.csv file
    transformations: list
        list of transformation names

    Return:
    -------
    """

    with open(results_file, "r+") as file:
        results = file.readlines()
        file.seek(0, 0)
        file.write("##############################################################\n")
        file.write("#                           Warning                          #\n")
        file.write("#                                                            #\n")
        file.write("#             The following lines contained NaNs:            #\n")
        file.write("#                                                            #\n")
        for line in nan_indices:
            transformation_index = line[0]
            repeat_index = line[1]
            line_width = len("############################################################")
            pad = "".rjust(12)
            middle_pad_length = line_width - 2 * len(pad) - len(transformations[transformation_index]) - len(f"repeat: {repeat_index}")
            middle_pad = "".ljust(middle_pad_length)
            file.write(f"#{pad}{transformations[transformation_index]}{middle_pad}repeat: {repeat_index}{pad}#\n")
        file.write("#                                                            #\n")
        file.write("##############################################################\n")    
        file.writelines(results)



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


def plot_individual_runs(protocol, experimental_free_energies_with_nans, experimental_errors_with_nans, results):
    """
    Plot correlation of individual AFE runs against experimental values

    Parameters:
    -----------
    protocol: dict
        protocol file as a dictionary
    experimental_free_energy: np.array
        experimental free energy values (converted from inhibition constants)
    experimental_error: np.array
        propagated error on the experimental free energy values
    results: 
        averaged free energy values and propagated errors

    Return:
    -------
    """
    means_from_file = results["average_ddg"].to_numpy() 
    calculated_nan_indices = functions.check_nan(means_from_file).flatten()

    errors_from_file = results["standard_deviation"].to_numpy()
    calculated_std_nan_indices = functions.check_nan(errors_from_file).flatten()

    experimental_nan_indices = functions.check_nan(experimental_free_energies_with_nans).flatten()
    experimental_error_nan_indices = functions.check_nan(experimental_errors_with_nans).flatten()

    calculated_free_energies_with_nans = means_from_file.copy()
    for i in experimental_nan_indices:
        calculated_free_energies_with_nans[i] = np.nan
    
    experimental_free_energy_array_with_nans = experimental_free_energies_with_nans.copy()
    for i in calculated_nan_indices:
        experimental_free_energy_array_with_nans[i] = np.nan
    
    calculated_std = errors_from_file.copy()
    for i in experimental_error_nan_indices:
        calculated_std[i] = np.nan
    
    experimental_errors_array = experimental_errors_with_nans.copy()
    for i in calculated_std_nan_indices:
        experimental_errors_array[i] = np.nan
    
    calculated_free_energies = calculated_free_energies_with_nans[~np.isnan(calculated_free_energies_with_nans)]
    calculated_errors = calculated_std[~np.isnan(calculated_std)]
    experimental_free_energy_array_with_nans = np.array(experimental_free_energy_array_with_nans)
    experimental_free_energies = experimental_free_energy_array_with_nans[~np.isnan(experimental_free_energy_array_with_nans)]
    experimental_errors_array = np.array(experimental_errors_array)
    experimental_errors = experimental_errors_array[~np.isnan(experimental_errors_array)]

    fig, ax = plt.subplots(1, 1, figsize=(16, 20))
    sns.set(context="notebook", palette="colorblind", style="ticks")
    outputs = protocol["outputs"]
    plots = protocol["plots directory"]
    repeats = functions.check_int(protocol["repeats"])
    engine = protocol["engine"]

    for i in range(1, repeats + 1):
        raw_datafile = outputs + "/" + engine + f"_{i}_raw.csv"
        dataframe = pd.read_csv(raw_datafile, index_col=False).sort_values(by=["transformation"])

        free_energies = dataframe["free-energy"].to_numpy()
        individual_run_ddg = free_energies.copy()

        individual_run_nan_indices = functions.check_nan(individual_run_ddg).flatten()

        for j in experimental_nan_indices:
            individual_run_ddg[j] = np.nan
        
        experimental_free_energies_for_individual_run = experimental_free_energies_with_nans.copy()
        for j in individual_run_nan_indices:
            experimental_free_energies_for_individual_run[j] = np.nan

        ddg = individual_run_ddg[~np.isnan(individual_run_ddg)]
        experimental_array = np.array(experimental_free_energies_for_individual_run)
        experimental_free_energies_for_individual_run = experimental_array[~np.isnan(experimental_array)]

        ax.scatter(experimental_free_energies_for_individual_run, ddg, s=50, label=f"{engine} {i}")


    (_, _, _) = plt.errorbar(experimental_free_energies,
                             calculated_free_energies,
                             color="#D0006F",
                             xerr=experimental_errors,
                             yerr=calculated_errors,
                             capsize=3,
                             linestyle="",
                             zorder=-1)
    

    max_calculated = max(np.absolute(calculated_free_energies)) + 1 
    max_experimental = max(np.absolute(experimental_free_energies)) + 1
    max_y = max(max_calculated, max_experimental)

    ax.plot([-max_y, max_y], [-max_y, max_y], color="#0099AB", linestyle=":", zorder=-1)
    ax.set_xlabel("$\Delta \Delta$ G$_\mathrm{EXP}$ (kcal mol⁻¹)", fontsize=14)
    ax.set_ylabel("$\Delta \Delta$ G$_\mathrm{AFE}$ (kcal mol⁻¹)", fontsize=14)
    ax.vlines(0, -max_y, max_y, color = "silver", linestyle="--", zorder=-1)
    ax.hlines(0, -max_y, max_y, color = "silver", linestyle="--", zorder=-1)
    ax.set_xlim(-max_y, max_y)
    ax.set_ylim(-max_y, max_y)

    labels = [transformation.strip().replace("_", "").replace("ligand", "").replace("~", " to ") for transformation in results["transformation"].tolist()]
    for i in experimental_nan_indices:
        labels[i] = np.nan
    for i in calculated_nan_indices:
        labels[i] = np.nan

    labels = [value for value in labels if str(value).lower() != "nan"]
    for i in range(len(labels)):
        ax.annotate(labels[i], (experimental_free_energies[i], calculated_free_energies[i]), fontsize=9)
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{plots}/individual_correlation.png", dpi=1000)


def plot_correlation(plots_directory, results, experimental_free_energies_with_nans, experimental_errors_with_nans, region=True):
    """
    Plot the correlation plot of experimental binding free energies vs calculated free energies

    Parameters:
    -----------
    plots_directory: str
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
    means_from_file = results["average_ddg"].to_numpy() 
    calculated_nan_indices = functions.check_nan(means_from_file).flatten()

    errors_from_file = results["standard_deviation"].to_numpy()
    calculated_std_nan_indices = functions.check_nan(errors_from_file).flatten()

    experimental_nan_indices = functions.check_nan(experimental_free_energies_with_nans).flatten()
    experimental_error_nan_indices = functions.check_nan(experimental_errors_with_nans).flatten()

    calculated_free_energies_with_nans = means_from_file.copy()
    for i in experimental_nan_indices:
        calculated_free_energies_with_nans[i] = np.nan
    
    experimental_free_energy_array_with_nans = experimental_free_energies_with_nans.copy()
    for i in calculated_nan_indices:
        experimental_free_energy_array_with_nans[i] = np.nan
    
    calculated_std = errors_from_file.copy()
    for i in experimental_error_nan_indices:
        calculated_std[i] = np.nan
    
    experimental_errors_array = experimental_errors_with_nans.copy()
    for i in calculated_std_nan_indices:
        experimental_errors_array[i] = np.nan
    
    calculated_free_energies = calculated_free_energies_with_nans[~np.isnan(calculated_free_energies_with_nans)]
    calculated_errors = calculated_std[~np.isnan(calculated_std)]
    experimental_free_energy_array_with_nans = np.array(experimental_free_energy_array_with_nans)
    experimental_free_energies = experimental_free_energy_array_with_nans[~np.isnan(experimental_free_energy_array_with_nans)]
    experimental_errors_array = np.array(experimental_errors_array)
    experimental_errors = experimental_errors_array[~np.isnan(experimental_errors_array)]


    fig, ax = plt.subplots(1, 1, figsize=(10,10))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)

    ax.scatter(experimental_free_energies, 
               calculated_free_energies, 
               s=50, 
               color=COLOURS["PINK"])

    ax.errorbar(experimental_free_energies,
                calculated_free_energies,
                color=COLOURS["PINK"],
                yerr=calculated_errors,
                xerr=experimental_errors,
                capsize=3,
                linestyle="",
                zorder=-1)

    max_calculated = max(np.absolute(calculated_free_energies)) + 1 
    max_experimental = max(np.absolute(experimental_free_energies)) + 1
    max_y = max(max_calculated, max_experimental)

    ax.plot([-max_y, max_y], [-max_y, max_y], color=COLOURS["BLUE"], linestyle=":", zorder=-1)
    ax.vlines(0, -max_y, max_y, color="silver", linestyle="--", zorder=-1)
    ax.hlines(0, -max_y, max_y, color="silver", linestyle="--", zorder=-1)

    if region:
        top = np.arange(-max_y+0.5, max_y+1.5)
        bottom = np.arange(-max_y-0.5, max_y+0.5)
        x = np.arange(-max_y, max_y+1)
        ax.fill_between(x, bottom, top, alpha=0.2, zorder=-1)

    ax.set_xlim(-max_y, max_y)
    ax.set_ylim(-max_y, max_y)
    ax.set_xlabel("$\Delta \Delta$ G$_\mathrm{EXP}$ (kcal mol \u207B \u00B9)")
    ax.set_ylabel("$\Delta \Delta$ G$_\mathrm{AFE}$ (kcal mol \u207B \u00B9)")
    fig.tight_layout()
    fig.savefig(f"{plots_directory}/meze_AFE_correlation.png", dpi=1000)
    

def plot_bar(plots_directory, afe_df, exp_free_energy, exp_error):
    """
    Plot the bar plot of experimental binding free energies vs calculated free energies

    Parameters:
    -----------
    plots_directory: str
        full path to plots directory
    afe_df: pd.DataFrame
        calculated free energy values and errors 
    exp_free_energy: np.array
        experimental free energy values (converted from inhibition constants)
    exp_error: np.array
        propagated error on the experimental free energy values 

    Return:
    -------

    """
    means = afe_df["average_ddg"].to_numpy() 
    std = afe_df["standard_deviation"].to_numpy()
    transformations = afe_df["transformation"].to_list()
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
    fig.savefig(f"{plots_directory}/meze_AFE_barplot.png", dpi=1000) 


def plot_rmsd_box_plot(protocol):
    """
    Read in RMSD from each lambda window in unbound & bound stages and plot in a single box plot for each transfor

    Parameters:
    -----------
    protocol: dict
        protocol file as a dictionary

    Return:
    -------
    """
    outputs = protocol["outputs"]
    engine = protocol["engine"]
    repeats = functions.check_int(protocol["repeats"])
    plots = protocol["plots directory"]

    for i in range(1, repeats + 1):

        path = outputs + "/" + engine + f"_{i}/"
        transformation_directories = functions.get_files(path + "ligand_*/")

        for j in range(len(transformation_directories)):
            transformation = transformation_directories[j].split("/")[-2]
            plot = plots + f"rmsd_repeat_{i}_{transformation}.png"
            if not os.path.isfile(plot):
                unbound_rmsd_files = functions.get_files(transformation_directories[j] + "unbound/lambda_*/*.npy")
                unbound_rmsds = [np.load(file)[1] for file in unbound_rmsd_files]

                bound_lambda_directories = functions.get_files(transformation_directories[j] + "/bound/lambda_*")
                lambda_windows = np.round(np.arange(0, len(bound_lambda_directories)/10, 0.1), 1)
                bound_rmsd_files = functions.get_files(transformation_directories[j] + "bound/lambda_*/*.npy")
                bound_rmsds = [np.load(file)[1] for file in bound_rmsd_files]

                fig, ax = plt.subplots(1, 1, figsize=(10, 10))
                sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2)
                median_line_properties = dict(linestyle="-", linewidth=2.4, color="k")
                xtick_positions = np.arange(1, len(lambda_windows)+1, 1)

                unbound_boxes = ax.boxplot(unbound_rmsds, positions=xtick_positions - 0.25/2, widths=0.25, patch_artist=True, medianprops=median_line_properties, manage_ticks=False)
                bound_boxes = ax.boxplot(bound_rmsds, positions=xtick_positions + 0.25/2, widths=0.25, patch_artist=True, medianprops=median_line_properties, manage_ticks=False)
                
                for patch in unbound_boxes["boxes"]:
                    patch.set_facecolor(COLOURS["RGBA_WHITE"])
                for patch in bound_boxes["boxes"]:
                    patch.set_facecolor(COLOURS["RGBA_PINK"])
                
                ax.set_xticks(ticks=xtick_positions, labels=lambda_windows)
                ax.legend(handles=[unbound_boxes["boxes"][0], bound_boxes["boxes"][1]], labels=["Unbound", "Bound"], frameon=False)
                ax.set_xlabel("$\lambda$ window")
                ax.set_ylabel("RMSD ($\AA$)")
                fig.tight_layout()
                fig.subplots_adjust(wspace=0.05)
                
                fig.savefig(plot, dpi=1000)
                plt.close(fig)


def plot_pairwise_lambda_rmsd(protocol):
    """
    Read in the pairwise lambda RMSD for bound and unbound and plot as heatmap.

    Parameters:
    -----------
    protocol: dict
        protocol file as a dictionary

    Return:
    -------
    """
    outputs = protocol["outputs"]
    engine = protocol["engine"]
    repeats = functions.check_int(protocol["repeats"])
    
    for i in range(1, repeats + 1):
        path = outputs + "/" + engine + f"_{i}/"
        transformation_directories = functions.get_files(path + "ligand_*/")

        for j in range(len(transformation_directories)):
            transformation_directory = transformation_directories[j]
            
            rmsd_files = functions.get_files(transformation_directory + "pairwise_*.npy")
            filenames = [functions.get_filename(file) for file in rmsd_files]
            pairwise_rmsds = [np.load(file) for file in rmsd_files]

            maximum_values = [array[np.unravel_index(array.argmax(), array.shape)[0], np.unravel_index(array.argmax(), array.shape)[1]] for array in pairwise_rmsds]
            transformation_max = max(maximum_values)

            for k in range(len(pairwise_rmsds)):
                plot_file = transformation_directory + filenames[k] + f"_repeat_{i}.png"
                if not os.path.isfile(plot_file):
                    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
                    sns.set(context="notebook", style="ticks", font_scale=2)
                    sns.heatmap(ax=ax, 
                                data=pairwise_rmsds[k], 
                                cmap="viridis", 
                                vmin=0, 
                                vmax=transformation_max, 
                                square=True, 
                                cbar_kws={"fraction": 0.460, "pad": 0.04, "label": r"Pairwise RMSD $\AA$"})
                    ax.xaxis.tick_top()
                    ax.tick_params(axis="y", rotation=360)
                    ax.set_title(r"$\lambda$ index")
                    ax.set_ylabel(r"$\lambda$ index")
                    fig.savefig(plot_file, dpi=1000)
                    plt.close(fig)
    

def plot_overlap_matrix(protocol):
    """
    Open overlap matrix arrays and plot overlap for each transformation in bound and unbound stages

    Parameters:
    -----------
    protocol: dict
        protocol file as a dictionary

    Return:
    -------
    """
    outputs = protocol["outputs"]
    engine = protocol["engine"]
    repeats = functions.check_int(protocol["repeats"])

    for i in range(1, repeats + 1):    
        path = outputs + "/" + engine + f"_{i}/"
        transformation_directories = functions.get_files(path + "ligand_*/")

        for j in range(len(transformation_directories)):
            transformation_directory = transformation_directories[j]

            overlap_matrix_files = functions.get_files(transformation_directory + "*_overlap_matrix.npy")
            filenames = [functions.get_filename(file) for file in overlap_matrix_files]
            overlap_matrices = [np.load(file) for file in overlap_matrix_files]

            colour_map = matplotlib.colors.ListedColormap(["#FBE8EB","#68246D","#61BF1A", "#154734"])
            n_colours = colour_map.N
            boundary_values = [0.0, 0.025, 0.1, 0.3, 0.8]
            norm_colours = matplotlib.colors.BoundaryNorm(boundary_values, n_colours, clip=False)
            colour_bar_args = dict(ticks=[0.025, 0.1, 0.3, 0.8],
                                   shrink=0.815)
            
            for k in range(len(overlap_matrix_files)):
            
                plot_file = transformation_directory + filenames[k] + ".png"
                if not os.path.isfile(plot_file):
                    try:    
                        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
                        sns.set(context="notebook", style="ticks", font_scale=2)  
                        sns.heatmap(ax=ax,
                                    data=overlap_matrices[k], 
                                    annot=True, 
                                    fmt=".1f", 
                                    linewidths=0.3, 
                                    annot_kws={"size": 14}, 
                                    square=True, 
                                    robust=True, 
                                    cmap=colour_map,
                                    norm=norm_colours, 
                                    cbar_kws=colour_bar_args,
                                    vmax=1)
                        ax.xaxis.tick_top()
                        ax.tick_params(axis="y", rotation=360)
                        ax.set_title(r"$\lambda$ index")
                        ax.set_ylabel(r"$\lambda$ index")
                        fig.savefig(plot_file, dpi=1000)
                        plt.close(fig)

                    except IndexError as e:
                        print(f"transformation {transformation_directory.split('/')[-2]} at repeat {i} raised error: {e}")


def main():

    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")

    parser.add_argument("protocol_file",
                        help="protocol file containing equilibration options",
                        type=str,
                        default=os.getcwd() + "/afe/protocol.dat")

    
    parser.add_argument("experimental_file",
                        default=os.getcwd() + "/afe/experimental_K_i.csv",
                        help="file containing experimental inhibition constants and errors")
    
    parser.add_argument("-s",
                        "--separator",
                        help="character separating the two ligand names",
                        default="~",
                        type=meze.character)
    
    arguments = parser.parse_args()
    protocol_file = functions.file_exists(arguments.protocol_file)
    protocol = functions.read_protocol(protocol_file)

    plot_rmsd_box_plot(protocol)
    plot_pairwise_lambda_rmsd(protocol)
    plot_overlap_matrix(protocol)

    transformations, _, _ = read_results(protocol) 
    results = combine_results(protocol)
    experimental_file = arguments.experimental_file
    experimental_free_energy, experimental_error = get_experimental_data(experimental_file, transformations)
    
    plot_bar(protocol["plots directory"], results, experimental_free_energy, experimental_error)
    plot_correlation(protocol["plots directory"], results, experimental_free_energy, experimental_error)
    plot_individual_runs(protocol, experimental_free_energy, experimental_error, results)

    statistics = output_statistics(experimental_free_energy, results)
    save_statistics_to_file(protocol["outputs"], statistics)
    

if __name__ == "__main__":
    main()