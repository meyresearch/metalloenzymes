import BioSimSpace as bss
import argparse
import functions
from argparse import RawTextHelpFormatter
import pandas as pd
import sys
from definitions import ADJUST_OPTIONS
import numpy as np


def create_network(ligands, ligand_names, ligand_path):
    """
    Create a transformation network with LOMAP
    Parameters:
    -----------
    ligands: list
        list containing BioSimSpace molecule objects of the ligands
    ligand_names: list
        list of ligand names, can be empty
    ligand_path: str
        path to ligand files
    Return:
    -------
    dictionary
        dict of transformation: score

    """
    transformations, lomap_scores = bss.Align.generateNetwork(ligands, plot_network=True, names=ligand_names, work_dir=ligand_path)
    network_dict = {}
    named_transformations = [(ligand_names[transformation[0]], ligand_names[transformation[1]]) for transformation in transformations]
    
    for transformation, score in zip(named_transformations, lomap_scores):
        network_dict[transformation] = score

    with open(ligand_path+f"/meze_lomap_network.csv", "w") as lomap_out:
        for key, value in network_dict.items():
            lomap_out.write(f"{key}: {value}\n")
    
    return network_dict


def network_to_df(dictionary):
    """
    Convert transformation network from dictionary to Pandas dataframe
    Parameters:
    -----------
    dictionary: dict
        LOMAP network
    Return:
    -------
    pd.DataFrame
        LOMAP network as a df
    """
    dataframe_from_dict = pd.DataFrame.from_dict(dictionary, "index")
    network_dataframe = dataframe_from_dict.reset_index().rename(columns={"index":"transformations",
                                                                          0:"score"})
    network_dataframe.index.name = "index"
    print("The LOMAP-generated perturbation network is:\n")
    print(network_dataframe)
    return network_dataframe


def adjust():
    """
    Ask user if they want to edit the LOMAP network
    Parameters:
    -----------
    Return:
    -------
    bool
        True: adjust network, False: do not adjust
    """
    adjust_network = False
    print("\nDo you want to edit this network? ([y]es/[n]o/[q]uit)")

    while True:
        edit = input("> ")
        if edit.lower().strip() == "yes" or edit.lower().strip() == "y":
            adjust_network = True
            break
        elif edit.lower().strip() == "no" or edit.lower().strip() == "n":
            adjust_network = False
            break
        elif edit.lower().strip() == "quit" or edit.lower().strip() == "q":
            print("Quitting.")
            sys.exit()
        elif edit == "":
            continue
        else: 
            print("Invalid option.")
            continue

    return adjust_network


def print_options():
    """
    Print available user options for adjusting the network
    Parameters:
    -----------
    Return:
    -------
    """
    print("Please choose an option:")
    for option in ADJUST_OPTIONS:
        print(option)
    

def get_user_options():
    """    
    Extract user options
    Parameters:
    -----------
    Return:
    -------
    list
        list of user specified option(s)
    """
    option = input("> ").lower()
    return option.replace(",", "").split()


def delete_transformation(edited_dataframe, options, first_index, last_index):
    """
    Delete transformation as specified by user
    Parameters:
    -----------
    edited_dataframe: pd.DataFrame
        dataframe of the network to be edited
    options: list
        list of user specified option(s)
    first_index: int
        first index of the original dataframe
    last_index: int
        last index of the original dataframe
    Return:
    -------
    pd.DataFrame
        updated network
    """
    try:
        delete_indices = [int(options[i]) for i in range(1, len(options))]
        try:
            edited_dataframe = edited_dataframe.drop(delete_indices, axis=0, inplace=False)
            print("\nThe NEW network is:\n")
            print(edited_dataframe)
        except KeyError:
            print(f"Error: delete index should be between {first_index} and {last_index}")
            pass
    except ValueError:
        print("Delete index should be an integer")
        pass
    return edited_dataframe


def edit_transformation(edited_dataframe, options, first_index, last_index):
    """
    Edit transformation score as specified by user
    Parameters:
    -----------
    edited_dataframe: pd.DataFrame
        dataframe of the network to be edited
    options: list
        list of user specified option(s)
    first_index: int
        first index of the original dataframe
    last_index: int
        last index of the original dataframe
    Return:
    -------
    pd.DataFrame
        updated network
    """
    try: 
        if len(options) != 3:
            raise ValueError
        index = int(options[1])
        score = np.float64(options[2])
        try:
            edited_dataframe.at[index, "score"] = score 
            print("\nThe NEW network is:\n")
            print(edited_dataframe)
        except KeyError:
            print(f"Error: delete index should be between {first_index} and {last_index}")
    except ValueError:
        print("edit existing LOMAP scores in format 'edit index score'")
        pass
    return edited_dataframe


def add_transformation(edited_dataframe, options, last_index):
    """
    Edit transformation score as specified by user
    Parameters:
    -----------
    edited_dataframe: pd.DataFrame
        dataframe of the network to be edited
    options: list
        list of user specified option(s)
    last_index: int
        last index of the original dataframe
    Return:
    -------
    pd.DataFrame
        updated network
    """
    try: 
        if len(options) != 4:
            raise ValueError
        transformation = tuple([options[1], options[2]])
        score = np.float64(options[3])
        added_row = pd.DataFrame([[transformation, score]], columns= ["transformations", "score"], index=[last_index+1])
        edited_dataframe = edited_dataframe.append(added_row, ignore_index=False)
        print("\nThe NEW network is:\n")
        print(edited_dataframe)
    except ValueError:
        print("Error: add transformations in format 'add ligand_1, ligand_2, score'")
    return edited_dataframe


def edit_network(network):

    network = network_to_df(network)   
    first_index = network.first_valid_index()
    last_index = network.last_valid_index()
    
    adjust_network = adjust()
    edited_dataframe = network.copy()

    while adjust_network:
        print_options()
        options = get_user_options()
        n_inputs = len(options)
        if options[0] == "del":
            edited_dataframe = delete_transformation(edited_dataframe, options, first_index, last_index)
        elif options[0] == "edit":
            edited_dataframe = edit_transformation(edited_dataframe, options, first_index, last_index)
        elif options[0] == "add":
            edited_dataframe = add_transformation(edited_dataframe, options, last_index)
            last_index = edited_dataframe.last_valid_index()
        elif options[0] == "s":
            pass
        elif n_inputs == 1 and options.lower() != "q":
            print("Error: expected index")
            continue
        elif options.lower() == "q":
            print("Quitting.")
            sys.exit()
        elif options == "":
            continue
        else:
            print("Invalid option.")
            continue
    
    
def main():

    parser = argparse.ArgumentParser(description="MEZE: MEtalloenZymE FF-builder for alchemistry\nCreate and edit a LOMAP network.",
                                     formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("-l",
                        "--ligand-path",
                        dest="ligand_path",
                        help="path to ligand files")
    arguments = parser.parse_args()

    ligand_path = functions.get_absolute_path(arguments.ligand_path)
    ligand_files = functions.get_ligand_files(ligand_path)
    ligands = [bss.IO.readMolecules(file)[0] for file in ligand_files]
    ligand_names = [functions.get_filenames(filepath) for filepath in ligand_files]

    network_dict = create_network(ligands, ligand_names, ligand_path)


if __name__ == "__main__":
    main()