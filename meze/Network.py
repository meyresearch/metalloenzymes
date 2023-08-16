import BioSimSpace as bss
bss.setVerbose(True)
import argparse
import functions
from argparse import RawTextHelpFormatter
import pandas as pd
import sys
from definitions import ADJUST_OPTIONS
import numpy as np
import shutil
import os
import logging
import Ligand
import csv
import Protein


def check_charge(value):
    """
    Check that given charge is an integer

    Parameters:
    -----------
    value: any
        ligand charge

    Return:
    -------
    charge: int
        ligand charge converted to integer
    """
    try:
        charge = int(value)
    except ValueError:
        print("Error: charge should be an integer")
    return charge 


def directories_exist(paths):
    """
    Check if LOMAP directories (images, inputs, outputs) already exist in ligand_path.
    Raises a warning in each case.

    Parameters:
    -----------
    paths: list
        list of ligand_path + images/inputs/outputs
        
    Return:
    -------
    exists: bool
        True: LOMAP directories exist, False: LOMAP directories do not exist
    """
    exists = False
    for directory in paths:
        if os.path.exists(directory) and os.path.isdir(directory):
            exists = True
            logging.warning(f"LOMAP output directory {directory} already exist.")
    return exists


def remove_lomap_directories(paths):
    """
    Remove LOMAP directories (images, inputs, outputs) in ligand_path.

    Parameters:
    -----------
    paths: list
        list of ligand_path + images/inputs/outputs

    Return:
    -------
    """
    for directory in paths:
        if os.path.exists(directory) and os.path.isdir(directory):
            print(f"Removing directory {directory}")
            shutil.rmtree(directory)


class Network(object):
    """
    Network class object

    Parameters:
    -----------
    object: 
        _description_

    Return:
    -------
    """
    def __init__(self, path, group_name, protein_file, protein_path, water_model, ligand_ff, protein_ff, ligand_charge, threshold=0.4, n_normal=11, n_difficult=17):
        """
        Class constructor
        """
        self.path = functions.path_exists(path)
        self.ligand_forcefield = ligand_ff
        self.water_model = water_model
        self.files = self.get_files()
        self.ligand_charge = check_charge(ligand_charge)
        self.ligands = [Ligand.Ligand(file) for file in self.files]
        self.ligand_molecules = [ligand.get_molecule() for ligand in self.ligands]
        self.names = [ligand.get_name() for ligand in self.ligands]
        
        self.protein_forcefield = protein_ff
        self.group_name = group_name
        self.protein_file = protein_file
        self.protein_path = protein_path
        self.protein = Protein.Protein(name=self.group_name,
                                       protein_file=self.protein_file,
                                       path=self.protein_path,
                                       forcefield=self.protein_forcefield,
                                       water_model=self.water_model)
        self.protein_water_complex = self.protein.create_complex()
        self.prepared_protein = self.protein.tleap(self.protein_water_complex)
        self.threshold = threshold
        self.n_normal = n_normal
        self.n_difficult = n_difficult
        self.n_ligands = self.get_n_ligands()
        self.bound = [None] * self.n_ligands



    def solvate_meze(self, idx, Network, AFE):
        """
        Solvate unbound and bound systems.

        Parameters:
        -----------
        idx: int
            Ligand index for sorting through Network.names
        Protein: 
            Protein class object
        Network: 
            Network class object
        AFE: 
            AlchemicalFreeEnergy class object

        Return:
        -------
        """
        ligand = Network.ligands[idx]
        names = Network.names
        ligand_number = Network.names[idx].split("_")[-1]
        print(f"Solvating unbound ligand {ligand_number}")
        ligand_parameters = ligand.parameterise(Network.forcefield, Network.charge)
        unbound_box, unbound_box_angles = AFE.create_box(ligand_parameters)
        solvated_ligand = bss.Solvent.solvate(model=Protein.water_model, 
                                                molecule=ligand_parameters, 
                                                box=unbound_box,
                                                angles=unbound_box_angles)


    def create_dictionary(self):
        """
        Create a transformation network with LOMAP

        Return:
        -------
        dictionary
            dict of transformation: score
        """
        lomap_work_directory = self.check_lomap_directory()
        transformations, lomap_scores = bss.Align.generateNetwork(self.ligand_molecules, plot_network=True, names=self.names, work_dir=lomap_work_directory)
        network_dict = {}
        named_transformations = [(self.names[transformation[0]], self.names[transformation[1]]) for transformation in transformations]
        
        for transformation, score in zip(named_transformations, lomap_scores):
            network_dict[transformation] = score

        with open(self.path+f"/meze_network.csv", "w") as lomap_out:
            for key, value in network_dict.items():
                lomap_out.write(f"{key}: {value}\n")
        self.dictionary = network_dict
        return network_dict
    

    def get_files(self):
        """
        Read in all .sdf and .mol2 files in the provided path.

        Return:
        -------
        list
            list of ligand filenames
        """
        # Adapted from dbmol.py:
        ligand_files = functions.read_files(f"{self.path}/*.sdf")
        ligand_files += functions.read_files(f"{self.path}/*.mol2")        
        if len(ligand_files) < 2:
            raise IOError(f"Path {self.path} must contain at least two sdf or mol2 files.")
        return ligand_files


    def get_n_ligands(self):
        """
        Get number of ligands

        Return:
        -------
        len(self.files): int
            number of ligands
        """
        return(len(self.files))


    def create_new_lomap_directory(self):
        """
        Create a lomap directory in ligand path 

        Return:
        -------
        str
            new LOMAP working directory 
        """
        new_lomap_directory = self.path + "/lomap/"
        print(f"Creating new directory {new_lomap_directory}")
        if not os.path.exists(new_lomap_directory):
            os.mkdir(new_lomap_directory)
        else:
            print(f"Directory {new_lomap_directory} already exists. Continuing.")
        return new_lomap_directory


    def check_lomap_directory(self):
        """
        This is work-around to avoid the network plotting failing if 
        the directories already exist. 

        Return:
        -------
        str
            updated working directory to be input into bss.generateNetwork
        """
        lomap_work_directory = self.path
        lomap_names = ["/images/", "/inputs/", "/outputs/"]
        lomap_directories = [self.path + name for name in lomap_names]

        exists = directories_exist(lomap_directories)

        if exists:
            remove_lomap_directories(lomap_directories)


        # while exists:
        #     delete = input("\nDo you want to over-write them? [y]es/[n]o: ").lower()
        #     if delete == "yes" or delete == "y":
        #         remove_lomap_directories(lomap_directories)
        #         break
        #     elif delete == "no" or delete == "n":
        #         self.create_new_lomap_directory()
        #         break
        #     else:
        #         print("Invalid option.")
        #         continue
        return lomap_work_directory


    def set_n_windows(self, lomap_score):
        """
        Set the number of lambda windows for given transformation depending 
        on lomap score

        Parameters:
        -----------
        lomap_score: float
            lomap score for transformation

        Return:
        -------
        n_windows: int
            number of lambda windows for given 
        """
        if lomap_score == None or lomap_score < float(self.threshold):
            self.n_windows = self.n_difficult
        else:
            self.n_windows = self.n_normal
        return self.n_windows


    def create_lambda_list_bash(self):
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
        lambda_list_numpy = list(np.linspace(0, 1, int(self.n_windows)))
        lambda_list = [format(item, ".4f") for item in lambda_list_numpy]
        return " ".join(lambda_list)


    def create_network_files(self, engine):
        """
        _summary_

        Parameters:
        -----------
        engine: str
            MD engine
        
        Return:
        -------
        tuple:
            forward and backward network.dat files
        """
        self.forward = self.path + "network_fwd.dat"
        self.backward = self.path + "network_bwd.dat"
        with open(self.forward, "w") as network_file:
            for transformation, lomap_score in self.dictionary.items():
                self.n_windows = self.set_n_windows(lomap_score)
                lambda_array_bash = self.create_lambda_list_bash()
                network_file.write(f"{transformation[0]}, {transformation[1]}, {self.n_windows}, {lambda_array_bash}, {engine}\n")
        with open(self.backward, "w") as network_file:
            for transformation, lomap_score in self.dictionary.items():
                self.n_windows = self.set_n_windows(lomap_score)
                lambda_array_bash = self.create_lambda_list_bash()
                network_file.write(f"{transformation[1]}, {transformation[0]}, {self.n_windows}, {lambda_array_bash}, {engine}\n")
        return self.forward, self.backward
    
    def create_ligand_dat_file(self, afe_directory):
        """
        Create ligands.dat file

        Parameters:
        -----------
        afe_directory: 
            AFE working directory

        Return:
        -------
        ligands_dat: str
            ligands datafile
        """
        self.ligands_dat = afe_directory + "ligands.dat"
        with open(self.ligands_dat, "w") as ligands_file:
            writer = csv.writer(ligands_file)
            for ligand in self.names:
                writer.writerow([ligand])  
        return self.ligands_dat


    def dict_to_df(self):
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
        dictionary = self.create_dictionary()
        dataframe_from_dict = pd.DataFrame.from_dict(dictionary, "index")
        network_dataframe = dataframe_from_dict.reset_index().rename(columns={"index":"transformations",
                                                                            0:"score"})
        network_dataframe.index.name = "index"
        print("The LOMAP-generated perturbation network is:\n")
        print(network_dataframe)
        print("\n")
        return network_dataframe


    def adjust(self):
        """
        Ask user if they want to edit the LOMAP network

        Parameters:
        -----------
        Return:
        -------
        bool
            True: adjust network, False: do not adjust
        """
        self.adjust_network = False
        print("\nDo you want to edit this network? ([y]es/[n]o/[q]uit)")

        while True:
            edit = input("> ")
            if edit.lower().strip() == "yes" or edit.lower().strip() == "y":
                self.adjust_network = True
                break
            elif edit.lower().strip() == "no" or edit.lower().strip() == "n":
                self.adjust_network = False
                break
            elif edit.lower().strip() == "quit" or edit.lower().strip() == "q":
                print("Quitting.")
                sys.exit()
            elif edit == "":
                continue
            else: 
                print("Invalid option.")
                continue

        return self.adjust_network


    def print_options(self):
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
        

    def get_user_options(self):
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


    def delete_transformation(self):
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
            delete_indices = [int(self.options[i]) for i in range(1, len(self.options))]
            try:
                self.edited_dataframe = self.edited_dataframe.drop(delete_indices, axis=0, inplace=False)
                print("\nThe NEW network is:\n")
                print(self.edited_dataframe)
            except KeyError:
                print(f"Error: delete index should be between {self.first_index} and {self.last_index}")
                pass
        except ValueError:
            print("Delete index should be an integer")
            pass
        return self.edited_dataframe


    def edit_transformation(self):
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
            if len(self.options) != 3:
                raise ValueError
            index = int(self.options[1])
            score = np.float64(self.options[2])
            try:
                self.edited_dataframe.at[index, "score"] = score 
                print("\nThe NEW network is:\n")
                print(self.edited_dataframe)
            except KeyError:
                print(f"Error: delete index should be between {self.first_index} and {self.last_index}")
        except ValueError:
            print("edit existing LOMAP scores in format 'edit index score'")
            pass
        return self.edited_dataframe


    def add_transformation(self):
        """
        Edit transformation score as specified by user

        Parameters:
        -----------

        Return:
        -------
        pd.DataFrame
            updated network
        """
        try: 
            if len(self.options) != 4:
                raise ValueError
            transformation = tuple([self.options[1], self.options[2]])
            score = np.float64(self.options[3])
            added_row = pd.DataFrame([[transformation, score]], columns= ["transformations", "score"], index=[self.last_index+1])
            self.edited_dataframe = self.edited_dataframe.append(added_row, ignore_index=False)
            print("\nThe NEW network is:\n")
            print(self.edited_dataframe)
        except ValueError:
            print("Error: add transformations in format 'add ligand_1, ligand_2, score'")
        return self.edited_dataframe


    def edit_network(self):
        """
        Edit LOMAP network

        Parameters:
        -----------

        Return:
        -------
        dict
            adjusted network
        """
        self.first_index = self.dataframe.first_valid_index()
        self.last_index = self.dataframe.last_valid_index()
        self.adjust_network = self.adjust()
        self.edited_dataframe = self.dataframe.copy()
        while self.adjust_network:
            self.print_options()
            self.options = self.get_user_options()
            n_inputs = len(self.options)
            if not self.options:
                continue
            elif self.options[0] == "del":
                self.edited_dataframe = self.delete_transformation()
            elif self.options[0] == "edit":
                edited_dataframe = self.edit_transformation()
            elif self.options[0] == "add":
                edited_dataframe = self.add_transformation()
                last_index = edited_dataframe.last_valid_index()
            elif self.options[0] == "s":
                edited_dataframe.to_csv(self.path+"/meze_adjusted_network.csv", sep=":", header=None, index=False)
                break
            elif n_inputs == 1 and self.options[0].lower() != "q":
                print("Error: expected index")
                continue
            elif self.options[0].lower() == "q":
                print("Quitting.")
                sys.exit()
            else:
                print("Invalid option.")
                continue
        return dict(zip(edited_dataframe["transformations"], edited_dataframe["score"]))
       


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
    #TODO add a thing where you can input a network and then edit it here
    network_dict = create_network(ligands, ligand_names, ligand_path)
    edit_network(ligand_path, network_dict)


if __name__ == "__main__":
    main()