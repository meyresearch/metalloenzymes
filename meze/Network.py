import BioSimSpace as bss
bss.setVerbose(True)
import argparse
import functions
from argparse import RawTextHelpFormatter
import pandas as pd
import sys
from definitions import ADJUST_OPTIONS, PICOSECOND, ANGSTROM, KELVIN, ATM
import numpy as np
import shutil
import os
import logging
import Ligand
import csv
import Protein
import pathlib
import multiprocessing.pool
import time
from BioSimSpace import _Exceptions


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


def run_process(system, protocol, process, working_directory, configuration=None):
    """
    Run a Gromacs minimisation or equilibration process 
    Adapted from https://tinyurl.com/BSSligprep

    Parameters:
    -----------
    system: bss.System
        run system
    protocol: bss.Protocol 
        minimisation or equilibration
    process: name 
        process name for saving process output
    working_directory: str
        save output into this directory

    Return:
    -------
    system: bss.System
        equilibrated or minimised system
    """
    process = bss.Process.Gromacs(system, protocol, name=process, work_dir=working_directory, #exe="/usr/local/gromacs/bin/gmx_mpi"
                                  )
    config = process.getConfig()
    if configuration:
        for setting in configuration:
            key = setting.split()[0]
            try:
                index = [i for i, string in enumerate(config) if key in string][0]
                config[index] = setting
                process.setConfig(config)
            except IndexError:
                process.addToConfig(setting)
                config = process.getConfig()
    # process.setArg("-ntmpi", 1)
    process.start()
    process.wait()
    if process.isError():
        print(process.stdout())
        print(process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")
    system = process.getSystem()
    return system


def combine_unbound_ligands(system_1, system_2):
    """
    Take two unbound bss.Systems and combine the ligands' systems

    Parameters:
    -----------
    system_1: bss.System 
    system_2: bss.System

    Return:
    -------
    system_1: bss.System
        system with combined ligand topologies
    """
    ligand_1, ligand_2 = system_1.getMolecule(0), system_2.getMolecule(0)
    merged_ligands = merge_ligands(ligand_1, ligand_2)
    system_1.removeMolecules(ligand_1)
    system_1.addMolecules(merged_ligands)
    return system_1


def combine_bound_ligands(system_1, system_2):
    """
    Take two bound bss.Systems and combine the ligands' systems

    Parameters:
    -----------
    system_1: bss.System 
    system_2: bss.System

    Return:
    -------
    system_1: bss.System
        system with combined ligand topologies
    """
    ligand_1 = None
    protein = None
    n_residues = [molecule.nResidues() for molecule in system_1]
    n_atoms = [molecule.nAtoms() for molecule in system_1]
    for j, (n_residues, n_atoms) in enumerate(zip(n_residues[:20], n_atoms[:20])):
        if n_residues == 1 and n_atoms > 5:  
            ligand_1 = system_1.getMolecule(j)
        elif n_residues > 1:
            protein = system_1.getMolecules(j)
    ligand_2 = None
    n_residues = [molecule.nResidues() for molecule in system_2]
    n_atoms = [molecule.nAtoms() for molecule in system_2]   
    for j, (n_residues, n_atoms) in enumerate(zip(n_residues, n_atoms)):
        if n_residues == 1 and n_atoms > 5:
            ligand_2 = system_2.getMolecule(j)
    if ligand_1 and ligand_2 and protein:
        pass
    else:
        raise _Exceptions.AlignmentError("Could not extract ligands or protein from input systems.")
    merged_ligands = merge_ligands(ligand_1, ligand_2)
    system_1.removeMolecules(ligand_1)
    system_1.addMolecules(merged_ligands)
    return system_1


def merge_ligands(ligand_1, ligand_2):
    """
    Take two ligands and merge their topologies

    Parameters:
    -----------
    ligand_1: bss.Molecule
    ligand_2: bss.Molecule

    Return:
    -------
    merged_ligands: bss.Molecule
    """
    mapping = bss.Align.matchAtoms(ligand_1, ligand_2, complete_rings_only=True)
    inverse_mapping = {value:key for key, value in mapping.items()}
    aligned_ligand_2 = bss.Align.rmsdAlign(ligand_2, ligand_1, inverse_mapping)
    return bss.Align.merge(ligand_1, aligned_ligand_2, mapping, allow_ring_breaking=True, allow_ring_size_change=True)



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
    def __init__(self, workdir, ligand_path, group_name, protein_file, protein_path, water_model, ligand_ff, protein_ff, ligand_charge, 
                 engine, sampling_time, box_edges, box_shape, min_steps, short_nvt, nvt, npt, 
                 min_dt, min_tol, temperature=300, pressure=1, threshold=0.4, n_normal=11, n_difficult=17):
        """
        Class constructor
        """
        self.ligand_path = functions.path_exists(ligand_path) # ligand directory (change name?)
        self.ligand_forcefield = ligand_ff
        self.water_model = water_model
        self.input_files = self.get_files()
        self.ligand_charge = check_charge(ligand_charge)
        self.ligands = [Ligand.Ligand(file) for file in self.input_files]
        self.ligand_molecules = [ligand.get_ligand() for ligand in self.ligands]
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
        self.n_windows = []
        self.n_ligands = self.get_n_ligands()
        self.bound_ligands = [None] * self.n_ligands

        self.workding_directory = functions.path_exists(workdir)
        self.md_engine = engine
        self.md_time = sampling_time,
        self.box_shape = box_shape
        self.box_edges = box_edges
        self.min_steps = min_steps
        self.short_nvt = functions.convert_to_units(short_nvt, PICOSECOND)
        self.nvt = functions.convert_to_units(nvt, PICOSECOND)
        self.npt = functions.convert_to_units(npt, PICOSECOND)
        self.min_dt = min_dt
        self.min_tol = min_tol
        self.temperature = functions.convert_to_units(temperature, KELVIN)
        self.pressure = functions.convert_to_units(pressure, ATM)
        self.afe_directory = self.create_directory("/afe/")
        self.equilibration_directory = self.create_directory("/equilibration/")


    def create_directory(self, name):
        """
        Create AFE working directory in path.

        Parameters:
        -----------
        name: str
            name of new directory
        Return:
        -------
        directory: str
            full path to afe directory
        """
        try:
            directory = self.workding_directory + str(name)
            pathlib.Path(directory).mkdir(parents=False, exist_ok=False)
        except FileNotFoundError as e:
            print(f"Could not create directory {directory}. Pathlib raised error: {e}")
        except FileExistsError as e:
            print(f"Could not create directory {directory}. Pathlib raised error: {e}")
        return directory   


    def prepare_meze(self):
        """
        Prepare AFE calculations by creating network dictionary, ligand and protocol dat files.

        Parameters:
        -----------

        Return:
        -------
        self: Network
            (prepared) Network object
        """
        self.dictionary_fwd, self.dictionary_bwd = self.create_dictionary()
        self.forward, self.backward = self.create_network_files() # possibly do not need
        self.ligands_dat_file = self.create_ligand_dat_file()
        self.protocol_file = self.create_protocol_file()
        return self


    def solvation(self):
        """
        Use multiprocessing to solvate unbound and bound legs

        Parameters:
        -----------

        Return:
        -------
        self: Network
            (solvated) Network object
        """
        # For testing only
        # ligs = []
        # mols = []
        # for i in range(self.n_ligands):
        #     ligs.append(self.solvate_bound(i))
        #     mols.append(ligs[i].get_system())

        with multiprocessing.pool.Pool() as pool:
            self.ligands = pool.map(self.solvate_unbound, range(2))
            self.ligand_molecules = [ligand.get_system() for ligand in self.ligands]
            self.ligands = pool.map(self.solvate_bound, range(2))
            self.bound_ligands = [ligand.get_system() for ligand in self.ligands] 
        return self
    

    def equilibration(self):
        """
        Use multiprocessing to equilibrate unbound and bound legs

        Parameters:
        -----------

        Return:
        -------
        self: Network
            (solvated) Network object
        """
        # multiprocessing - takes up a lot of RAM!
        # start_equil = time.time()
        # with multiprocessing.pool.Pool() as pool:
        #     self.ligand_molecules = pool.map(self.heat_unbound, range(self.n_ligands))
        # print(f"\t Heating unbound took {time.time() - start_equil} s")
        # start_equil = time.time()
        # with multiprocessing.pool.Pool() as pool:
        #     self.bound_ligands = pool.map(self.heat_bound, range(self.n_ligands))
        # print(f"\t Heating bound took {time.time() - start_equil} s")

        start_equil = time.time()
        self.ligands = [self.heat_unbound(i) for i in range(2)]
        self.ligand_molecules = [ligand.get_system() for ligand in self.ligands]
        print(f"\t Heating unbound took {time.time() - start_equil} s")
        start_equil = time.time()
        self.ligands = [self.heat_bound(i) for i in range(2)]
        self.bound_ligands = [ligand.get_system() for ligand in self.ligands]
        print(f"\t Heating bound took {time.time() - start_equil} s")
        return self


    def afe_prep(self):

        for transformation, lomap_score in self.dictionary_fwd.items():
            ligand_1_name, ligand_2_name = transformation[0], transformation[1]
            # get_ligand_number = lambda ligand_name: ligand_name.split("_")[-1]
            # ligand_1_number = get_ligand_number(ligand_1_name)
            # ligand_2_number = get_ligand_number(ligand_2_name)   
            # self.ligands: list of Ligand objects
            unbound_ligand_1_system = self.match_name_to_system(self.ligands, ligand_1_name)
            unbound_ligand_2_system = self.match_name_to_system(self.ligands, ligand_2_name)
            bound_ligand_1_system = self.match_name_to_system(self.bound_ligands, ligand_1_name)
            bound_ligand_2_system = self.match_name_to_system(self.bound_ligands, ligand_2_name)

            unbound = combine_unbound_ligands(unbound_ligand_1_system, unbound_ligand_2_system)
            bound = combine_bound_ligands(bound_ligand_1_system, bound_ligand_2_system)

            free_energy_protocol = bss.Protocol.FreeEnergy()





            # print(f"lig {ligand_1_number} ----> lig {ligand_2_number}")

        # for transformation, lomap_score in self.dictionary_bwd.items():
        #     ligand_1, ligand_2 = transformation[0], transformation[1]
        #     get_ligand_number = lambda ligand_name: ligand_name.split("_")[-1]
        #     ligand_1_number = get_ligand_number(ligand_1)
        #     ligand_2_number = get_ligand_number(ligand_2)            

        #     print(f"lig {ligand_1_number} ----> lig {ligand_2_number}")

    def solvate_unbound(self, index):
        """
        Solvate unbound systems.

        Parameters:
        -----------
        index: list
            Ligand indices for sorting through Network.names and Network.ligands
        Return:
        -------
        solvated_ligand: Ligand
            (solvated) Ligand object whose file attribute is the prm7 and rst7 files 
        """ 
        ligand = self.ligands[index]
        names = self.names
        ligand_number = self.names[index].split("_")[-1]
        print(f"Solvating unbound ligand {ligand_number}")
        ligand_parameters = ligand.parameterise(self.ligand_forcefield, self.ligand_charge)
        unbound_box, unbound_box_angles = self.create_box(ligand_parameters)
        solvated_molecule = bss.Solvent.solvate(model=self.protein.water_model, 
                                              molecule=ligand_parameters, 
                                              box=unbound_box,
                                              angles=unbound_box_angles)
        ligand_savename = self.ligand_path + "ligand_" + ligand_number + "_solvated"
        solvated_files = bss.IO.saveMolecules(ligand_savename, solvated_molecule, ["PRM7", "RST7"])
        solvated_ligand = Ligand.Ligand(file=solvated_files, parameterised=True)
        return solvated_ligand
        

    def solvate_bound(self, index):
        """
        Solvate unbound systems.

        Parameters:
        -----------
        index: int
            Ligand indices for sorting through Network.names and Network.ligands
        Return:
        -----
        solvated_ligand: Ligand
            (solvated) Ligand object whose file attribute is the prm7 and rst7 files         
        """
        ligand = self.ligands[index]
        names = self.names
        ligand_number = self.names[index].split("_")[-1]
        ligand_parameters = ligand.parameterise(self.ligand_forcefield, self.ligand_charge)
        print(f"Solvating bound ligand {ligand_number}")        
        system_parameters = ligand_parameters + self.protein.get_prepared_protein()
        bound_box, bound_box_angles = self.create_box(system_parameters)
        solvated_molecules = bss.Solvent.solvate(model=self.protein.water_model,
                                                molecule=system_parameters,
                                                box=bound_box,
                                                angles=bound_box_angles)
        
        system_savename = self.protein_path + "system_" + ligand_number + "_solvated"
        solvated_files = bss.IO.saveMolecules(system_savename, solvated_molecules, ["PRM7", "RST7"])
        solvated_system = Ligand.Ligand(file=solvated_files, parameterised=True)
        return solvated_system


    def minimise(self, system, workdir):
        """
        Minimise the system using Gromacs

        Parameters:
        -----------
        system: bss.System
            system to be minimised
        working_directory: str
            current working dir

        Return:
        -------
        minimsed_system: bss.System
            minimised system
        """
        protocol = bss.Protocol.Minimisation(steps=self.min_steps)
        configuration = [f"emstep = {self.min_dt}", f"emtol = {self.min_tol}"]
        minimised_system = run_process(system, protocol, "min", workdir, configuration=configuration)
        return minimised_system        


                                                     # for start_t, end_t need to convert the default vals to Kelvin
    def equilibrate(self, system, name, workdir, time, start_t=300, end_t=300, temperature=None, pressure=None, configuration=None, restraints=None):
        """
        Run NVT or NPT equilibration

        Parameters:
        -----------
        system: bss.System
            system to be equilibrated

        Return:
        -------
        equilibrated_system: bss.System
            equilibrated system
        """
        protocol = bss.Protocol.Equilibration(runtime=time,
                                              temperature_start=start_t,
                                              temperature_end=end_t,
                                              temperature=temperature,
                                              pressure=pressure,
                                              restraint=restraints)
        equilibrated_system = run_process(system, protocol, name, workdir, configuration)
        return equilibrated_system


    def match_name_to_system(self, objects, name):
        """
        Take ligand name and return bss.System with that name

        Parameters:
        -----------
        name: str
            ligand name
        objects: list
            list of unbound/bound objects
        
        Return:
        -------
        bss.System: 
            system with the given name
        """
        for i in range(2):
            if objects[i].name == name:
                return objects[i].get_system()


    def heat_unbound(self, index):
        """
        Perform minimisation and NVT and NPT equilibrations on ligand

        Parameters:
        -----------
        index: int
            ligand index

        Return:
        -------
        equilibrated_ligand: bss.System
            equilibrated ligand object
        """
        ligand_number = self.names[index].split("_")[-1]
        directory = functions.mkdir(self.equilibration_directory + f"/unbound/ligand_{ligand_number}/")
        solvated_ligand = self.ligand_molecules[index]
        print(f"Equilibrating unbound ligand {ligand_number}")
        directories = lambda step: functions.mkdir(directory+step)
        min_directory = directories("min")
        r_nvt_directory = directories("r_nvt")
        nvt_directory = directories("nvt")
        r_npt_directory = directories("r_npt")
        npt_directory = directories("npt")

        minimised_ligand = self.minimise(system=solvated_ligand, workdir=min_directory)
        start_temp = functions.convert_to_units(0, KELVIN)
        restrained_nvt = self.equilibrate(system=minimised_ligand,
                                          name="r_nvt",
                                          workdir=r_nvt_directory,
                                          time=self.short_nvt,
                                          start_t=start_temp, end_t=self.temperature,
                                          configuration=["dt = 0.0005"], # need to be able to change
                                          restraints="all")
        nvt = self.equilibrate(system=restrained_nvt,
                               name="nvt",
                               workdir=nvt_directory,
                               time=self.nvt,
                               temperature=self.temperature)
        restrained_npt = self.equilibrate(system=nvt,
                                          name="r_npt",
                                          workdir=r_npt_directory,
                                          time=self.npt,
                                          pressure=self.pressure,
                                          temperature=self.temperature,
                                          restraints="heavy")
        equilibrated_molecule = self.equilibrate(system=restrained_npt,
                                               name="npt",
                                               workdir=npt_directory,
                                               time=self.npt,
                                               pressure=self.pressure,
                                               temperature=self.temperature)
        unbound_savename = npt_directory + f"/ligand_{ligand_number}"
        equilibrated_files = bss.IO.saveMolecules(filebase=unbound_savename, system=equilibrated_molecule, fileformat=["PRM7", "RST7"])        
        equilibrated_ligand = Ligand.Ligand(file=equilibrated_files, parameterised=True)
        return equilibrated_ligand
    

    def heat_bound(self, index):
        """
        Perform minimisation and NVT and NPT equilibrations on bound ligand 

        Parameters:
        -----------
        index: int
            ligand index

        Return:
        -------
        equilibrated_system: bss.System
            equilibrated system object
        """        
        ligand_number = self.names[index].split("_")[-1]
        solvated_system = self.bound_ligands[index]
        directory = functions.mkdir(self.equilibration_directory+f"/bound/ligand_{ligand_number}/")
        directories = lambda step: functions.mkdir(directory+step)
        min_dir = directories("min")
        r_nvt_dir = directories("r_nvt")
        bb_r_nvt_dir = directories("bb_r_nvt")
        nvt_dir = directories("nvt")
        r_npt_dir = directories("r_npt")
        npt_dir = directories("npt")     
        start_temp = functions.convert_to_units(0, KELVIN)
        print(f"Equilibrating bound ligand {ligand_number}")
        minimised_system = self.minimise(system=solvated_system, workdir=min_dir)
        restrained_nvt = self.equilibrate(system=minimised_system,
                                          workdir=r_nvt_dir,
                                          name="r_nvt",
                                          time=self.short_nvt,
                                          start_t=start_temp, end_t=self.temperature,
                                          restraints="all",
                                          configuration=["dt = 0.0005"]) # need to be able to change
        backbone_restrained_nvt = self.equilibrate(system=restrained_nvt,
                                                   name="bb_r_nvt",
                                                   workdir=bb_r_nvt_dir,
                                                   time=self.nvt,
                                                   temperature=self.temperature,
                                                   restraints="backbone")
        nvt = self.equilibrate(system=backbone_restrained_nvt,
                               name="nvt",
                               workdir=nvt_dir,
                               time=self.nvt,
                               temperature=self.temperature)
        restrained_npt = self.equilibrate(system=nvt,
                                          name="r_npt",
                                          workdir=r_npt_dir,
                                          time=self.npt,
                                          pressure=self.pressure,
                                          temperature=self.temperature,
                                          restraints="heavy")
        equilibrated_protein = self.equilibrate(system=restrained_npt,
                                               name="npt",
                                               workdir=npt_dir,
                                               time=self.npt,
                                               pressure=self.pressure,
                                               temperature=self.temperature)
        bound_savename = npt_dir + f"/system_{ligand_number}"
        equilibrated_files = bss.IO.saveMolecules(filebase=bound_savename, system=equilibrated_protein, fileformat=["PRM7", "RST7"])     
        equilibrated_system = Ligand.Ligand(file=equilibrated_files, parameterised=True)
        return equilibrated_system


    def create_box(self, molecule):
        """
        Create a bss.Box object for solvation.

        Parameters:
        -----------
        molecule: 
            bss.Molecule: usually either a protein or a ligand

        Return:
        -------
        tuple: 
            bss.Box and angles
        """

        box_min, box_max = molecule.getAxisAlignedBoundingBox()
        box_size = [y - x for x, y in zip(box_min, box_max)]
        box_area = [x + int(self.box_edges) * ANGSTROM for x in box_size]
        self.box, self.box_angles = None, None
        if self.box_shape == "cubic":
            self.box, self.box_angles = bss.Box.cubic(max(box_area))
        elif self.box_shape == "rhombicDodecahedronHexagon":
            self.box, self.box_angles = bss.Box.rhombicDodecahedronHexagon(max(box_area))
        elif self.box_shape == "rhombicDodecahedronSquare":
            self.box, self.box_angles = bss.Box.rhombicDodecahedronSquare(max(box_area))
        elif self.box_shape == "truncatedOctahedron":
            self.box, self.box_angles = bss.Box.truncatedOctahedron(max(box_area))
        else:
            print(f"Box shape {self.box_shape} not supported.")
        return self.box, self.box_angles
    

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
        network_dict_fwd, network_dict_bwd = {}, {}
        named_transformations_fwd = [(self.names[transformation[0]], self.names[transformation[1]]) for transformation in transformations]
        named_transformations_bwd = [(self.names[transformation[1]], self.names[transformation[0]]) for transformation in transformations]

        for transformation, score in zip(named_transformations_fwd, lomap_scores):
            network_dict_fwd[transformation] = score
        for transformation, score in zip(named_transformations_bwd, lomap_scores):
            network_dict_bwd[transformation] = score

        with open(self.ligand_path+f"/meze_network.csv", "w") as lomap_out:
            for key, value in network_dict_fwd.items():
                lomap_out.write(f"{key}: {value}\n")
        return network_dict_fwd, network_dict_bwd
    

    def get_files(self):
        """
        Read in all .sdf and .mol2 files in the provided path.

        Return:
        -------
        list
            list of ligand filenames
        """
        # Adapted from dbmol.py:
        ligand_files = functions.read_files(f"{self.ligand_path}/*.sdf")
        ligand_files += functions.read_files(f"{self.ligand_path}/*.mol2")        
        if len(ligand_files) < 2:
            raise IOError(f"Path {self.ligand_path} must contain at least two sdf or mol2 files.")
        return ligand_files


    def get_n_ligands(self):
        """
        Get number of ligands

        Return:
        -------
        len(self.files): int
            number of ligands
        """
        return(len(self.input_files))


    def create_new_lomap_directory(self):
        """
        Create a lomap directory in ligand path 

        Return:
        -------
        str
            new LOMAP working directory 
        """
        new_lomap_directory = self.ligand_path + "/lomap/"
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
        lomap_work_directory = self.ligand_path
        lomap_names = ["/images/", "/inputs/", "/outputs/"]
        lomap_directories = [self.ligand_path + name for name in lomap_names]

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
            number of lambda windows for given transformation
        """
        if lomap_score == None or lomap_score < float(self.threshold):
            n_windows = self.n_difficult
        else:
            n_windows = self.n_normal
        return n_windows


    def create_lambda_list_bash(self, n_windows): 
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


    def create_network_files(self):
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
        forward = self.workding_directory + "network_fwd.dat"
        backward = self.workding_directory + "network_bwd.dat"
        with open(forward, "w") as network_file:
            for transformation, lomap_score in self.dictionary_fwd.items():
                n_windows = self.set_n_windows(lomap_score)
                self.n_windows.append(n_windows)
                lambda_array_bash = self.create_lambda_list_bash(n_windows)
                network_file.write(f"{transformation[0]}, {transformation[1]}, {n_windows}, {lambda_array_bash}, {self.md_engine}\n")
        with open(backward, "w") as network_file:
            for transformation, lomap_score in self.dictionary_bwd.items():
                n_windows = self.set_n_windows(lomap_score)
                lambda_array_bash = self.create_lambda_list_bash(n_windows)
                network_file.write(f"{transformation[0]}, {transformation[1]}, {n_windows}, {lambda_array_bash}, {self.md_engine}\n")
        return forward, backward
    

    def create_ligand_dat_file(self):
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
        ligands_dat = self.afe_directory + "ligands.dat"
        with open(ligands_dat, "w") as ligands_file:
            writer = csv.writer(ligands_file)
            for ligand in self.names:
                writer.writerow([ligand])  
        return ligands_dat


    def create_protocol_file(self):
        """
        Create protocol.dat file for AFE runs

        Parameters:
        -----------
        Network: Network
            Network class object

        Return:
        -------
        protocol_file: str
            protocol datafile
        """
        protocol = [f"ligand forcefield = {self.ligand_forcefield}", 
                    f"protein forcefield = {self.protein_forcefield}", 
                    f"solvent = {self.water_model}", 
                    f"box edges = {self.box_edges}*angstrom", 
                    f"box shape = {self.box_shape}", 
                    f"protocol = default",
                    f"sampling = {self.md_time}*ns",
                    f"engine = {self.md_engine}"]
        protocol_file = self.workding_directory + "/protocol.dat"

        with open(protocol_file, "w") as file:
            writer = csv.writer(file)
            for protocol_line in protocol:
                writer.writerow([protocol_line])
        return protocol_file
    

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
                edited_dataframe.to_csv(self.ligand_path+"/meze_adjusted_network.csv", sep=":", header=None, index=False)
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