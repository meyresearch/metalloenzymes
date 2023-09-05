import BioSimSpace as bss
bss.setVerbose(True)
import functions
import pandas as pd
import sys
from definitions import ADJUST_OPTIONS, PICOSECOND, NANOSECOND, ANGSTROM, KELVIN, ATM
import numpy as np
import shutil
import os
import logging
import Ligand
import csv
import Protein
import pathlib
import multiprocessing.pool
import tqdm
import istarmap
import shutil
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


def combine_unbound_ligands(system_a, system_b):
    """
    Take two unbound bss.Systems and combine the ligands' systems

    Parameters:
    -----------
    a: Ligand() 
    b: Ligand()

    Return:
    -------
    system_1: bss.System
        system with combined ligand topologies
    """
    # system_a, system_b = a.get_system(), b.get_system()
    ligand_1, ligand_2 = system_a.getMolecule(0), system_b.getMolecule(0)
    merged_ligands = merge_ligands(ligand_1, ligand_2)
    system_a.removeMolecules(ligand_1)
    system_a.addMolecules(merged_ligands)
    return system_a


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
            protein = system_1.getMolecule(j)
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


def create_lambda_windows(n_windows):
    """
    Create evenly spaced list of lambda windows 

    Parameters:
    -----------
    n_windows: int
        number of windows

    Return:
    -------
    list: 
        list of lambdas
    """
    return list(np.linspace(0, 1, int(n_windows)))


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
    lambda_list_numpy = create_lambda_windows(n_windows)
    lambda_list = [format(item, ".4f") for item in lambda_list_numpy]
    return " ".join(lambda_list)


def create_minimisation_configs(files):
    """
    Open FreeEnergy configuration file and convert it into a minimisation configuration

    Parameters:
    -----------
    files: list
        list of configuration files in each lambda minimisation window

    Return:
    -------
    """
    minimisation_config = ["minimise = True\n", "minimise maximum iterations = 10000\n"]
    for i in range(len(files)):
        with open(files[i], "r") as f:
            old_config = f.readlines()
        with open(files[i], "w") as f:
            replaced_config = [setting.replace("ncycles = 5\n", "ncycles = 1").replace("nmoves = 200000", "nmoves = 50000") for setting in old_config]
            for setting in minimisation_config:
                replaced_config.append(setting)
            f.writelines(replaced_config)


def construct_relative_afe(system, protocol, engine, workdir):
    """
    Wrap the BSS.FreeEnergy.Relative so it is easier to run with multiprocessing

    Parameters:
    -----------
    system: bss.System
    protocol: bss.Protocol
    engine: str
    workdir: str
        
    Return:
    -------
    bss.FreeEnergy.Relative: 
        class constructor to setup free energy calculations
    """
    return bss.FreeEnergy.Relative(system=system, protocol=protocol, engine=engine, work_dir=workdir, setup_only=True)


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
                 min_dt, min_tol, repeats, temperature=300, pressure=1, threshold=0.4, n_normal=11, n_difficult=17):
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
        self.lambdas = []
        self.n_ligands = self.get_n_ligands()
        self.bound_ligands = [None] * self.n_ligands
        self.bound_ligand_molecules = [None] * self.n_ligands

        self.workding_directory = functions.path_exists(workdir)

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
        
        self.md_engine = engine
        self.md_time = functions.convert_to_units(sampling_time, NANOSECOND)
        self.n_repeats = repeats
        self.log_directory = self.create_directory("/logs/")
        self.afe_input_directory = self.create_directory("/afe/")
        self.equilibration_directory = self.create_directory("/equilibration/")
        self.output_directories = self.create_output_directories()


    def create_directory(self, name, create_parents=False):
        """
        Create AFE working directory in path.

        Parameters:
        -----------xantham gum
        name: str
            name of new directory
        Return:
        -------
        directory: str
            full path to afe directory
        """
        try:
            directory = self.workding_directory + str(name)
            pathlib.Path(directory).mkdir(parents=create_parents, exist_ok=False)
        except FileNotFoundError as e:
            print(f"Could not create directory {directory}. Pathlib raised error: {e}")
        except FileExistsError as e:
            print(f"Could not create directory {directory}. Pathlib raised error: {e}")
        return directory   


    def create_output_directories(self):
        """
        Create bound and unbound output directories for each repeat AFE run

        Parameters:
        -----------

        Return:
        -------
        output, unbound, bound: tuple
            parent directories list, unbound dirs, bound dirs
        """
        output_directories = []
        for i in range(1, self.n_repeats + 1, 1):
            output_directories.append(self.create_directory(f"/outputs/{self.md_engine}_{i}/", create_parents=True))
        return output_directories


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
        self.transformations = self.set_transformations()
        self.n_transformations = len(self.transformations)
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
        self.ligands, self.bound_ligands = [], []
        with multiprocessing.pool.Pool() as pool:
            for result in tqdm.tqdm(pool.imap(self.solvate_unbound, range(self.n_ligands)), desc="Solvate unbound\n", total=self.n_ligands):
                self.ligands.append(result)
            self.ligand_molecules = [ligand.get_system() for ligand in self.ligands]
            for result in tqdm.tqdm(pool.imap(self.solvate_bound, range(self.n_ligands)), desc="Solvate bound\n", total=self.n_ligands):
                self.bound_ligands.append(result)
            self.bound_ligand_molecules = [ligand.get_system() for ligand in self.ligands] 
        return self


    def get_equilibrated(self):
        """
        Read in equilibrated files and update Network class object

        Parameters:
        -----------

        Return:
        -------
        self: Network
            (equilibrated) Network object
        """
        print("\n")
        unbound_paths = functions.read_files(self.equilibration_directory+"/unbound/ligand_*/npt/")
        paths = list(map(lambda x: x + "ligand_*", unbound_paths))
        unbound_equilibrated = []
        with multiprocessing.pool.Pool() as pool:
            for result in tqdm.tqdm(pool.imap(functions.read_files, paths), desc="Equilbrated unbound", total=len(paths)):
                unbound_equilibrated.append(result)
        print("\n")
        self.ligand_molecules = []
        with multiprocessing.pool.Pool() as pool:
            for result in tqdm.tqdm(pool.imap(bss.IO.readMolecules, unbound_equilibrated), desc="Unbound objects", total=len(unbound_equilibrated)):
                self.ligand_molecules.append(result)
        print("\n")
        bound_paths = functions.read_files(self.equilibration_directory+"/bound/ligand_*/npt/")
        paths = list(map(lambda x: x + "system_*", bound_paths))
        bound_equilibrated = []
        with multiprocessing.pool.Pool() as pool:
            for result in tqdm.tqdm(pool.imap(functions.read_files, paths), desc="Equilibrated bound", total=len(paths)):
                bound_equilibrated.append(result)
        print("\n")
        self.bound_ligand_molecules = []
        with multiprocessing.pool.Pool() as pool:
            for result in tqdm.tqdm(pool.imap(bss.IO.readMolecules, bound_equilibrated), desc="Bound objects", total=len(bound_equilibrated)):
                self.bound_ligand_molecules.append(result)
        print("\n")
        return self


    def afe_prep(self):
        """
        Prepare minimisation and free energy lambda windows directory tree

        Parameters:
        -----------

        Return:
        -------
        """
        # forward
        columns_to_list = lambda column: self.transformations[column].tolist()

        ligands_a, ligands_b = columns_to_list("ligand_a"), columns_to_list("ligand_b")

        indices_a, indices_b = columns_to_list("index_a"), columns_to_list("index_b")
        lambdas = self.transformations["lambdas"].to_numpy()

        arguments = [(self.ligand_molecules[indices_a[i]], self.ligand_molecules[indices_b[i]]) for i in range(self.n_transformations)]
        unbound_systems = []
        with multiprocessing.pool.Pool() as pool:
            for result in tqdm.tqdm(pool.istarmap(combine_unbound_ligands, arguments), desc="Merging unbound", total=len(arguments)):
                unbound_systems.append(result)
        print("\n")
        arguments = [(self.bound_ligand_molecules[indices_a[i]], self.bound_ligand_molecules[indices_b[i]]) for i in range(self.n_transformations)]
        bound_systems = []
        with multiprocessing.pool.Pool() as pool:
            for result in tqdm.tqdm(pool.istarmap(combine_bound_ligands, arguments), desc="Merging bound", total=len(arguments)):
                bound_systems.append(result)
        print("\n")
        free_energy_protocols = [bss.Protocol.FreeEnergy(lam_vals=lambdas[i], runtime=self.md_time) for i in tqdm.tqdm(range(self.n_transformations), desc="Free energy protocols")]
        print("\n")
        first_run_directory = self.output_directories[0]
        transformation_directories = [first_run_directory + f"/{ligands_a[i]}~{ligands_b[i]}/" for i in range(self.n_transformations)]
        bound_directories = [directory + "/bound/" for directory in transformation_directories]
        unbound_directories = [directory + "/unbound/" for directory in transformation_directories]

        _ = [bss.FreeEnergy.Relative(system=unbound_systems[i], protocol=free_energy_protocols[i], engine=self.md_engine, work_dir=unbound_directories[i]) for i in tqdm.tqdm(range(self.n_transformations), desc="Unbound AFE")]
        print("\n")
        _ = [bss.FreeEnergy.Relative(system=bound_systems[i], protocol=free_energy_protocols[i], engine=self.md_engine, work_dir=bound_directories[i]) for i in tqdm.tqdm(range(self.n_transformations), desc="Bound AFE")]
        print("\n")

        _ = [shutil.copytree(unbound_directories[i], self.output_directories[j], dirs_exist_ok=True) for i in tqdm.tqdm(range(len(unbound_directories)), desc="Copy unbound") for j in range(self.n_repeats)]
        print("\n")
        _ = [shutil.copytree(bound_directories[i], self.output_directories[j], dirs_exist_ok=True) for i in tqdm.tqdm(range(len(bound_directories)), desc="Copy bound") for j in range(self.n_repeats)]
        print("\n")
        bound_lambda_minimisation_directories = [functions.mkdir(directory + "/minimisation/") for directory in bound_directories]
        unbound_lambda_minimisation_directories = [functions.mkdir(directory + "/minimisation/") for directory in unbound_directories]


        print("\n")
        _ = [shutil.copytree(unbound_directories[i], unbound_lambda_minimisation_directories[i], dirs_exist_ok=True) for i in tqdm.tqdm(range(len(unbound_lambda_minimisation_directories)), desc="Unbound minimisation")]
        print("\n")
        _ = [shutil.copytree(bound_directories[i], bound_lambda_minimisation_directories[i], dirs_exist_ok=True) for i in tqdm.tqdm(range(len(bound_lambda_minimisation_directories)), desc="Bound minimisation")]
        
        bound_configurations = list(map(lambda x: x + "/*/*.cfg", bound_lambda_minimisation_directories))
        unbound_configurations = list(map(lambda x: x + "/*/*.cfg", unbound_lambda_minimisation_directories))

        bound_configuration_files = [functions.read_files(file) for file in bound_configurations]
        unbound_configuration_files = [functions.read_files(file) for file in unbound_configurations]

        _ = [create_minimisation_configs(config_files) for config_files in bound_configuration_files]
        _ = [create_minimisation_configs(config_files) for config_files in unbound_configuration_files]

        print("\n")
        _ = [shutil.copytree(unbound_lambda_minimisation_directories[i], self.output_directories[j], dirs_exist_ok=True) for i in tqdm.tqdm(range(len(unbound_lambda_minimisation_directories)), desc="Copy unbound") for j in range(1, self.n_repeats)]
        print("\n")
        _ = [shutil.copytree(bound_lambda_minimisation_directories[i], self.output_directories[j], dirs_exist_ok=True) for i in tqdm.tqdm(range(len(bound_lambda_minimisation_directories)), desc="Copy bound") for j in range(1, self.n_repeats)]
        

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
    

    def set_transformations(self):
        """
        Create a transformation network with LOMAP

        Return:
        -------
        transformations: pd.DataFrame
            forward and backward transformations with their associated lomap scores, number of windows and lambda list
        """
        lomap_work_directory = self.check_lomap_directory()
        transformations, lomap_scores = bss.Align.generateNetwork(self.ligand_molecules, plot_network=True, names=self.names, work_dir=lomap_work_directory)        
        start_ligand = [self.names[transformation[0]] for transformation in transformations]
        end_ligand = [self.names[transformation[1]] for transformation in transformations] 
        start_indices = [transformation[0] for transformation in transformations]
        end_indices = [transformation[1] for transformation in transformations]
        dataframe = pd.DataFrame()
        dataframe["ligand_a"] = start_ligand
        dataframe["index_a"] = start_indices
        dataframe["ligand_b"] = end_ligand
        dataframe["index_b"] = end_indices
        dataframe["score"] = lomap_scores
        for score in lomap_scores:
            n_windows = self.set_n_windows(score)
            self.n_windows.append(n_windows)
            lambda_windows = create_lambda_windows(n_windows)
            self.lambdas.append(lambda_windows)
        dataframe["n_windows"] = self.n_windows
        dataframe["lambdas"] = self.lambdas
        dataframe.to_csv(self.afe_input_directory+f"/meze_network.csv", sep="\t")
        return dataframe
    

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
        afe_input_directory: 
            dir for afe input files

        Return:
        -------
        ligands_dat: str
            ligands datafile
        """
        ligands_dat = self.afe_input_directory + "ligands.dat"
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
        protocol_file = self.afe_input_directory + "/protocol.dat"

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
        dictionary = self.set_transformations()
        dataframe_from_dict = pd.DataFrame.from_dict(dictionary, "index")
        network_dataframe = dataframe_from_dict.reset_index().rename(columns={"index":"transformations", 0:"score"})
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
   pass

if __name__ == "__main__":
    main()