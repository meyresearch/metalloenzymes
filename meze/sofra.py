import re
import BioSimSpace as bss
import functions
import Protein
import Ligand
import pandas as pd
import sys
from definitions import ADJUST_OPTIONS, PICOSECOND, NANOSECOND, ANGSTROM, KELVIN, ATM
import numpy as np
import shutil
import os
import logging
import csv
import pathlib
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



def combine_unbound_ligands(system_a, system_b, flexible_align=False):
    """
    Take two unbound bss.Systems and combine the ligands' systems

    Parameters:
    -----------
    system_a: bss.Molecule
    system_b: bss.Molecule

    Return:
    -------
    system_1: bss.System
        system with combined ligand topologies
    """
    ligand_1, ligand_2 = system_a.getMolecule(0), system_b.getMolecule(0)
    merged_ligands = merge_ligands(ligand_1, ligand_2, flexible_align=flexible_align)
    system_a.removeMolecules(ligand_1)
    system_a.addMolecules(merged_ligands)
    return system_a


def merge_ligands(ligand_1, ligand_2, flexible_align=False):
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
    if flexible_align:
        aligned_ligand_2 = bss.Align.flexAlign(ligand_2, ligand_1, inverse_mapping)
    else:
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
    str: 
        string of lambdas as strings formatted to 4 decimal places
    """
    return " ".join([format(item, ".4f") for item in np.linspace(0, 1, int(n_windows))])


def fix_afe_configurations(files):
    """
    Open FreeEnergy configuration file and remove the gpu id specification line
    Fixes https://github.com/OpenBioSim/biosimspace/issues/180 

    Parameters:
    -----------
    files: list
        list of configuration files in each lambda minimisation window
    
    Return:
    -------
    """
    for i in range(len(files)):
        with open(files[i], "r") as f:
            old_config = f.readlines()
        for line in old_config:
            if "gpu" in line:
                idx = old_config.index(line)
                del old_config[idx]
        replaced_config = old_config
        with open(files[i], "w") as f:
            f.writelines(replaced_config)


def create_minimisation_configurations(files, min_cycles=1, min_moves=50000):
    """
    Open FreeEnergy configuration file and convert it into a minimisation configuration

    Parameters:
    -----------
    files: list
        list of configuration files in each lambda minimisation window
    min_cycles: int
        number of cycles for the SOMD minimisation
    min_moves: int
        number of moves for the SOMD minimisation
    
    Return:
    -------
    """

    for i in range(len(files)):
        with open(files[i], "r") as f:
            old_config = f.readlines()
        for line in old_config:
            if "ncycles" in line and "ncycles_per_snap" not in line:
                idx = old_config.index(line)
                old_config[idx] = f"ncycles = {min_cycles}\n"
            elif "nmoves" in line:
                idx = old_config.index(line)
                old_config[idx] = f"nmoves = {min_moves}\n"
            elif "buffered coordinates frequency" in line:
                idx = old_config.index(line)
                frames = min_moves // 100
                old_config[idx] = f"buffered coordinates frequency = {frames}\n"
            elif "minimise" in line and "minimise maximum iterations" not in line:
                idx = old_config.index(line)
                old_config[idx] = "minimise = True\n"
        replaced_config = old_config

        with open(files[i], "w") as f:
            f.writelines(replaced_config)


class Sofra(object):
    """
    Network class object

    Parameters:
    -----------
    object: 
        _description_

    Return:
    -------
    """
    def __init__(self, protein_file, workdir=os.getcwd(), prepared=False, is_md=False, md_input_directory=None,
                 afe_input_path=os.getcwd()+"/afe/", equilibration_path=os.getcwd()+"/equilibration/", outputs=os.getcwd()+"/outputs/", log_directory=os.getcwd()+"/logs/",
                 ligand_path=os.getcwd()+"/inputs/ligands/", ligand_charge=0, ligand_ff="gaff2", 
                 group_name=None, protein_path=os.getcwd()+"/inputs/protein/", water_model="tip3p", protein_ff="ff14SB", 
                 engine="SOMD", sampling_time=4, box_edges=20, box_shape="cubic", min_steps=5000, short_nvt=5, nvt=50, npt=200, 
                 min_dt=0.01, min_tol=1000, repeats=3, temperature=300, pressure=1, threshold=0.4, n_normal=11, n_difficult=17,
                 cutoff_scheme="rf", solvation_method="gromacs", solvent_closeness=1.0, only_save_end_states=False):
        """
        Class constructor
        """
        self.working_directory = functions.path_exists(workdir)
        self.only_save_end_states = only_save_end_states
        self.md_engine = engine
        self.md_time = functions.convert_to_units(sampling_time, NANOSECOND)
        self.is_md = is_md
        if not self.is_md:
            self.n_repeats = repeats
        else:
            self.n_repeats = 0

        self.prepared = prepared 
        if self.prepared: 
            self.prepared_protein = functions.get_files(protein_file + ".*")
            self.protein_file = self.prepared_protein
            if not is_md:
                self.afe_input_directory = functions.path_exists(afe_input_path)
                self.output_directories = functions.get_files(outputs + f"/{engine}_*/")
            else: 
                self.md_input_directory = functions.path_exists(md_input_directory)
            self.equilibration_directory = functions.path_exists(equilibration_path)  
            self.outputs = functions.path_exists(outputs) 
            self.log_directory = functions.path_exists(log_directory)
            self.plots = functions.get_files(f"{self.outputs}/plots/")
        else:
            self.protein_file = protein_file
            if not self.is_md:
                self.afe_input_directory = self.create_directory(afe_input_path)
                self.output_directories = self.create_output_directories()
            else:
                self.md_input_directory = self.create_directory(md_input_directory)
            self.equilibration_directory = self.create_directory(equilibration_path)
            self.outputs = self.create_directory(outputs)
            self.plots = self.create_directory(f"{self.outputs}/plots/")
            self.log_directory = self.create_directory(log_directory)
        self.solvation_method = solvation_method
        self.solvent_closeness = functions.check_positive(solvent_closeness)

        self.ligand_path = functions.path_exists(ligand_path)
        self.ligand_forcefield = ligand_ff
        self.water_model = water_model.lower()
        self.input_files = self.get_files()
        self.ligand_charge = check_charge(ligand_charge)
        self.ligands = [Ligand.Ligand(file) for file in self.input_files]
        self.ligand_molecules = [ligand.get_ligand() for ligand in self.ligands]
        self.names = [ligand.get_name() for ligand in self.ligands]

        self.protein_forcefield = protein_ff
        self.protein_path = functions.path_exists(protein_path)
        self.group_name = self.get_name(group_name)
        self.protein = Protein.Protein(name=self.group_name,
                                       protein_file=self.protein_file,
                                       path=self.protein_path,
                                       forcefield=self.protein_forcefield,
                                       water_model=self.water_model,
                                       parameterised=self.prepared)
        
        self.cutoff_scheme = cutoff_scheme.lower()
        if not is_md:
            self.threshold = threshold
            self.n_normal = n_normal
            self.n_difficult = n_difficult
            self.n_windows = []
            self.lambdas = []
        self.short_nvt = functions.convert_to_units(short_nvt, PICOSECOND)
        self.nvt = functions.convert_to_units(nvt, PICOSECOND)
        self.npt = functions.convert_to_units(npt, PICOSECOND)
        self.n_ligands = self.get_n_ligands()
        self.bound_ligands = [None] * self.n_ligands
        self.bound_ligand_molecules = [None] * self.n_ligands

        self.box_shape = box_shape
        self.box_edges = box_edges
        self.min_steps = min_steps
        self.min_dt = min_dt
        self.min_tol = min_tol
        self.temperature = functions.convert_to_units(temperature, KELVIN)
        self.pressure = functions.convert_to_units(pressure, ATM)
        

    def get_ligand_by_name(self, name):
        """
        Extract a single ligand from the list of network ligands using a name

        Parameters:
        -----------
        name: str
            name of the ligand to be extracted

        Return:
        -------
        Ligand
            Ligand object whose name matches given name
        """
        return [ligand for ligand in self.ligands if ligand.name == name][0]



    def create_directory(self, name, create_parents=False):
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
            directory = str(name)
            pathlib.Path(directory).mkdir(parents=create_parents, exist_ok=False)
        except FileNotFoundError as e:
            print(f"Could not create directory {directory}. Pathlib raised error: {e}")
        except FileExistsError as e:
            print(f"Could not create directory {directory}. Pathlib raised error: {e}")
        # https://chat.openai.com/share/6bf91255-5a4c-4a4f-b7d8-b489e784397f 
        return re.sub(r"/+", "/", directory)


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
            output_directories.append(self.create_directory(f"{self.outputs}/{self.md_engine}_{i}/"))
        return output_directories


    def prepare_network(self):
        """
        Prepare AFE calculations by creating network dictionary, ligand and protocol dat files.

        Parameters:
        -----------

        Return:
        -------
        self: Network
            (prepared) Network object
        """
        self.protein_water_complex = self.protein.create_complex()
        self.prepared_protein = self.protein.tleap(self.protein_water_complex)
        self.log_directory = self.create_directory(f"{self.working_directory}/logs/")
        self.transformations, self.network_file = self.set_transformations()
        self.n_transformations = len(self.transformations)
        self.ligands_dat_file = self.create_ligand_dat_file()
        self.protocol_file = self.create_protocol_file()
        self.prepared = True
        return self


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
        elif self.box_shape == "truncatedOctahedron" or "octahedron" in self.box_shape:
            self.box, self.box_angles = bss.Box.truncatedOctahedron(max(box_area))
        else:
            print(f"Box shape {self.box_shape} not supported.")
        return self.box, self.box_angles


    def get_equilibrated(self, ligand_a, ligand_b):
        """
        Read in equilibrated ligand files and update Network class object

        Parameters:
        -----------
        ligand_a: str
            name of ligand a 
        ligand_b: str
            name of ligand b
        Return:
        -------
        self: Network
            (equilibrated) Network object
        """
        print("\n")
        print(f"Getting equilibrated ligands {ligand_a} and {ligand_b}")
        unbound_files = [functions.get_files(self.equilibration_directory+f"/unbound/{ligand}/npt/{ligand}.*") for ligand in [ligand_a, ligand_b]]
        self.ligands = [Ligand.Ligand(files, parameterised=True) for files in unbound_files]
        self.names = [ligand.get_name() for ligand in self.ligands]
        self.ligand_molecules = [bss.IO.readMolecules(files) for files in unbound_files]
        bound_files = [functions.get_files(self.equilibration_directory+f"/bound/{ligand}/npt/bound_{ligand}.*") for ligand in [ligand_a, ligand_b]]
        self.bound_ligands = [Ligand.Ligand(files, parameterised=True) for files in bound_files]
        self.bound_ligand_molecules = [bss.IO.readMolecules(files) for files in bound_files]
        return self


    def combine_bound_ligands(self, flexible_align=False):
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
        system_1 = self.bound_ligand_molecules[0]
        system_2 = self.bound_ligand_molecules[1]
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
        merged_ligands = merge_ligands(ligand_1, ligand_2, flexible_align=flexible_align)
        system_1.removeMolecules(ligand_1)
        system_1.addMolecules(merged_ligands)
        return system_1


    def prepare_afe(self, ligand_a_name, ligand_b_name, extra_edges=None, only_save_end_states=False, flexible_align=False):
        """
        Prepare minimisation and free energy lambda windows directory tree

        Parameters:
        -----------
        ligand_a_name: str
            name of first ligand
        ligand_b_name: str
            name of second ligand

        Return:
        -------
        """
        self.transformations = self.get_transformations(extra_edges)
        self.n_transformations = len(self.transformations)
        dataframe = self.transformations.copy()
        condition_1 = dataframe["ligand_a"] == ligand_a_name
        condition_2 = dataframe["ligand_b"] == ligand_b_name

        row_index = dataframe[condition_1 & condition_2].index.tolist()[0]

        lambda_strings = dataframe.iloc[row_index]["lambdas"].split(" ")
        lambda_values = list(map(float, lambda_strings))
   
        ligand_a = self.ligand_molecules[0]
        ligand_b = self.ligand_molecules[1]
        unbound = combine_unbound_ligands(ligand_a, ligand_b, flexible_align=flexible_align)
        bound = self.combine_bound_ligands(flexible_align=flexible_align)

        # Create ligand transformation directory tree in the first repeat directory, e.g. SOMD_1/lig_a~lig_b/ for bound/unbound
        # Only construct the BioSimSpace Relative AFE objects in the first repeat directory, to save on computation
        first_run_directory = self.output_directories[0]
        transformation_directory = first_run_directory + f"/{ligand_a_name}~{ligand_b_name}/"
        unbound_directory = transformation_directory + "/unbound/"
        bound_directory = transformation_directory + "/bound/"
        
        restart_interval = 500 #TODO add to config? or add function
        report_interval = 500 #TODO add to config? or add function

        free_energy_protocol = bss.Protocol.FreeEnergy(lam_vals=lambda_values, runtime=self.md_time, restart_interval=restart_interval, report_interval=report_interval)

        n_cycles = int(self.md_time._value) * 5 #TODO make editable; do 1 cycle per every 0.5 ns 
        n_moves = self.set_n_moves(number_of_cycles=n_cycles) # n_cycles * n_moves * timestep = runtime in ps

        n_frames = 250 #TODO add to config or add function
        buffered_coordinates_frequency = max(int(n_moves / n_frames), 10000) # https://github.com/OpenBioSim/sire/issues/113#issuecomment-1834317501

        cycles_per_saved_frame = max(1, restart_interval // n_moves) #Credit: Anna Herz https://github.com/michellab/BioSimSpace/blob/feature-amber-fep/python/BioSimSpace/_Config/_somd.py 

        config_options = {"ncycles": n_cycles, 
                          "nmoves": n_moves, 
                          "buffered coordinates frequency": buffered_coordinates_frequency, 
                          "ncycles_per_snap": cycles_per_saved_frame,
                          "minimal coordinate saving": only_save_end_states,
                        #   "cutoff distance": "8 angstrom", # Make editable? 
                          "minimise": True,
                          "minimise maximum iterations": self.min_steps}

        if self.cutoff_scheme == "pme":
            config_options["cutoff type"] = "PME"

        bss.FreeEnergy.Relative(system=unbound, protocol=free_energy_protocol, engine=self.md_engine, work_dir=unbound_directory, extra_options=config_options, setup_only=True)
        bss.FreeEnergy.Relative(system=bound, protocol=free_energy_protocol, engine=self.md_engine, work_dir=bound_directory, extra_options=config_options, setup_only=True)

        unbound_configurations = functions.get_files(unbound_directory + "/*/*.cfg")
        bound_configurations = functions.get_files(bound_directory + "/*/*.cfg")
        fix_afe_configurations(unbound_configurations)
        fix_afe_configurations(bound_configurations)

        # Copy lambda transformation directories from first repeat directory to the rest of the repeat directories, e.g. SOMD_2, SOMD_3
        _ = [os.system(f"cp -r {transformation_directory.rstrip('/')} {self.output_directories[i]}") for i in range(1, self.n_repeats)]
        

    def set_n_moves(self, number_of_cycles, stepsize=2):
        """
        Set SOMD nmoves in a reasonable way for better performance.
        See reference to: https://github.com/michellab/BioSimSpace/issues/258 and
        https://github.com/OpenBioSim/biosimspace/issues/18 

        Parameters:
        -----------
        stepsize: int
            MD stepsize in femtoseconds, default 2
        number_of_cycles: int
            number of MD cycles, default 5
        Return:
        -------
        number of moves: int
            number of moves for SOMD 
        """
         
        time = float(self.md_time._value) * 1000000
        number_of_steps = int(time / stepsize)
        return  number_of_steps // number_of_cycles


    def write_afe_run_script(self):
        """
        Write a AFE run script for given engine 

        Parameters:
        -----------

        Return:
        -------
        output: str
            full path to AFE run script
        """
        output = self.afe_input_directory + f"05_run_{self.md_engine}.sh"
        meze = os.environ["MEZEHOME"]
        
        template = meze + "/05_run_afe.sh"
        with open(template, "r") as file:
            lines = file.readlines()
        
        options = {"PATH_TO_LOGS": self.log_directory,
                   "N_TASKS": str(1), # change
                   "N_GPUS": str(1), # change
                   "N_CPUS": str(10), # change
                   "MEMORY": str(4069), # change
                #    "JOB": f"{self.md_engine}_afe", 
                   "ENGINE": self.md_engine,
                   "N_REPEATS": str(self.n_repeats),
                   "OUTPUTS_DIR": self.working_directory + "/outputs/"}

        # Credit: https://stackoverflow.com/a/51240945
        with open(output, "w") as file:
            for line in lines:
                for key, value in options.items():
                    line = line.replace(key, value)
                file.write(line)
        os.system(f"chmod +x {output}")
        return output
          

    def set_transformations(self, links_file=None):
        """
        Create a transformation network with LOMAP

        Return:
        -------
        transformations: pd.DataFrame
            forward and backward transformations with their associated lomap scores, number of windows and lambda list
        """
        if links_file:
            links_filename = functions.get_filename(links_file)
            save_name = self.afe_input_directory + f"{links_filename}.csv"
            links = pd.read_csv(links_file, sep="\s", names=["ligA", "ligB"])
            first_column = links["ligA"].tolist()
            second_column = links["ligB"].tolist()
            lomap_work_directory = self.check_lomap_directory()
            scores = [0.5] * len(first_column)
            start_ligand = [ligand for ligand in first_column]
            end_ligand = [ligand for ligand in second_column] 

            dataframe = pd.DataFrame()
            dataframe["ligand_a"] = start_ligand
            dataframe["ligand_b"] = end_ligand
            dataframe["score"] = scores
            for score in scores:
                n_windows = self.set_n_windows(score)
                self.n_windows.append(n_windows)
                lambda_windows = create_lambda_windows(n_windows)
                self.lambdas.append(lambda_windows)
            dataframe["n_windows"] = self.n_windows
            dataframe["lambdas"] = self.lambdas
            dataframe.to_csv(save_name)
            
        else:
            save_name = self.afe_input_directory+f"/meze_network.csv"
            lomap_work_directory = self.check_lomap_directory()
            transformations, lomap_scores = bss.Align.generateNetwork(self.ligand_molecules, plot_network=True, names=self.names, work_dir=lomap_work_directory)   
        
        if not os.path.isfile(save_name):
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
            dataframe.to_csv(save_name)
        else:
            dataframe = pd.read_csv(save_name)
        return dataframe, save_name
    

    def get_transformations(self, extra_edges=None):
        """
        Getter for transformations csv file

        Parameters:
        -----------

        Return:
        -------
        df: pd.DataFrame: 
            meze network as a csv file
        """
        if not extra_edges:
            df = pd.read_csv(self.afe_input_directory+f"/meze_network.csv", header=0, index_col=0)
        else:
            df = pd.read_csv(extra_edges, header=0, index_col=0)
        return df
    

    def get_files(self):
        """
        Read in all .sdf and .mol2 files in the provided path.

        Return:
        -------
        list
            list of ligand filenames
        """
        # Adapted from dbmol.py:
        ligand_files = functions.get_files(f"{self.ligand_path}/*.sdf")
        ligand_files += functions.get_files(f"{self.ligand_path}/*.mol2")        
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
        return(len(self.ligands))


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

        if os.path.isfile(lomap_work_directory + "/outputs/network.png"):
            os.rename(lomap_work_directory + "/outputs/network.png", lomap_work_directory + "initial_network.png")

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
        if not self.is_md:
            strip = self.output_directories[0].split("/")[-2]
            path_to_outputs = self.output_directories[0].replace(strip, "")
        else: 
            path_to_outputs = self.outputs

        protocol = [f"group name = {self.group_name}",
                    f"ligand forcefield = {self.ligand_forcefield}", 
                    f"ligand charge = {self.ligand_charge}",
                    f"prepared protein file = {self.prepared_protein}",
                    f"protein input file = {self.protein_file}",
                    f"protein forcefield = {self.protein_forcefield}", 
                    f"cutoff scheme = {self.cutoff_scheme}",
                    f"water model = {self.water_model}", 
                    f"box edges = {self.box_edges}", # in angstrom 
                    f"box shape = {self.box_shape}", 
                    f"solvation method = {self.solvation_method}",
                    f"solvent closeness = {self.solvent_closeness}",
                    f"minimisation steps = {self.min_steps}",
                    f"minimisation stepsize = {self.min_dt}",
                    f"minimisation tolerance = {self.min_tol}",

                    f"temperature = {self.temperature._value}",
                    f"pressure = {self.pressure._value}",
                    f"sampling time = {self.md_time._value}",
                    f"engine = {self.md_engine}",
                    f"only save end states = {self.only_save_end_states}",
                    f"outputs = {path_to_outputs}",
                    f"repeats = {self.n_repeats}",

                    f"project directory = {self.working_directory}",
                    f"equilibration directory = {self.equilibration_directory}",
                    f"ligand directory = {self.ligand_path}",
                    f"protein directory = {self.protein_path}",
                    f"log directory = {self.log_directory}",
                    f"plots directory = {self.plots}"]
        
        if not self.is_md:
            protocol.append(f"afe input directory = {self.afe_input_directory}",
                            f"short nvt = {self.short_nvt._value}",
                            f"nvt = {self.nvt._value}",
                            f"npt = {self.npt._value}",
                            f"network file = {self.network_file}",)
            protocol_file = self.afe_input_directory + "/protocol.dat"
        else:
            protocol_file = self.md_input_directory + "/protocol.dat"

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
       

    def get_name(self, given_group_name):
        """
        Get group name. If user group name is not given, infer name from input protein filename.

        Parameters:
        -----------

        Return:
        -------
        str
            group name
        """
        if given_group_name:
            return given_group_name
        else:
            if len(self.protein_file) > 1:
                return self.protein_file[0].split("/")[-1].split(".")[0]
            else:
                return self.protein_file.split("/")[-1].split(".")[0]


def main():
   pass

if __name__ == "__main__":
    main()