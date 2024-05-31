import argparse
import pathlib
import sofra
import os
import BioSimSpace as bss
import numpy as np
import functions
from definitions import PICOSECOND, KELVIN, NANOSECOND
import csv
import MDAnalysis as mda
import MDAnalysis.analysis.distances 
import shutil
from BioSimSpace import _Exceptions


def residue_restraint_mask(residue_ids):
    """
    Overriden method from BioSimSpace: Instead of atom indices (1 indexed),
    get a restraint mask for residue ids.
    Internal helper function to create an AMBER restraint mask from a
    list of residue ids.
    Reference: https://github.com/OpenBioSim/biosimspace/blob/devel/python/BioSimSpace/_Config/_amber.py#L203

    Parameters
    ----------
    residue_ids: [int]
        A list of residue ids.

    Returns
    -------
    restraint_mask: str
        The AMBER restraint mask.
    """

    residue_ids = list(set(residue_ids))
    residue_ids.sort()

    # Handle single residue restraints differently.
    if len(residue_ids) == 1:
        restraint_mask = f"{residue_ids[0]}"

    else:
        # Start the mask with the first residue id.
        restraint_mask = f"{residue_ids[0]}"

        # Store the current id.
        previous_id = residue_ids[0]

        # Store the lead id for this block.
        lead_id = previous_id

        # Loop over all other indices.
        for idx in residue_ids[1:]:
            # There is a gap in the indices.
            if idx - previous_id > 1:
                if previous_id != lead_id:
                    restraint_mask += f"{previous_id},{idx+1}"
                else:
                    restraint_mask += f",{idx}"
                lead_id = idx
            else:
                # This is the first index beyond the lead.
                if idx - lead_id == 1:
                    restraint_mask += "-"
            # Set the value of the previous index.
            previous_id = idx

        # Add the final residue to the mask.
        if idx - residue_ids[-2] == 1:
            restraint_mask += f"{idx}"
        else:
            if idx != lead_id:
                restraint_mask += f",{idx}"

    return restraint_mask


class Meze(sofra.Sofra):

    #TODO change the afe_input_path=os.getcwd()+"/afe/" (etc.) somehow 
    # it should not be taken from the cwd as this raises an error if we run any of the run scripts from the afe dir
    # it's annoying to have to put this in as an argument to the constructor everytime

    def __init__(self, protein_file, prepared=False, metal="ZN", cut_off=2.6, force_constant_0=100, water_file=None, 
                 workdir=os.getcwd(), is_md=False, md_input_directory=None, restraints=True,
                 afe_input_path=os.getcwd()+"/afe/", equilibration_path=os.getcwd()+"/equilibration/", outputs=os.getcwd()+"/outputs/", log_directory=os.getcwd()+"/logs/",
                 ligand_path=os.getcwd()+"/inputs/ligands/", ligand_charge=0, ligand_ff="gaff2", 
                 group_name=None, protein_path=os.getcwd()+"/inputs/protein/", water_model="tip3p", protein_ff="ff14SB", 
                 engine="SOMD", sampling_time=4, box_edges=20, box_shape="cubic", min_steps=5000, short_nvt=50, nvt=1, npt=200, 
                 min_dt=0.01, min_tol=1000, repeats=3, temperature=300, pressure=1, threshold=0.4, n_normal=11, n_difficult=17,
                 cutoff_scheme="rf", solvation_method="gromacs", solvent_closeness=1.0, only_save_end_states=False):
        
        self.protein_file = protein_file

        super().__init__(prepared=prepared, workdir=workdir, ligand_path=ligand_path, group_name=group_name, protein_file=protein_file, protein_path=protein_path,
                         water_model=water_model, ligand_ff=ligand_ff, protein_ff=protein_ff, ligand_charge=ligand_charge, is_md=is_md,md_input_directory=md_input_directory,
                         afe_input_path=afe_input_path, equilibration_path=equilibration_path, outputs=outputs, log_directory=log_directory,
                         engine=engine, sampling_time=sampling_time, box_edges=box_edges, box_shape=box_shape, min_steps=min_steps, short_nvt=short_nvt, nvt=nvt, npt=npt, 
                         min_dt=min_dt, min_tol=min_tol, repeats=repeats, temperature=temperature, pressure=pressure, threshold=threshold, n_normal=n_normal, n_difficult=n_difficult,
                         cutoff_scheme=cutoff_scheme, solvation_method=solvation_method, solvent_closeness=solvent_closeness, only_save_end_states=only_save_end_states)
        
        self.restraints = restraints
        
        if self.prepared:
            extensions = [functions.get_file_extension(file) for file in self.protein_file]
            for extension in extensions:
                if extension not in ["rst7", "prm7", "gro", "top"]:
                    raise ValueError("The prepared file must have a topology file in either PARM7 or Gro format!")
            
            if "prm7" in extensions:
                topology_format = "PARM7"
                topology_file = [file for file in self.protein_file if functions.get_file_extension(file) == "prm7"][0]

                coordinate_file = [file for file in self.protein_file if functions.get_file_extension(file) == "rst7"][0]

                self.universe = mda.Universe(topology=topology_file, topology_format=topology_format,
                                             coordinates=coordinate_file)
            elif "gro" in extensions:
                file = [file for file in self.protein_file if functions.get_file_extension(file) == "gro"][0]
                self.universe = mda.Universe(file)
            # for file in self.protein_file:
            #     extension = functions.get_file_extension(file)
            #     if extension.lower() == "prm7":
            #         topology_format = "PARM7"
            #         topology_file = file
            #         # self.universe = mda.Universe(topology=topology_file, topology_format=topology_format)
            #     elif extension.lower() == "gro":
            #         coordinate_format = extension
            #         coordinate_file = file

            # self.universe = mda.Universe(topology=topology_format, topology_format=topology_format, coordinates=coordinate_file, format=coordinate_format)
        else:
            self.universe = mda.Universe(self.protein_file, format="pdb")
        self.force_constant_0 = functions.check_float(force_constant_0)
        self.cut_off = functions.check_positive(functions.check_float(cut_off))
        
        if water_file:
            self.xtal_water = bss.IO.readMolecules(water_file)

        if metal.upper() == "ZN":
            self.metal_resname = metal
            self.metals = self.universe.select_atoms(f"resname {self.metal_resname}")
            self.metal_resids = self.metals.resids
            self.metal_atomids = self.metals.atoms.ids
        else:
            print(f"Metal {metal} is not supported yet.")
        
    
    def prepare_meze(self):
        """
        Prepare QM/MM calculations by creating ligand and protocol dat files.

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
        self.protocol_file = self.create_qm_protocol_file()
        self.prepared = True
        return self


    def create_protocol_file(self):
        """
        Add metal protein information to the AFE protocol file 

        Parameters:
        -----------
        Network: Network
            Network class object

        Return:
        -------
        protocol_file: str
            protocol datafile
        """
        protocol_file = super().create_protocol_file()
        metal_protocol = [f"metal = {self.metal_resname}",
                          f"cutoff = {self.cut_off}",
                          f"force constant = {self.force_constant_0}"]

        with open(protocol_file, "a") as file:
            writer = csv.writer(file)
            for protocol_line in metal_protocol:
                writer.writerow([protocol_line])
        return protocol_file
            

    def add_somd_restraints(self, directory):
        
        ligand_a = self.ligands[0].name
        amber_restraints_file = functions.get_files(self.equilibration_directory + f"bound/{ligand_a}/*.RST")[0]
        somd_restraints_dict = {}
 
        with open(amber_restraints_file, "r") as f:
            for line in f:
                if "#" not in line:
                    iat1 = int(line.split("iat=")[1].split(",")[0]) - 1
                    iat2 = int(line.split("iat=")[1].split(",")[1]) - 1 
                    r2 = float(line.split("r2=")[1].split(",")[0])
                    r3 = float(line.split("r3=")[1].split(",")[0])
                    rk2 = float(line.split("rk2=")[1].split(",")[0])
                    flat_bottom_radius = np.round((r3 - r2)/2, decimals=2)
                    equilibrium_distance = np.round(r2 + flat_bottom_radius, decimals=2)
                    force_constant = np.round(rk2, decimals=2)

                    atom_key = (iat1, iat2)
                    restraint_value = (equilibrium_distance, force_constant, flat_bottom_radius)
                    somd_restraints_dict[atom_key] = restraint_value

        somd_restraints_file = self.somd_restraints(directory, somd_restraints_dict)
        with open(somd_restraints_file, "r") as file:
            restraints = file.readlines()

        with open(directory + "somd.cfg", "a") as file:
            file.writelines(restraints)


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
        merged_ligands = sofra.merge_ligands(ligand_1, ligand_2, flexible_align=flexible_align)

        removal_system = system_1.copy()
        removal_system.removeWaterMolecules()
        removal_system.removeMolecules(ligand_1)
        zn_and_salts = removal_system - protein

        ion_molecules = []
        for res in zn_and_salts.getResidues():
            if res.name() != self.metal_resname:
                removal_system.removeMolecules(res.toMolecule())
                ion_molecules.append(res.toMolecule())

        water_molecules = system_1.getWaterMolecules()
        final_system = removal_system.copy()
        final_system.addMolecules(merged_ligands)
        final_system.addMolecules(water_molecules)
        final_system.addMolecules(ion_molecules)
        return final_system


    def prepare_afe(self, ligand_a_name, ligand_b_name, extra_edges=None, only_save_end_states=False, flexible_align=False):
        """
        Inherited method from sofra.Network; adds restraints to somd config files

        Parameters:
        -----------
        ligand_a_name: str
            name of the first ligand of the transformation
        ligand_b_name: str
            name of the second ligand of the transformation

        Return:
        -------
        """
        
        super().prepare_afe(ligand_a_name, ligand_b_name, extra_edges=extra_edges, only_save_end_states=only_save_end_states, flexible_align=flexible_align)
        
        if self.restraints:
            first_run_directory = self.output_directories[0]
            transformation_directory = first_run_directory + f"/{ligand_a_name}~{ligand_b_name}/"
            bound_directory = transformation_directory + "/bound/" # unbound doesn't have restraints
            lambda_directories = functions.get_files(bound_directory + "lambda_*/")
            for i in range(len(lambda_directories)):
                self.add_somd_restraints(lambda_directories[i])
            _ = [os.system(f"cp -r {transformation_directory.rstrip('/')} {self.output_directories[i]}") for i in range(1, self.n_repeats)]


    def set_universe(self, file_name):
        """
        Take a file name and create an MDAnalysis Universe from prm7 and rst7 files.
        Set the universe attribute.

        Parameters:
        -----------
        file_name: str
            file name without extension

        Return:
        -------
        self: Meze.Meze()
            Meze class object
        """
        self.universe = mda.Universe(file_name + ".prm7", file_name + ".rst7",
                                     topology_format="PARM7", format="RESTRT")
        self.metals = self.universe.select_atoms(f"resname {self.metal_resname}")
        self.metal_resids = self.metals.resids
        self.metal_atomids = self.metals.atoms.ids
        return self
    

    def get_active_site(self):
        """
        Find active site residues defined as around metal ions defined by the given cutoff value.

        Parameters:
        -----------

        Return:
        -------
        mda.AtomGroup
            active site atoms
        """
        selection = f"resname {self.metal_resname} or around {self.cut_off} (resname {self.metal_resname})"
        return self.universe.select_atoms(selection)


    def get_metal_ligands(self):
        """
        Create a dictionary of metal ion(s) and their coordinating ligand(s)

        Parameters:
        -----------

        Return:
        -------
        metal_ligands: dict
            dictionary where keys are metal ions (differentiated by subscripts) and 
            values are the atom groups containing coordinating ligands for that ion
        """
        metal_ligands = {}
        for i in range(len(self.metal_resids)):
            ligands = self.universe.select_atoms(f"not (name {self.metal_resname} or element H) and sphzone {self.cut_off} (resid {self.metal_resids[i]})")
            key = self.metal_atomids[i] 
            metal_ligands[key] = ligands
        return metal_ligands


    def amber_restraints(self, workdir, restraints):
        """
        Write amber style restraints

        Parameters:
        -----------
        workdir: str
            output path for saving restraints file
        restraints: dict
            dictionary where keys are metal ions (differentiated by subscripts) and 
            values are the atom groups containing coordinating ligands for that ion
            
        Return:
        -------
        output_file: str
            amber style restraints file
        """
        atom_ids = []
        r1, r2, r3, r4 = [], [], [], []

        for key, value in restraints.items():
            metal_id = key[0]
            ligating_atom_id = key[1]
            distance = value[0]
            flat_bottom_radius = value[-1]
            atom_ids.append(f"iat={metal_id},{ligating_atom_id}")

            linear_response_region_lower_bound = np.round(distance - flat_bottom_radius, decimals=2)
            flat_region_lower_bound = np.round(distance - (flat_bottom_radius / 2), decimals=2)
            flat_region_upper_bound = np.round(distance + (flat_bottom_radius / 2), decimals=2)
            linear_response_region_upper_bound = np.round(distance + flat_bottom_radius, decimals=2)
            r1.append("r1="+str(linear_response_region_lower_bound))
            r2.append("r2="+str(flat_region_lower_bound))
            r3.append("r3="+str(flat_region_upper_bound))
            r4.append("r4="+str(linear_response_region_upper_bound))

        rk2 = ["rk2="+str(self.force_constant_0)] * len(atom_ids)
        rk3 = ["rk3="+str(self.force_constant_0)] * len(atom_ids)

        output_file = workdir + "restraints_0.RST"
        with open(output_file, "w") as file:
            file.write(f"# Harmonic bond restraints between {self.metal_resname} and coordinating protein residues\n")
            for i in range(len(atom_ids)):
                line = f"&rst {atom_ids[i]}, {r1[i]}, {r2[i]}, {r3[i]}, {r4[i]}, {rk2[i]}, {rk3[i]},/\n"
                file.write(line)
        return output_file


    def build_restraints(self):
        """
        Build a dictionary of restraints for model 0

        Parameters:
        -----------

        Return:
        -------
        restraints: dict    
            keys: tuple(metal_id, ligating_atom_id), values: tuple(eq_distance, force_constant, flat_bottom_radius)
        """
        metal_ligands = self.get_metal_ligands()
        protein = self.universe.select_atoms("protein")
        restraints = {}
        for i in range(len(self.metal_atomids)):
            metal_id = self.metal_atomids[i]
            for ligating_atom in metal_ligands[metal_id]:
                if ligating_atom in protein or ligating_atom.resname != "WAT" and ligating_atom.resname != "MOL": # don't restrain ligand or water!
                    key = (metal_id, ligating_atom.id)
                    atom_group_1 = self.universe.select_atoms(f"resid {self.metal_resids[i]}")
                    atom_group_2 = self.universe.select_atoms(f"resid {ligating_atom.resid} and name {ligating_atom.name}")
                    distance = MDAnalysis.analysis.distances.dist(atom_group_1, atom_group_2)[-1][0]
                    equilibrium_distance = np.round(distance, decimals=2)
                    force_constant = np.round(self.force_constant_0, decimals=2)
                    flat_bottom_radius = 1.00 #TODO Make editable?
                    value = (equilibrium_distance, force_constant, flat_bottom_radius)
                    restraints[key] = value
        return restraints


    def write_restraints_file_0(self, workdir, engine="amber"):
        """
        Write an Amber-compatible restraint file for model 0. 

        Parameters:
        -----------
        engine: str
            name of the MD engine used; current choices include amber and somd

        Return:
        -------
        output_file: str
            amber-format restraints file for model 0
        """
        restraints = self.build_restraints()
        if engine == "amber":
            output_file = self.amber_restraints(workdir, restraints)
        elif engine == "somd":
            output_file = self.somd_restraints(workdir, restraints)  
        return output_file
    

    def somd_restraints(self, workdir, restraints):
        """
        Write somd restraints to be added to somd.cfg

        Parameters:
        -----------
        restraints: dict
            dictionary where keys are metal ions (differentiated by subscripts) and 
            values are the atom groups containing coordinating ligands for that ion

        Return:
        -------
        output_file: str
            somd style restraints file to be added to somd.cfg
        """
        output_file = workdir + "somd.restraints"

        with open(output_file, "w") as file:
            file.write("use permanent distance restraints = True\n")
            file.write(f"permanent distance restraints dictionary = {str(restraints)}")
        return output_file


def clean_arguments(arguments):
    """
    Check arguments and clean them

    Parameters:
    -----------
    arguments: Namespace
        command-line arguments
    
    Return:
    -------
    cleaned Namespace of arguments
    
    """
    check_integers = lambda arg: functions.check_positive(functions.check_int(arg))
    check_integers(arguments.min_steps)
    check_floats = lambda arg: functions.check_positive(functions.check_float(arg))
    check_floats(arguments.short_nvt)
    check_floats(arguments.nvt)
    check_floats(arguments.npt)
    return arguments


def character(string):
    if len(string) > 1:
        message = f"Input must be a single character, not {string}"
        raise argparse.ArgumentTypeError(message)
    if string == "_":
        message = "The character _ is not an allowed separator"
        raise argparse.ArgumentTypeError(message)


def main():

    parser = argparse.ArgumentParser(description="MEZE: MetalloEnZymE FF-builder for alchemistry")

    parser.add_argument("protocol_file",
                        help="protocol file containing equilibration options",
                        type=str,
                        default=os.getcwd() + "/afe/protocol.dat")
    
    parser.add_argument("transformation",
                        help="the pair of ligands undergoing AFE transformation, e.g. ligand_1~ligand_2",
                        type=str)
    
    parser.add_argument("-f",
                        "--flexible-align",
                        dest="flexible_align",
                        help="use bss.Align.flexAlign to merge ligands",
                        action=argparse.BooleanOptionalAction)

    parser.add_argument("-et",
                        "--extra-transformations",
                        dest="extra_transformations_file",
                        help="file containing additional transformations in format: lig A lig B; equivalent to BioSimSpace.generateNetwork links file",
                        type=str)
    
    parser.add_argument("-s",
                        "--separator",
                        help="character separating the two ligand names",
                        default="~",
                        type=character)
    
    parser.add_argument("--no-restraints",
                    dest="no_restraints",
                    help="do not apply harmonic restraints on metal-coordinating residues",
                    action="store_true")

    
    arguments = parser.parse_args()
    
    protocol = functions.input_to_dict(arguments.protocol_file)

    keys = list(protocol.keys())
    if "metal" in keys:
        metal = True
    else:
        metal = False

    if arguments.no_restraints:
        apply_restraints = False
    else:
        apply_restraints = True

    if metal:
        meze = Meze(prepared=True,
                    restraints=apply_restraints,
                    protein_file=protocol["prepared protein file"],
                    cut_off=protocol["cutoff"],
                    force_constant_0=protocol["force constant"],
                    workdir=protocol["project directory"],
                    equilibration_path=protocol["equilibration directory"],
                    afe_input_path=protocol["afe input directory"],
                    log_directory=protocol["log directory"],
                    outputs=protocol["outputs"],
                    ligand_path=protocol["ligand directory"],
                    ligand_charge=protocol["ligand charge"],
                    ligand_ff=protocol["ligand forcefield"],
                    group_name=protocol["group name"],        
                    protein_path=protocol["protein directory"],
                    water_model=protocol["water model"],
                    protein_ff=protocol["protein forcefield"],
                    engine=protocol["engine"],
                    sampling_time=protocol["sampling time"],
                    box_edges=protocol["box edges"],
                    box_shape=protocol["box shape"],
                    min_steps=protocol["minimisation steps"],
                    short_nvt=protocol["short nvt"],
                    nvt=protocol["nvt"],
                    npt=protocol["npt"],
                    min_dt=protocol["minimisation stepsize"],
                    min_tol=protocol["minimisation tolerance"],
                    repeats=protocol["repeats"],
                    temperature=protocol["temperature"],
                    pressure=protocol["pressure"],
                    only_save_end_states=protocol["only save end states"])
        
    elif not metal:

        meze = sofra.Sofra(prepared=True,
                           equilibration_path=protocol["equilibration directory"],
                           outputs=protocol["outputs"],
                           workdir=protocol["project directory"],
                           log_directory=protocol["log directory"],
                           afe_input_path=protocol["afe input directory"],
                           ligand_path=protocol["ligand directory"],
                           group_name=protocol["group name"],
                           protein_file=protocol["prepared protein file"],
                           protein_path=protocol["protein directory"],
                           water_model=protocol["water model"],
                           ligand_ff=protocol["ligand forcefield"],
                           protein_ff=protocol["protein forcefield"],
                           ligand_charge=protocol["ligand charge"],
                           engine=protocol["engine"],
                           sampling_time=protocol["sampling time"],
                           box_edges=protocol["box edges"],
                           box_shape=protocol["box shape"],
                           min_steps=protocol["minimisation steps"],
                           short_nvt=protocol["short nvt"],
                           nvt=protocol["nvt"],
                           npt=protocol["npt"],
                           min_dt=protocol["minimisation stepsize"],
                           min_tol=protocol["minimisation tolerance"],
                           repeats=protocol["repeats"],
                           temperature=protocol["temperature"],
                           pressure=protocol["pressure"],
                           only_save_end_states=protocol["only save end states"])
          
    ligand_a, ligand_b = functions.separate(arguments.transformation)

    equilibrated_network = meze.get_equilibrated(ligand_a, ligand_b)

    equilibrated_network.prepare_afe(ligand_a, ligand_b, extra_edges=arguments.extra_transformations_file, only_save_end_states=meze.only_save_end_states, flexible_align=arguments.flexible_align) 

if __name__ == "__main__":
    main()
