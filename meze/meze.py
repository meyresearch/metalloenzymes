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
                 workdir=os.getcwd(), afe_input_path=os.getcwd()+"/afe/", equilibration_path=os.getcwd()+"/equilibration/", outputs=os.getcwd()+"/outputs/",
                 ligand_path=os.getcwd()+"/inputs/ligands/", ligand_charge=0, ligand_ff="gaff2", 
                 group_name=None, protein_path=os.getcwd()+"/inputs/protein/", water_model="tip3p", protein_ff="ff14SB", 
                 engine="SOMD", sampling_time=4, box_edges=20, box_shape="cubic", min_steps=5000, short_nvt=50, nvt=1, npt=200, 
                 min_dt=0.01, min_tol=1000, repeats=3, temperature=300, pressure=1, threshold=0.4, n_normal=11, n_difficult=17,
                 solvation_method="gromacs", solvent_closeness=1.0):
        
        self.protein_file = protein_file

        super().__init__(prepared=prepared, workdir=workdir, ligand_path=ligand_path, group_name=group_name, protein_file=protein_file, protein_path=protein_path, 
                         water_model=water_model, ligand_ff=ligand_ff, protein_ff=protein_ff, ligand_charge=ligand_charge, 
                         afe_input_path=afe_input_path, equilibration_path=equilibration_path, outputs=outputs,
                         engine=engine, sampling_time=sampling_time, box_edges=box_edges, box_shape=box_shape, min_steps=min_steps, short_nvt=short_nvt, nvt=nvt, npt=npt, 
                         min_dt=min_dt, min_tol=min_tol, repeats=repeats, temperature=temperature, pressure=pressure, threshold=threshold, n_normal=n_normal, n_difficult=n_difficult,
                         solvation_method=solvation_method, solvent_closeness=solvent_closeness)
        
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
        # def convert_amber_restraints_to_somd(restraints_file):
        #   pass        
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


    def combine_bound_ligands(self):
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
        merged_ligands = sofra.merge_ligands(ligand_1, ligand_2)

        removal_system = system_1.copy()
        removal_system.removeWaterMolecules()
        protein = removal_system.getMolecules()[0]
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


    def prepare_afe(self, ligand_a_name, ligand_b_name, extra_edges=None):
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
        
        super().prepare_afe(ligand_a_name, ligand_b_name, extra_edges=extra_edges)
        first_run_directory = self.output_directories[0]
        transformation_directory = first_run_directory + f"/{ligand_a_name}~{ligand_b_name}/"
        bound_directory = transformation_directory + "/bound/" # unbound doesn't have restraints
        lambda_directories = functions.get_files(bound_directory + "lambda_*/")
        minimisation_directories = functions.get_files(bound_directory + "minimisation/*/")
        # lambda_config_path = self.outputs + f"*/{ligand_a_name}~{ligand_b_name}/*/*/*.cfg"
        # lambda_minimisation_config_path =  self.outputs + f"*/{ligand_a_name}~{ligand_b_name}/*/*/*/*.cfg"
        for i in range(len(lambda_directories)):
            self.add_somd_restraints(lambda_directories[i])
            self.add_somd_restraints(minimisation_directories[i])   
            
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


    def model_0(self, ligand_name):

        if self.is_qm:
            self.qmmm_minimisation(ligand_name)
            self.qmmm_equilibration(ligand_name)
            self.qmmm_production(ligand_name)

        # elif not self.is_qm:
        #     minimised_system = self.minimisation_0(ligand_name)
        #     equilibrated_system = self.equilibration_0(ligand_name, minimised_system)
        #     self.production_0(ligand_name, equilibrated_system)


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
                if ligating_atom in protein or ligating_atom.resname == "WAT" and ligating_atom.resname != "MOL": # don't restrain ligand
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


    # Potentially to be depracated:

    def get_qm_region(self):
        """
        Take an atom group of the active site and return a dictionary of the QM region selections.

        Parameters:
        -----------
        active_site: mda.AtomGroup
            group of atoms contained in the active site

        Return:
        -------
        qm_region: dict
            qm region selection for QM/MM
        """
        active_site = self.get_active_site()
        protein = self.universe.select_atoms("protein")
        qm_region = {}
        whole_residues = []
        atom_ids = []
        for residue in active_site.residues:
            if residue in protein.residues:
                selection = get_side_chain_selection(self.universe, residue)
                atom_ids.append(selection)
            else:
                whole_residues.append(residue.resid)
        formatted_whole_residues = residue_restraint_mask(whole_residues)
        qm_region["whole_residues"] = formatted_whole_residues
        qm_region["atom_ids"] = atom_ids
        return qm_region
    

    def get_qm_charge(self):
        """
        Determine the net charge of the QM region.

        Parameters:
        -----------

        Return:
        -------
        int:
            net integer charge of QM region
        """
        active_site = self.get_active_site()
        protein = self.universe.select_atoms("protein")
        qm_region = []
        for residue in active_site.residues:
            if residue in protein.residues:
                n_terminus = "name N or name H"
                alpha_carbon = "name CA or name HA"
                c_terminus = "name C or name O"
                atoms_in_residue = self.universe.select_atoms(f"resid {residue.resid}")
                qm_region_for_residue = atoms_in_residue.select_atoms(f"not ({n_terminus} or {alpha_carbon} or {c_terminus})")
                qm_region.append(qm_region_for_residue)
            else:
                qm_residue = self.universe.select_atoms(f"resid {residue.resid}")
                qm_region.append(qm_residue)
        charges = [round(atom_group.charges.sum()) for atom_group in qm_region]
        return sum(charges)

    
    def get_dftb3_options(self, qm_region):
        """
        Create QM/MM options for DFTB3 to be used in QM/MM input files

        Parameters:
        -----------
        qm_region: dict
            dictionary of the qm region selections

        Return:
        -------
        options: dict
            DFTB3 options
        """
        whole_residues = qm_region["whole_residues"]
        atom_ids = ",".join(qm_region["atom_ids"])
        options = {"qmmask": f"':{whole_residues}|(@{atom_ids})'",
                   "writepdb": "1",
                   "qmcharge": str(self.get_qm_charge()),
                   "qm_theory": "'DFTB3'",
                   "qmshake": "0",
                   "qm_ewald": "1",
                   "qm_pme": "1"}
        return options
    
    
    def create_qm_protocol_file(self):
            """
            Create protocol.dat file for QM/MM runs

            Parameters:
            -----------
            Network: Network
                Network class object

            Return:
            -------
            protocol_file: str
                protocol datafile
            """
            protocol = [f"metal = {self.metal_resname}",
                        f"cut-off = {self.cut_off}",
                        f"force constant 0 = {self.force_constant_0}",
                        f"group name = {self.group_name}",
                        f"ligand forcefield = {self.ligand_forcefield}", 
                        f"ligand charge = {self.ligand_charge}",
                        f"prepared protein file = {self.prepared_protein}",
                        f"protein input file = {self.protein_file}",
                        f"protein forcefield = {self.protein_forcefield}", 
                        f"water model = {self.water_model}", 
                        f"box edges = {self.box_edges}", # in angstrom 
                        f"box shape = {self.box_shape}", 
                        f"minimisation steps = {self.min_steps}",
                        f"minimisation stepsize = {self.min_dt}",
                        f"minimisation tolerance = {self.min_tol}",
                        f"short nvt = {self.short_nvt._value}",
                        f"nvt = {self.nvt._value}",
                        f"npt = {self.npt._value}",
                        f"temperature = {self.temperature._value}",
                        f"pressure = {self.pressure._value}",
                        f"sampling time = {self.md_time._value}",
                        f"engine = {self.md_engine}",
                        f"outputs = {self.output_directory}",
                        f"repeats = {self.n_repeats}",
                        f"project directory = {self.working_directory}",
                        f"equilibration directory = {self.equilibration_directory}",
                        f"ligand directory = {self.ligand_path}",
                        f"protein directory = {self.protein_path}",
                        f"log directory = {self.log_directory}",
                        f"input directory = {self.input_directory}"]

            protocol_file = self.input_directory + "/protocol.dat"
            with open(protocol_file, "w") as file:
                writer = csv.writer(file)
                for protocol_line in protocol:
                    writer.writerow([protocol_line])
            return protocol_file

    def production_0(self, ligand_name, equilibrated_system, nonbonded_cut_off=9.0, dt=0.002, runtime=50):

        directory = self.output_directory+f"{ligand_name}/"
        production_directory = functions.mkdir(directory)

        template_restraints_file = self.write_restraints_file_0()
        restraints_file = shutil.copy(template_restraints_file, production_directory).split("/")[-1]
        output_frequency = runtime // dt
        production_options = {"ioutfm": 1,
                              "cut": nonbonded_cut_off,
                              "gamma_ln": 1.0,
                              "ntp": 1,
                              "ntb": 2, 
                              "ntc": 2,
                              "ntf": 2, 
                              "iwrap": 1,  
                              "nmropt": 1}
        
        runtime_ns = functions.convert_to_units(runtime, NANOSECOND)
        namelist = ["&wt TYPE='DUMPFREQ', istep1=1 /"]

        protocol = bss.Protocol.Production(timestep=dt*PICOSECOND,
                                           runtime=runtime_ns,
                                           temperature=self.temperature,
                                           report_interval=output_frequency,
                                           restart_interval=output_frequency,
                                           restart=True)
        
        process = bss.Process.Amber(system=equilibrated_system, 
                                    protocol=protocol, 
                                    name="md",
                                    work_dir=production_directory,
                                    extra_options=production_options,
                                    extra_lines=namelist)
        config = production_directory + "/*.cfg"
        config_file = functions.get_files(config)[0]

        with open(config_file, "a") as file:
            file.write("\n")
            file.write(f"DISANG={restraints_file}\n")
            file.write(f"DUMPAVE=distances.out\n")

        process.start()
        process.wait()
        if process.isError():
            print(process.stdout())
            print(process.stderr())
            raise bss._Exceptions.ThirdPartyError("The process exited with an error!")

    
    def equilibration_0(self, ligand_name, minimised_system, nonbonded_cut_off=8.0, dt=0.001, runtime=1,
                        start_temperature=100): # in ns
        
        directory = self.equilibration_directory+f"{ligand_name}/"
        directories = lambda step: functions.mkdir(directory + step)
        
        heat_02_dir = directories("02_heat")
        relax_03_dir = directories("03_relax")
        lower_04_dir = directories("04_lower")
        bb_min_05_dir = directories("05_bb_min")
        relax_06_dir = directories("06_relax")
        reduce_07_dir = directories("07_reduce")
        continue_08_dir = directories("08_continue")
        relax_09_dir = directories("09_relax")

        template_restraints_file = self.write_restraints_file_0()

        restraints_files = []
        for directory in [heat_02_dir, relax_03_dir, lower_04_dir, bb_min_05_dir, relax_06_dir, reduce_07_dir, continue_08_dir, relax_09_dir]:
            restraints_files.append(shutil.copy(template_restraints_file, directory).split("/")[-1])

        max_cycles = self.min_steps        
        output_frequency = max_cycles // 10

        heat_02_options = {"ioutfm": 1,
                         "cut": nonbonded_cut_off,
                         "tol": 0.00001,
                         "nscm": 0,
                         "barostat": 2,
                         "taup": 1.0,
                         "gamma_ln": 1.0,
                         "ntp": 0,
                         "ntb": 1, 
                         "ntc": 2,
                         "ntf": 2,   
                         "nmropt": 1}
        
        runtime_ns = functions.convert_to_units(runtime, NANOSECOND)
        namelist = ["&wt TYPE='DUMPFREQ', istep1=1 /"]
        heat_02_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                    runtime=runtime_ns, 
                                                    temperature_start=start_temperature*KELVIN, 
                                                    temperature_end=self.temperature, 
                                                    report_interval=output_frequency, 
                                                    restart_interval=output_frequency,
                                                    restraint="heavy",
                                                    force_constant=self.force_constant_0)
        amber_home = os.environ["AMBERHOME"]
        heat_02_process = bss.Process.Amber(system=minimised_system, 
                                          protocol=heat_02_protocol, 
                                          name="02_heat", 
                                          work_dir=heat_02_dir, 
                                          extra_options=heat_02_options,
                                          extra_lines=namelist,
                                          exe=amber_home + "/bin/pmemd.cuda")
        heat_02_config = heat_02_dir + "/*.cfg"
        heat_02_config_file = functions.get_files(heat_02_config)[0]        

        with open(heat_02_config_file, "a") as file:
            file.write("\n")
            file.write(f"DISANG={restraints_files[0]}\n")
            file.write(f"DUMPAVE=distances.out\n")

        # heat_02_process.start()
        # heat_02_process.wait()
        # if heat_02_process.isError():
        #     print(heat_02_process.stdout())
        #     print(heat_02_process.stderr())
        #     raise bss._Exceptions.ThirdPartyError("The process exited with an error!")

        heat_02_system = heat_02_process.getSystem()

        relax_03_options = {"ioutfm": 1,
                         "cut": nonbonded_cut_off,
                         "nscm": 0,
                         "barostat": 2,
                         "ntp": 1,
                         "taup": 1.0,
                         "gamma_ln": 1.0,
                         "ntp": 1,
                         "ntb": 2, # constant pressure 
                         "irest": 1,
                         "ntx": 5,
                         "ntr": 1,
                         "nmropt": 1}
        
        relax_03_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                    runtime=runtime_ns, 
                                                    temperature_start=0.0*KELVIN, 
                                                    temperature_end=self.temperature, 
                                                    report_interval=output_frequency, 
                                                    restart_interval=output_frequency,
                                                    restraint="heavy",
                                                    force_constant=self.force_constant_0)
        
        relax_03_process = bss.Process.Amber(system=heat_02_system,  
                                          protocol=relax_03_protocol, 
                                          name="03_relax", 
                                          work_dir=relax_03_dir, 
                                          extra_options=relax_03_options,
                                          extra_lines=namelist,
                                          exe=amber_home + "/bin/pmemd.cuda")
        

        relax_03_config = relax_03_dir + "/*.cfg"
        relax_03_config_file = functions.get_files(relax_03_config)[0]

        with open(relax_03_config_file, "r") as file:
            config = file.readlines()

        new_config = [line for line in config if "tempi=" not in line and "TEMP0" not in line]
        
        with open(relax_03_config_file, "w") as file:
            file.writelines(new_config)
            file.write("\n")
            file.write(f"DISANG={restraints_files[1]}\n")
            file.write(f"DUMPAVE=distances.out\n")

        # relax_03_process.start()
        # relax_03_process.wait()
        # if relax_03_process.isError():
        #     print(relax_03_process.stdout())
        #     print(relax_03_process.stderr())
        #     raise bss._Exceptions.ThirdPartyError("The process exited with an error!")

        relax_03_system = relax_03_process.getSystem()

        lower_04_options = relax_03_options
        lower_04_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                       runtime=runtime_ns, 
                                                       temperature_start=0.0*KELVIN, 
                                                       temperature_end=self.temperature, 
                                                       report_interval=output_frequency, 
                                                       restart_interval=output_frequency,
                                                       restraint="heavy",
                                                       force_constant=self.force_constant_0 / 10)

        lower_04_process = bss.Process.Amber(system=relax_03_system,  #CHANGE
                                             protocol=lower_04_protocol, 
                                             name="04_lower", 
                                             work_dir=lower_04_dir, 
                                             extra_options=lower_04_options,
                                             extra_lines=namelist,
                                             exe=amber_home + "/bin/pmemd.cuda")
        lower_04_config = lower_04_dir + "/*.cfg"
        lower_04_config_file = functions.get_files(lower_04_config)[0]

        with open(lower_04_config_file, "r") as file:
            config = file.readlines()

        new_config = [line for line in config if "tempi=" not in line and "TEMP0" not in line]
        
        with open(lower_04_config_file, "w") as file:
            file.writelines(new_config)
            file.write("\n")
            file.write(f"DISANG={restraints_files[2]}\n")
            file.write(f"DUMPAVE=distances.out\n")

        # lower_04_process.start()
        # lower_04_process.wait()
        # if lower_04_process.isError():
        #     print(lower_04_process.stdout())
        #     print(lower_04_process.stderr())
        #     raise bss._Exceptions.ThirdPartyError("The process exited with an error!")

        lower_04_system = lower_04_process.getSystem()


        max_cycles = self.min_steps 
        n_steep_cycles = max_cycles // 33

        bb_min_05_options = {"ntmin": 1,
                                "ntpr": output_frequency,
                                "ntwx": output_frequency,
                                "ntpr": output_frequency,
                                "ntp": 0,
                                "ioutfm": 1,
                                "cut": nonbonded_cut_off, 
                                "ncyc": n_steep_cycles,
                                "ntc": 2, # constrain H bonds
                                "ntf": 2, # exlcude H bonds from force calc.
                                "ntb": 1, # constant volume
                                "ntr": 1,
                                "nmropt": 1}   
        
        bb_min_05_protocol = bss.Protocol.Minimisation(steps=max_cycles,
                                                       restraint="backbone",
                                                       force_constant=self.force_constant_0 / 10)

        bb_min_05_process = bss.Process.Amber(system=lower_04_system, 
                                                 protocol=bb_min_05_protocol, 
                                                 name="05_bb_min", 
                                                 work_dir=bb_min_05_dir, 
                                                 extra_options=bb_min_05_options,
                                                 extra_lines=namelist,
                                                 exe=amber_home + "/bin/pmemd.cuda")
        
        bb_min_05_config = bb_min_05_dir + "/*.cfg"
        bb_min_05_config_file = functions.get_files(bb_min_05_config)[0]

        with open(bb_min_05_config_file, "r") as file:
            config = file.readlines()
        new_config = [line for line in config if "ig=" not in line]
        
        with open(bb_min_05_config_file, "w") as file:
            file.writelines(new_config)
            file.write("\n")
            file.write(f"DISANG={restraints_files[3]}\n")
            file.write(f"DUMPAVE=distances.out\n")

        # bb_min_05_process.start()
        # bb_min_05_process.wait()
        # if bb_min_05_process.isError():
        #     print(bb_min_05_process.stdout())
        #     print(bb_min_05_process.stderr())
        #     raise bss._Exceptions.ThirdPartyError("The process exited with an error!")
        
        bb_min_05_system = bb_min_05_process.getSystem()        

        relax_06_options = lower_04_options
        relax_06_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                    runtime=runtime_ns, 
                                                    temperature_start=self.temperature, 
                                                    temperature_end=self.temperature, 
                                                    report_interval=output_frequency, 
                                                    restart_interval=output_frequency,
                                                    restraint="backbone",
                                                    force_constant=self.force_constant_0 / 10)

        relax_06_process = bss.Process.Amber(system=bb_min_05_system,  #CHANGE
                                          protocol=relax_06_protocol, 
                                          name="06_relax", 
                                          work_dir=relax_06_dir, 
                                          extra_options=relax_06_options,
                                          extra_lines=namelist,
                                          exe=amber_home + "/bin/pmemd.cuda")        

        relax_06_config = relax_06_dir + "/*.cfg"
        relax_06_config_file = functions.get_files(relax_06_config)[0]   

        with open(relax_06_config_file, "r") as file:
            config = file.readlines()

        new_config = [line for line in config if "tempi=" not in line and "TEMP0" not in line]
                

        with open(relax_06_config_file, "w") as file:
            file.writelines(new_config)
            file.write("\n")
            file.write(f"DISANG={restraints_files[4]}\n")
            file.write(f"DUMPAVE=distances.out\n")

        # relax_06_process.start()
        # relax_06_process.wait()
        # if relax_06_process.isError():
        #     print(relax_06_process.stdout())
        #     print(relax_06_process.stderr())
        #     raise bss._Exceptions.ThirdPartyError("The process exited with an error!")        

        relax_06_system = relax_06_process.getSystem()

        reduce_07_options = relax_06_options

        reduce_07_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                    runtime=runtime_ns, 
                                                    temperature_start=self.temperature, 
                                                    temperature_end=self.temperature, 
                                                    report_interval=output_frequency, 
                                                    restart_interval=output_frequency,
                                                    restraint="backbone",
                                                    force_constant=self.force_constant_0 / self.force_constant_0)
        
        reduce_07_process = bss.Process.Amber(system=relax_06_system,  #CHANGE
                                          protocol=reduce_07_protocol, 
                                          name="07_reduce", 
                                          work_dir=reduce_07_dir, 
                                          extra_options=reduce_07_options,
                                          extra_lines=namelist,
                                          exe=amber_home + "/bin/pmemd.cuda")

        reduce_07_config = reduce_07_dir + "/*.cfg"
        reduce_07_config_file = functions.get_files(reduce_07_config)[0]          

        with open(reduce_07_config_file, "r") as file:
            config = file.readlines()

        new_config = [line for line in config if "tempi=" not in line and "TEMP0" not in line]
                
        with open(reduce_07_config_file, "w") as file:
            file.writelines(new_config)
            file.write("\n")
            file.write(f"DISANG={restraints_files[5]}\n")
            file.write(f"DUMPAVE=distances.out\n")

        # reduce_07_process.start()
        # reduce_07_process.wait()
        # if relax_06_process.isError():
        #     print(reduce_07_process.stdout())
        #     print(reduce_07_process.stderr())
        #     raise bss._Exceptions.ThirdPartyError("The process exited with an error!")        

        reduce_07_system = reduce_07_process.getSystem()

        continue_08_options = reduce_07_options

        continue_08_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                    runtime=runtime_ns, 
                                                    temperature_start=self.temperature, 
                                                    temperature_end=self.temperature, 
                                                    report_interval=output_frequency, 
                                                    restart_interval=output_frequency,
                                                    restraint="backbone",
                                                    force_constant=(self.force_constant_0 / self.force_constant_0) / 10)
        
        continue_08_process = bss.Process.Amber(system=reduce_07_system,
                                                protocol=continue_08_protocol,
                                                name="08_continue",
                                                work_dir=continue_08_dir,
                                                extra_lines=namelist, 
                                                exe=amber_home + "/bin/pmemd.cuda")
        
        continue_08_config = continue_08_dir + "/*.cfg"
        continue_08_config_file = functions.get_files(continue_08_config)

        with open(continue_08_config_file, "r") as file:
            config = file.readlines()

        new_config = [line for line in config if "tempi=" not in line and "TEMP0" not in line]
                    
        with open(continue_08_config_file, "w") as file:
            file.writelines(new_config)
            file.write("\n")
            file.write(f"DISANG={restraints_files[6]}\n")
            file.write(f"DUMPAVE=distances.out\n")

        # continue_08_process.start()
        # continue_08_process.wait()
        # if continue_08_process.isError():
        #     print(continue_08_process.stdout())
        #     print(continue_08_process.stderr())
        #     raise bss._Exceptions.ThirdPartyError("The process exited with an error!")        

        continue_08_system = continue_08_process.getSystem()        

        relax_09_options = {"ioutfm": 1,
                            "cut": nonbonded_cut_off,
                            "tol": 0.00001,
                            "nscm": max_cycles,
                            "barostat": 2,
                            "ntp": 1,
                            "taup": 1.0,
                            "gamma_ln": 1.0,
                            "ntb": 2, 
                            "ntc": 2,
                            "ntf": 2,
                            "iwrap": 0,  
                            "irest": 1,
                            "ntx": 5, 
                            "nmropt": 1}

        relax_09_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                    runtime=runtime_ns, 
                                                    temperature=self.temperature, 
                                                    report_interval=output_frequency, 
                                                    restart_interval=output_frequency)

        relax_09_process = bss.Process.Amber(system=continue_08_system,
                                             protocol=relax_09_protocol,
                                             name="09_relax",
                                             work_dir=relax_09_dir,
                                             extra_options=relax_09_options,
                                             extra_lines=namelist,
                                             exe=amber_home + "/bin/pmemd.cuda")

        relax_09_config = relax_09_dir + "/*.cfg"
        relax_09_config_file = functions.get_files(relax_09_config)[0]        

        with open(relax_09_config_file, "a") as file:
            file.write("\n")
            file.write(f"DISANG={restraints_files[7]}\n")
            file.write(f"DUMPAVE=distances.out\n")        

        # relax_09_process.start()
        # relax_09_process.wait()
        # if relax_09_process.isError():
        #     print(relax_09_process.stdout())
        #     print(relax_09_process.stderr())
        #     raise bss._Exceptions.ThirdPartyError("The process exited with an error!")
        
        equilibrated_system = relax_09_process.getSystem()

        return equilibrated_system  


    def minimisation_0(self, ligand_name, nonbonded_cut_off=10.0):

        directory = functions.mkdir(self.equilibration_directory+f"{ligand_name}/")
        files = functions.get_files(f"{self.protein_path}/bound_{ligand_name}_solvated.*")
        solvated_system = bss.IO.readMolecules(files)
        directories = lambda step: functions.mkdir(directory + step)
        min_dir = directories("01_min")

        template_restraints_file = self.write_restraints_file_0()
        restraints_file = shutil.copy(template_restraints_file, min_dir).split("/")[-1]

        max_cycles = self.min_steps        
        n_steep_cycles = max_cycles // 20
        output_frequency = max_cycles // 20

        minimisation_options = {"ntmin": 1,
                                "ntpr": output_frequency,
                                "ntwx": output_frequency,
                                "ntpr": output_frequency,
                                "ntp": 0,
                                "ioutfm": 1,
                                "cut": nonbonded_cut_off, 
                                "ncyc": n_steep_cycles,
                                "ntc": 2, # constrain H bonds
                                "ntf": 2, # exlcude H bonds from force calc.
                                "ntb": 1, # constant volume
                                "nmropt": 1}       
        
        namelist = ["&wt TYPE='DUMPFREQ', istep1=1 /"]

        minimisation_protocol = bss.Protocol.Minimisation(steps=max_cycles,
                                                          restraint="heavy",
                                                          force_constant=self.force_constant_0)

        amber_home = os.environ["AMBERHOME"]
        minimisation_process = bss.Process.Amber(system=solvated_system, 
                                                 protocol=minimisation_protocol, 
                                                 name="01_min", 
                                                 work_dir=min_dir, 
                                                 extra_options=minimisation_options,
                                                 extra_lines=namelist,
                                                 exe=amber_home + "/bin/pmemd.cuda")
        
        min_config = min_dir + "/*.cfg"
        config_file = functions.get_files(min_config)[0]

        with open(config_file, "r") as file:
            config = file.readlines()
        new_config = [line for line in config if "ig=" not in line]
        
        with open(config_file, "w") as file:
            file.writelines(new_config)
            file.write("\n")
            file.write(f"DISANG={restraints_file}\n")
            file.write(f"DUMPAVE=distances.out\n")

        # minimisation_process.start()
        # minimisation_process.wait()
        # if minimisation_process.isError():
        #     print(minimisation_process.stdout())
        #     print(minimisation_process.stderr())
        #     raise bss._Exceptions.ThirdPartyError("The process exited with an error!")
        system = minimisation_process.getSystem()
        return system


    def qmmm_minimisation(self, ligand_name, nonbonded_cut_off=12.0):

        directory = functions.mkdir(self.equilibration_directory+f"{ligand_name}/")
        files = functions.get_files(f"{self.protein_path}/bound_{ligand_name}_solvated.*")
        solvated_system = bss.IO.readMolecules(files)
        directories = lambda step: functions.mkdir(directory + step)
        min_dir = directories("min")

        max_cycles = self.min_steps
        output_frequency = max_cycles // 20
        minimisation_options = {"ntmin": 1,
                                "ntpr": output_frequency,
                                "ntwx": output_frequency,
                                "ntpr": output_frequency,
                                "ntb": 1,
                                "ioutfm": 1,
                                "cut": nonbonded_cut_off,
                                "iwrap": 0,
                                "ifqnt": 1}
        
        qm_region = self.get_qm_region()
        qm_options = self.get_dftb3_options(qm_region)
        qm_namelist = [f"  {key}={value}" for key, value in qm_options.items()]
        qm_namelist.insert(0, "&qmmm")
        qm_namelist.append("/")

        minimisation_protocol = bss.Protocol.Minimisation(steps=max_cycles)

        minimisation_process = bss.Process.Amber(system=solvated_system, 
                                                 protocol=minimisation_protocol, 
                                                 name="min", 
                                                 work_dir=min_dir, 
                                                 extra_options=minimisation_options,
                                                 extra_lines=qm_namelist)
        
        min_config = min_dir + "/*.cfg"
        config_file = functions.get_files(min_config)[0]

        with open(config_file, "r") as file:
            config = file.readlines()
        new_config = [line for line in config if "ig=" not in line]
        
        with open(config_file, "w") as file:
            file.writelines(new_config)


    def qmmm_equilibration(self, ligand_name, nonbonded_cut_off=12.0, dt=0.001, runtime=1): # runtime in ps 

        directory = functions.mkdir(self.equilibration_directory+f"{ligand_name}/")
        files = [f"{directory}/min/min.prm7", f"{directory}/min/min.rst7"]
        minimised_system = bss.IO.readMolecules(files) 
        directories = lambda step: functions.mkdir(directory + step)
        heat_dir = directories("heat")

        max_cycles = self.min_steps
        output_frequency = max_cycles // 20
        
        equilibration_options = {"ioutfm": 1,
                                 "cut": nonbonded_cut_off,
                                 "iwrap": 0,
                                 "nscm": int(runtime / dt),
                                 "barostat": 2,
                                 "ntp": 1,
                                 "taup": 1.0,
                                 "gamma_ln": 5.0,
                                 "ntc": 1,
                                 "ntf": 1,   
                                 "ifqnt": 1}
        
        qm_region = self.get_qm_region()
        qm_options = self.get_dftb3_options(qm_region)
        qm_namelist = [f"  {key}={value}" for key, value in qm_options.items()]
        qm_namelist.insert(0, "&qmmm")
        qm_namelist.append("/")

        runtime_ns = functions.convert_to_units(runtime / 1000, NANOSECOND)
        equilibration_protocol = bss.Protocol.Equilibration(timestep=dt*PICOSECOND, 
                                                            runtime=runtime_ns, 
                                                            temperature_start=0.0*KELVIN, 
                                                            temperature_end=self.temperature, 
                                                            report_interval=output_frequency, 
                                                            restart_interval=output_frequency)

        equilibration_process = bss.Process.Amber(system=minimised_system, 
                                                  protocol=equilibration_protocol, 
                                                  name="heat", 
                                                  work_dir=heat_dir, 
                                                  extra_options=equilibration_options,
                                                  extra_lines=qm_namelist)
        
        heat_config = heat_dir + "/*.cfg"
        config_file = functions.get_files(heat_config)[0]

        with open(config_file, "r") as file:
            config = file.readlines()

        new_config = [line for line in config if "tempi=" not in line and "TEMP0" not in line]
        
        with open(config_file, "w") as file:
            file.writelines(new_config)
        

    def qmmm_production(self, ligand_name, nonbonded_cut_off=12.0, dt=0.001, runtime=10): # runtime in ps 

        directory = functions.mkdir(self.output_directory+f"{ligand_name}/")
        files = [f"{directory}/heat/heat.prm7", f"{directory}/heat/heat.rst7"]
        heated_system = bss.IO.readMolecules(files) 
      
        max_cycles = self.min_steps
        output_frequency = max_cycles // 20
        ncsm = runtime // 10
        production_options = {"ioutfm": 1,
                              "cut": nonbonded_cut_off,
                              "iwrap": 0,
                              "nscm": ncsm,
                              "barostat": 2,
                              "ntp": 1,
                              "taup": 1.0,
                              "gamma_ln": 5.0,
                              "ntc": 1,
                              "ntf": 1,   
                              "ifqnt": 1}
        
        qm_region = self.get_qm_region()
        qm_options = self.get_dftb3_options(qm_region)
        qm_namelist = [f"  {key}={value}" for key, value in qm_options.items()]
        qm_namelist.insert(0, "&qmmm")
        qm_namelist.append("/")

        runtime_ns = functions.convert_to_units(runtime / 1000, NANOSECOND)
        production_protocol = bss.Protocol.Production(timestep=dt*PICOSECOND, 
                                                      runtime=runtime_ns, 
                                                      temperature=self.temperature, 
                                                      restart=True,
                                                      report_interval=output_frequency, 
                                                      restart_interval=output_frequency)

        production_process = bss.Process.Amber(system=heated_system, 
                                               protocol=production_protocol, 
                                               name="qmmm", 
                                               work_dir=directory, 
                                               extra_options=production_options,
                                               extra_lines=qm_namelist)
        
        production_config = directory + "/*.cfg"
        config_file = functions.get_files(production_config)[0]

        with open(config_file, "r") as file:
            config = file.readlines()

        new_config = [line for line in config if "tempi=" not in line and "TEMP0" not in line and "pres0" not in line]
        
        with open(config_file, "w") as file:
            file.writelines(new_config)






    
def get_side_chain_selection(universe, residue):
    """
    Take a residue and return the atom indices of the side chain of that residue.

    Parameters:
    -----------
    universe: mda.Universe
    residue: mda.residue

    Return:
    -------
    str 
        atom selection in the format "{first_atom}-{last_atom}"
    """
    n_terminus = "name N or name H"
    alpha_carbon = "name CA or name HA"
    c_terminus = "name C or name O"
    atoms_in_residue = universe.select_atoms(f"resid {residue.resid}")
    qm_region_for_residue = list(atoms_in_residue.select_atoms(f"not ({n_terminus} or {alpha_carbon} or {c_terminus})").ids)
    first_atom = qm_region_for_residue[0]
    last_atom = qm_region_for_residue[-1]
    return f"{first_atom}-{last_atom}"





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
    
    arguments = parser.parse_args()
    
    protocol = functions.input_to_dict(arguments.protocol_file)

    keys = list(protocol.keys())
    if "metal" in keys:
        metal = True
    else:
        metal = False
    
    if metal:
        meze = Meze(prepared=True,
                    protein_file=protocol["protein input file"],
                    cut_off=protocol["cutoff"],
                    force_constant_0=protocol["force constant"],
                    workdir=protocol["project directory"],
                    equilibration_path=protocol["equilibration directory"],
                    afe_input_path=protocol["afe input directory"],
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
                    pressure=protocol["pressure"])
        
    elif not metal:
      
        meze = sofra.Sofra(prepared=True,
                               equilibration_path=protocol["equilibration directory"],
                               outputs=protocol["outputs"],
                               workdir=protocol["project directory"],
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
                               pressure=protocol["pressure"])
          
    ligand_a, ligand_b = functions.separate(arguments.transformation)
    
    equilibrated_network = meze.get_equilibrated(ligand_a, ligand_b)

    equilibrated_network.prepare_afe(ligand_a, ligand_b, extra_edges=arguments.extra_transformations_file) 

if __name__ == "__main__":
    main()
