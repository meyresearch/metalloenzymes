from meze.network import Network
import os
import BioSimSpace as bss
import numpy as np
import functions
from definitions import PICOSECOND, KELVIN, NANOSECOND
import csv
import MDAnalysis as mda
import MDAnalysis.analysis.distances 
import shutil


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



class Meze(Network):

    def __init__(self, protein_file, prepared=False, metal="ZN", cut_off=2.6, force_constant_0=100,
                 water_file=None, workdir=os.getcwd(), is_qm=False, qmmm_inputs=None, equilibration_path=None, output=None,
                 ligand_path=os.getcwd()+"/inputs/ligands/", ligand_charge=0, ligand_ff="gaff2",
                 group_name=None, protein_path=os.getcwd()+"/inputs/protein/", water_model="tip3p", protein_ff="ff14SB", 
                 engine=None, sampling_time=10, box_edges=20, box_shape="cubic", min_steps=1000, short_nvt=0, nvt=1, npt=1, 
                 min_dt=None, min_tol=None, repeats=0, temperature=300, pressure=1, threshold=0.4, n_normal=11, n_difficult=17):
        
        super().__init__(prepared=prepared, workdir=workdir, ligand_path=ligand_path, group_name=group_name, protein_file=protein_file, protein_path=protein_path, 
                         water_model=water_model, ligand_ff=ligand_ff, protein_ff=protein_ff, ligand_charge=ligand_charge, equilibration_path=equilibration_path, outputs=output,
                         engine=engine, sampling_time=sampling_time, box_edges=box_edges, box_shape=box_shape, min_steps=min_steps, short_nvt=short_nvt, nvt=nvt, npt=npt, 
                         min_dt=min_dt, min_tol=min_tol, repeats=repeats, temperature=temperature, pressure=pressure, threshold=threshold, n_normal=n_normal, n_difficult=n_difficult)
        
        self.md_time = functions.convert_to_units(sampling_time, PICOSECOND)
        self.universe = mda.Universe(self.protein_file, format="pdb")
        self.force_constant_0 = functions.check_float(force_constant_0)
        self.is_qm = is_qm
        if qmmm_inputs and self.is_qm:
            self.input_directory = functions.path_exists(qmmm_inputs)
        elif not qmmm_inputs and self.is_qm:
            os.rmdir(self.afe_input_directory)
            self.input_directory = self.create_directory(f"/qmmm_input_files/")
        # To be depracated
        # elif not self.is_qm:
        #     os.rmdir(self.afe_input_directory)
        #     self.input_directory = self.create_directory(f"/md_input_files/")

        if output:
            self.output_directory = functions.path_exists(output)
        else:
            self.output_directory = self.create_directory(f"/outputs/")
            
        if water_file:
            self.xtal_water = bss.IO.readMolecules(water_file)

        if metal.upper() == "ZN":
            self.metal_resname = metal
            self.metals = self.universe.select_atoms(f"resname {self.metal_resname}")
            self.metal_resids = self.metals.resids
            self.metal_atomids = self.metals.atoms.ids
        else:
            print(f"Metal {metal} is not supported yet.")
        
        self.cut_off = functions.check_positive(functions.check_float(cut_off))


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
        self.log_directory = self.create_directory("/logs/")
        self.protocol_file = self.create_qm_protocol_file()
        self.prepared = True
        return self


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
        strip = self.output_directories[0].split("/")[-2]
        path_to_outputs = self.output_directories[0].replace(strip, "")
        protocol = [f"group name = {self.group_name}",
                    f"metal = {self.metal_resname}",
                    f"cutoff = {self.cut_off}",
                    f"force constant = {self.force_constant_0}",
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
                    f"outputs = {path_to_outputs}",
                    f"repeats = {self.n_repeats}",
                    f"network file = {self.network_file}",
                    f"project directory = {self.workding_directory}",
                    f"equilibration directory = {self.equilibration_directory}",
                    f"ligand directory = {self.ligand_path}",
                    f"protein directory = {self.protein_path}",
                    f"log directory = {self.log_directory}",
                    f"afe input directory = {self.afe_input_directory}"]

        protocol_file = self.afe_input_directory + "/protocol.dat"

        with open(protocol_file, "w") as file:
            writer = csv.writer(file)
            for protocol_line in protocol:
                writer.writerow([protocol_line])
        return protocol_file
            

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
    
 
    def get_metal_ligands(self):
        """
        Create a dictionary of metal ion(s) and their coordinating ligand(s)

        Parameters:
        -----------

        Return:
        -------
        metal_ligands: dict
            dictionary where keys are metal ions (differentiated by subscripts) and values are the atom groups containing coordinating ligands for that ion
        """
        metal_ligands = {}
        for i in range(len(self.metal_resids)):
            ligands = self.universe.select_atoms(f"not (name {self.metal_resname} or element H) and sphzone {self.cut_off} (resid {self.metal_resids[i]})")
            key = self.metal_atomids[i] 
            metal_ligands[key] = ligands
        return metal_ligands


    def model_0(self, ligand_name):

        if self.is_qm:
            self.qmmm_minimisation(ligand_name)
            self.qmmm_equilibration(ligand_name)
            self.qmmm_production(ligand_name)
        # To be depracated
        # elif not self.is_qm:
        #     minimised_system = self.minimisation_0(ligand_name)
        #     equilibrated_system = self.equilibration_0(ligand_name, minimised_system)
        #     self.production_0(ligand_name, equilibrated_system)


    def write_restraints_file_0(self):
        """
        Write an Amber-compatible restraint file for model 0. 

        Parameters:
        -----------

        Return:
        -------
        output_file: str
            amber-format restraints file for model 0
        """
        metal_ligands = self.get_metal_ligands()
        protein = self.universe.select_atoms("protein")
        atom_ids = []
        r1, r2, r3, r4 = [], [], [], []
        for i in range(len(self.metal_atomids)):
            metal_id = self.metal_atomids[i]
            for ligating_atom in metal_ligands[metal_id]:
                if ligating_atom in protein or ligating_atom.resname == "WAT":
                    atom_ids.append(f"iat={metal_id},{ligating_atom.id}")
                    atom_group_1 = self.universe.select_atoms(f"resid {self.metal_resids[i]}")
                    atom_group_2 = self.universe.select_atoms(f"resid {ligating_atom.resid} and name {ligating_atom.name}")
                    distance = MDAnalysis.analysis.distances.dist(atom_group_1, atom_group_2)[-1][0]
                    linear_response_region_lower_bound = np.round(distance - 1.0, decimals=2)
                    flat_region_lower_bound = np.round(distance - 0.5, decimals=2)
                    flat_region_upper_bound = np.round(distance + 0.5, decimals=2)
                    linear_response_region_upper_bound = np.round(distance + 1.0, decimals=2)
                    r1.append("r1="+str(linear_response_region_lower_bound))
                    r2.append("r2="+str(flat_region_lower_bound))
                    r3.append("r3="+str(flat_region_upper_bound))
                    r4.append("r4="+str(linear_response_region_upper_bound))
        rk2 = ["rk2="+str(self.force_constant_0)] * len(atom_ids)
        rk3 = ["rk3="+str(self.force_constant_0)] * len(atom_ids)
        output_file = self.input_directory + "restraints_0.RST"
        with open(output_file, "w") as file:
            file.write(f"# Harmonic bond restraints between {self.metal_resname} and coordinating protein residues\n")
            for i in range(len(atom_ids)):
                line = f"&rst {atom_ids[i]}, {r1[i]}, {r2[i]}, {r3[i]}, {r4[i]}, {rk2[i]}, {rk3[i]},/\n"
                file.write(line)
        return output_file


    # Potentially to be depracated:
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
                        f"project directory = {self.workding_directory}",
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
        ncsm = runtime // 10 # change this
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
