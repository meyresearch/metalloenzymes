import meze


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





def run_qmmm(ligand_name):

    qmmm_minimisation(ligand_name)
    qmmm_equilibration(ligand_name)
    qmmm_production(ligand_name)
