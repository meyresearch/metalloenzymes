"""
Minimise and equilibrate bound and unbound stages.
"""
import shutil
from definitions import FEMTOSECOND, PICOSECOND, KELVIN, ATM
import functions
import BioSimSpace as bss
import argparse
import os 
import meze
import warnings


class coldMeze(meze.Meze):

    def __init__(self, group_name, ligand_name, outputs, input_protein_file,
                 protein_directory, ligand_directory, log_directory, equilibration_directory, afe_input_directory, 
                 min_steps, short_nvt, nvt, npt, min_dt, min_tol, temperature, pressure, short_timestep=0.5, 
                 is_md=False, is_metal=True, restraints=True, prepared=True, 
                 force_constant_0=100, restraint_weight=10, restart_write_steps=100, coordinate_write_steps=500):

        self.is_metal = is_metal
        self.prepared = prepared
        if self.is_metal:
            super().__init__(protein_file=input_protein_file, prepared=prepared, group_name=group_name, is_md=is_md, log_directory=log_directory,
                             equilibration_path=equilibration_directory, afe_input_path=afe_input_directory, outputs=outputs,
                             protein_path=protein_directory, ligand_path=ligand_directory, force_constant_0=force_constant_0)
        else:

            self.protein_file = input_protein_file
            super(meze.Meze, self).__init__(protein_file=input_protein_file, prepared=prepared, group_name=group_name, is_md=is_md, log_directory=log_directory,
                                            equilibration_path=equilibration_directory, afe_input_path=afe_input_directory, outputs=outputs,
                                            protein_path=protein_directory, ligand_path=ligand_directory,)
            restraints = False

        self.restraints = restraints
            
        self.ligand_name = ligand_name

        self.min_steps = min_steps
        self.min_dt = min_dt
        self.min_tol = min_tol
        if not prepared:
            self.temperature = functions.convert_to_units(temperature, KELVIN)
            self.pressure = functions.convert_to_units(pressure, ATM)
            self.short_nvt = functions.convert_to_units(short_nvt, PICOSECOND)
            self.nvt = functions.convert_to_units(nvt, PICOSECOND)
            self.npt = functions.convert_to_units(npt, PICOSECOND)
            self.short_timestep = functions.convert_to_units(short_timestep, FEMTOSECOND)
        else:
            self.temperature = temperature
            self.pressure = pressure
            self.short_nvt = short_nvt
            self.nvt = nvt
            self.npt = npt
            self.short_timestep = short_timestep       
        self.restraint_weight = functions.check_float(restraint_weight)
        self.restart_write_steps = functions.check_positive(functions.check_int(restart_write_steps))
        self.coordinate_write_steps = functions.check_positive(functions.check_int(coordinate_write_steps))


    def run(self, system, protocol, name, working_directory, configuration={}, checkpoint=None, restraints_file=None):
        """
        Run a minimisation or equilibration process 
        Adapted from https://tinyurl.com/BSSligprep

        Parameters:
        -----------
        system: bss.System
            run system
        protocol: bss.Protocol 
            minimisation or equilibration
        name: str 
            process name for saving process output
        working_directory: str
            save output into this directory
        configuration: dict 
            extra options to pass to the process, overwrites defaults
        checkpoint: str
            path to a checkpoint file from a previous run; corresponds to the -t flag for gmx grompp

        Return:
        -------
        system: bss.System
            equilibrated or minimised system
        """
        if checkpoint:
            configuration["irest"] = 1
            configuration["ntx"] = 5 


        if self.is_metal and restraints_file: # restraints for bound with metalloenzymes
            try:
                restraints_file = shutil.copy(restraints_file, working_directory).split("/")[-1]
            except shutil.SameFileError as e:
                print(e)
                restraints_file = restraints_file.split("/")[-1]
                pass
            namelist = ["&wt TYPE='DUMPFREQ', istep1=1 /"]
            amber_path = os.environ["AMBERHOME"] + "/bin/pmemd.cuda"

            process = bss.Process.Amber(system=system, protocol=protocol, name=name, work_dir=working_directory, extra_options=configuration, extra_lines=namelist, exe=amber_path)             
            config = working_directory + "/*.cfg"
            config_file = functions.get_files(config)[0]



            with open(config_file, "a") as file:
                file.write("\n")
                file.write(f"DISANG={restraints_file}\n")
                file.write(f"DUMPAVE=distances.out\n")

        elif not restraints_file: 
            amber_path = os.environ["AMBERHOME"] + "/bin/pmemd.cuda"
            process = bss.Process.Amber(system=system, protocol=protocol, name=name, work_dir=working_directory, extra_options=configuration, exe=amber_path)

        with open(config_file, "r") as file:
            config_lines = file.readlines()

        config_lines = list(map(lambda key_word: key_word.replace('   restraintmask="@N,CA,C,O",\n', '   restraintmask="@N,CA,C",\n'), config_lines))

        with open(config_file, "w") as file:
            file.writelines(config_lines)


        process.start()
        process.wait()
        if process.isError():

            amberhome = os.environ["AMBERHOME"]
            try:
                sander_mpi = functions.file_exists(f"{amberhome}/bin/sander.MPI")

                readme = working_directory + "/" + "README.txt"
                with open(readme, "r") as file:
                    all_lines = file.readlines()
                    for line in all_lines:
                        if "#" not in line:
                            cuda_command = line
                amber_command = " ".join(cuda_command.split()[1:])
                sander_run_command = sander_mpi + " " + amber_command
                code_dir = os.getcwd()
                os.chdir(working_directory) 
                os.system(f"mpirun {sander_run_command}")
                os.chdir(code_dir)
                with open(readme, "a") as file:
                    file.write(f"# Process {name} was run with the below command for {working_directory}:\n")
                    file.write(sander_run_command + "\n")
                system = bss.IO.readMolecules([working_directory + "/" + name + ".prm7", working_directory + "/" + name + ".rst7"])

            except argparse.ArgumentTypeError:

                print("No sander.MPI detected, running with regular sander.")
                sander_mpi = functions.file_exists(f"{amberhome}/bin/sander")
                readme = working_directory + "README.txt"
                with open(readme, "r") as file:
                    all_lines = file.readlines()
                    for line in all_lines:
                        if "#" not in line:
                            cuda_command = line
                amber_command = " ".join(cuda_command.split()[1:])
                sander_run_command = sander_mpi + " " + amber_command
                code_dir = os.getcwd()
                os.chdir(working_directory) 
                os.system(sander_run_command)
                os.chdir(code_dir)
                system = bss.IO.readMolecules([working_directory + "/" + name + ".prm7", working_directory + "/" + name + ".rst7"])                    
                with open(readme, "a") as file:
                    file.write("\n")
                    file.write(f"# Step {name} was rerun using sander.MPI for {working_directory}:\n")
                    file.write(f"{sander_run_command}\n")
        else:
            system = process.getSystem()
        return system


    def minimise(self, system, working_directory, process_name="min", configuration={}, restraints_file=None, position_restraints=None, checkpoint=None):
        """
        Minimise the system using Gromacs

        Parameters:
        -----------
        system: bss.System
            system to be minimised
        working_directory: str
            current working dir
        configuration

        Return:
        -------
        minimsed_system: bss.System
            minimised system
        """

        protocol = bss.Protocol.Minimisation(steps=self.min_steps, restraint=position_restraints, force_constant=self.restraint_weight)

        minimised_system = self.run(system, protocol, process_name, working_directory, configuration=configuration, restraints_file=restraints_file, checkpoint=checkpoint)
        return minimised_system   


    def heat(self, system, process_name, working_directory, time, timestep=2, start_t=300, end_t=300, temperature=None, pressure=None, configuration={}, position_restraints=None, checkpoint=None, restraints_file=None):
        """
        Run NVT or NPT equilibration

        Parameters:
        -----------
        system: bss.System
            system to be equilibrated
        name: str
            name of equilibration process
        working_directory: str
            path to the equilibration step directory
        time: bss.Units.Time
            runtime for equilibration

        Return:
        -------
        equilibrated_system: bss.System
            
        """
        if type(timestep) != bss.Types._time.Time:
            timestep = functions.convert_to_units(timestep, FEMTOSECOND)
        protocol = bss.Protocol.Equilibration(timestep=timestep, 
                                              runtime=time,
                                              temperature_start=start_t,
                                              temperature_end=end_t,
                                              temperature=temperature,
                                              pressure=pressure,
                                              restraint=position_restraints,
                                              force_constant=self.restraint_weight,
                                              restart_interval=self.restart_write_steps,
                                              report_interval=self.coordinate_write_steps)
        equilibrated_system = self.run(system, protocol, process_name, working_directory, configuration, checkpoint, restraints_file) 
        return equilibrated_system


    def heat_unbound(self):
        """
        Perform minimisation and NVT and NPT equilibrations on ligand

        Parameters:
        -----------

        Return:
        -------

        """

        directory = functions.mkdir(self.equilibration_directory + f"/unbound/{self.ligand_name}/")
        files = functions.get_files(f"{self.ligand_path}/{self.ligand_name}_solvated.*")
        solvated_ligand = bss.IO.readMolecules(files)
        print(f"Equilibrating unbound ligand {self.ligand_name}")
        directories = lambda step: functions.mkdir(directory+step)
        min_directory = directories("min")
        r_nvt_directory = directories("r_nvt")
        nvt_directory = directories("nvt")
        r_npt_directory = directories("r_npt")
        npt_directory = directories("npt")

        minimised_ligand = self.minimise(system=solvated_ligand, working_directory=min_directory)
        start_temp = functions.convert_to_units(0, KELVIN)
        restrained_nvt = self.heat(system=minimised_ligand,
                                   process_name="r_nvt",
                                   working_directory=r_nvt_directory,
                                   time=self.short_nvt,
                                   start_t=start_temp, end_t=self.temperature,
                                   timestep=self.short_timestep,
                                   position_restraints="all")
        nvt = self.heat(system=restrained_nvt,
                        process_name="nvt",
                        working_directory=nvt_directory,
                        time=self.nvt,
                        temperature=self.temperature,
                        checkpoint=r_nvt_directory + "/r_nvt.cpt")
        restrained_npt = self.heat(system=nvt,
                                   process_name="r_npt",
                                   working_directory=r_npt_directory,
                                   time=self.npt,
                                   pressure=self.pressure,
                                   temperature=self.temperature,
                                   position_restraints="heavy",
                                   checkpoint=nvt_directory + "/nvt.cpt")
        equilibrated_molecule = self.heat(system=restrained_npt,
                                          process_name="npt",
                                          working_directory=npt_directory,
                                          time=self.npt,
                                          pressure=self.pressure,
                                          temperature=self.temperature,
                                          checkpoint=r_npt_directory + "/r_npt.cpt")
        unbound_savename = npt_directory + f"/{self.ligand_name}"
        bss.IO.saveMolecules(filebase=unbound_savename, system=equilibrated_molecule, fileformat=["PRM7", "RST7"])        

    
    def heat_bound(self):
        """
        Perform minimisation and NVT and NPT equilibrations on bound ligand 

        Parameters:
        -----------

        Return:
        -------
        
        """ 
            
        directory = functions.mkdir(self.equilibration_directory+f"/bound/{self.ligand_name}/")
        filename = f"bound_{self.ligand_name}_solvated"
        files = functions.get_files(f"{self.protein_path}/{filename}" + ".*")
        solvated_system = bss.IO.readMolecules(files)
        
        self.set_universe(self.protein_path + "/" + filename)

        directories = lambda step: functions.mkdir(directory+step)
        min_dir = directories("min")
        r_nvt_dir = directories("r_nvt")
        bb_r_nvt_dir = directories("bb_r_nvt")
        nvt_dir = directories("nvt")
        r_npt_dir = directories("r_npt")
        npt_dir = directories("npt")     
        start_temp = functions.convert_to_units(0, KELVIN)
        print(f"Equilibrating bound ligand {self.ligand_name}")
        restraints_file = None
        configuration = {}
        if self.is_metal and self.restraints:
            configuration = {"nmropt": 1}
            restraints_file = self.write_restraints_file_0(workdir=directory)
        elif not self.is_metal:
            configuration = {"emstep": self.min_dt, "emtol": self.min_tol}
        minimised_system = self.minimise(system=solvated_system, working_directory=min_dir, configuration=configuration, restraints_file=restraints_file)

        restrained_nvt = self.heat(system=minimised_system,
                                   working_directory=r_nvt_dir,
                                   process_name="r_nvt",
                                   time=self.short_nvt,
                                   start_t=start_temp, end_t=self.temperature,
                                   position_restraints="all",
                                   timestep=self.short_timestep, 
                                   restraints_file=restraints_file,
                                   configuration=configuration) 
        backbone_restrained_nvt = self.heat(system=restrained_nvt,
                                            process_name="bb_r_nvt",
                                            working_directory=bb_r_nvt_dir,
                                            time=self.nvt,
                                            temperature=self.temperature,
                                            position_restraints="backbone",
                                            checkpoint=r_nvt_dir + "/r_nvt.cpt",
                                            restraints_file=restraints_file,
                                            configuration=configuration)
        nvt = self.heat(system=backbone_restrained_nvt,
                        process_name="nvt",
                        working_directory=nvt_dir,
                        time=self.nvt,
                        temperature=self.temperature,
                        checkpoint=bb_r_nvt_dir + "/bb_r_nvt.cpt", 
                        restraints_file=restraints_file,
                        configuration=configuration)
        restrained_npt = self.heat(system=nvt,
                                   process_name="r_npt",
                                   working_directory=r_npt_dir,
                                   time=self.npt,
                                   pressure=self.pressure,
                                   temperature=self.temperature,
                                   position_restraints="heavy",
                                   checkpoint=nvt_dir + "/nvt.cpt", 
                                   restraints_file=restraints_file,
                                   configuration=configuration)
        equilibrated_protein = self.heat(system=restrained_npt,
                                         process_name="npt",
                                         working_directory=npt_dir,
                                         time=self.npt,
                                         pressure=self.pressure,
                                         temperature=self.temperature,
                                         checkpoint=r_npt_dir + "/r_npt.cpt", 
                                         restraints_file=restraints_file,
                                         configuration=configuration)
        bound_savename = npt_dir + f"/bound_{self.ligand_name}"
        bss.IO.saveMolecules(filebase=bound_savename, system=equilibrated_protein, fileformat=["PRM7", "RST7"])     
    

    def heat_md(self):
        # heat a bound system for md simulation
        directory = functions.mkdir(self.equilibration_directory+f"/{self.ligand_name}/")
        filename = f"bound_{self.ligand_name}_solvated"
        files = functions.get_files(f"{self.protein_path}/{filename}" + ".*")
        solvated_system = bss.IO.readMolecules(files)
        
        if self.is_metal:
            self.set_universe(self.protein_path + "/" + filename)
        
        start_temp = functions.convert_to_units(0, KELVIN)
        
        directories = lambda step: functions.mkdir(directory+step)

        min_dir = directories("01_min")
        heat_02_dir = directories("02_heat")
        relax_03_dir = directories("03_relax")
        lower_04_dir = directories("04_lower")
        relax_05_dir = directories("05_relax")
        reduce_06_dir = directories("06_reduce")
        continue_07_dir = directories("07_continue")
        relax_08_dir = directories("08_relax")   

        restraints_file = None
        configuration = {}
        if self.is_metal:
            configuration = {"nmropt": 1}
            restraints_file = self.write_restraints_file_0(workdir=directory)

        minimised_system = self.minimise(system=solvated_system, 
                                         process_name="01_min", 
                                         working_directory=min_dir, 
                                         position_restraints="heavy", 
                                         configuration=configuration, 
                                         restraints_file=restraints_file)
        
        heat_02_system = self.heat(system=minimised_system,
                                   working_directory=heat_02_dir,
                                   process_name="02_heat",
                                   time=self.nvt,
                                   start_t=start_temp, end_t=self.temperature,
                                   position_restraints="heavy",
                                   timestep=self.short_timestep,
                                   restraints_file=restraints_file,
                                   configuration=configuration) 
        
        relax_03_system = self.heat(system=heat_02_system,
                                    working_directory=relax_03_dir,
                                    process_name="03_relax",
                                    time=self.nvt,
                                    temperature=self.temperature,
                                    position_restraints="heavy",
                                    timestep=self.short_timestep,
                                    restraints_file=restraints_file,
                                    configuration=configuration,
                                    checkpoint=heat_02_dir + "/02_heat")
        
        self.restraint_weight = self.restraint_weight / 10

        lower_04_system = self.heat(system=relax_03_system,
                                    working_directory=lower_04_dir,
                                    process_name="04_lower",
                                    time=self.nvt,
                                    temperature=self.temperature,
                                    position_restraints="heavy",
                                    timestep=self.short_timestep,
                                    restraints_file=restraints_file,
                                    configuration=configuration,
                                    checkpoint=relax_03_dir + "/03_relax")
        # Debugging
        # lower_04_system = bss.IO.readMolecules(["/home/jguven/projects/alchemistry/vim2/md/equilibration/ligand_10/04_lower/04_lower.prm7",
        #                                         "/home/jguven/projects/alchemistry/vim2/md/equilibration/ligand_10/04_lower/04_lower.rst7"])
        
        relax_05_system = self.heat(system=lower_04_system,
                                    working_directory=relax_05_dir,
                                    process_name="05_relax",
                                    time=self.nvt,
                                    temperature=self.temperature,
                                    position_restraints="backbone",
                                    timestep=self.short_timestep,
                                    restraints_file=restraints_file,
                                    configuration=configuration,
                                    checkpoint=lower_04_dir + "/04_lower")
        
        self.restraint_weight = self.restraint_weight / 10

        reduce_06_system = self.heat(system=relax_05_system,
                                     working_directory=reduce_06_dir,
                                     process_name="06_reduce",
                                     time=self.nvt,
                                     temperature=self.temperature,
                                     position_restraints="backbone",
                                     timestep=self.short_timestep,
                                     restraints_file=restraints_file,
                                     configuration=configuration,
                                     checkpoint=relax_05_dir + "/05_relax")
        
        self.restraint_weight = self.restraint_weight / 10

        continue_07_system = self.heat(system=reduce_06_system,
                                       working_directory=continue_07_dir,
                                       process_name="07_continue",
                                       time=self.nvt,
                                       temperature=self.temperature,
                                       position_restraints="backbone",
                                       timestep=self.short_timestep,
                                       restraints_file=restraints_file,
                                       configuration=configuration,
                                       checkpoint=reduce_06_dir + "/06_reduce")

        # Debugging
        # continue_07_system = bss.IO.readMolecules([continue_07_dir + "/07_continue.top", continue_07_dir + "/07_continue.gro"])
        
        relax_08_system = self.heat(system=continue_07_system,
                                    working_directory=relax_08_dir,
                                    process_name="08_relax",
                                    time=self.nvt,
                                    temperature=self.temperature,
                                    timestep=self.short_timestep,
                                    restraints_file=restraints_file,
                                    configuration=configuration,
                                    checkpoint=continue_07_dir + "/07_continue")
        
        savename = relax_08_dir + f"/bound_{self.ligand_name}"
        bss.IO.saveMolecules(filebase=savename, system=relax_08_system, fileformat=["PRM7", "RST7"])
        return self, relax_08_system
        

def main():

    parser = argparse.ArgumentParser(description="minimisation and equilibration for meze workflow")

    parser.add_argument("ligand_name",
                        help="ligand name",
                        type=str)
    
    parser.add_argument("protocol_file",
                        help="protocol file containing equilibration options",
                        type=str,
                        default=os.getcwd() + "/afe/protocol.dat")

    parser.add_argument("--no-restraints",
                        dest="no_restraints",
                        help="do not apply harmonic restraints on metal-coordinating residues",
                        action="store_true")
    
    arguments = parser.parse_args()
    protocol = functions.input_to_dict(file=arguments.protocol_file)

    keys = list(protocol.keys())
    if "metal" in keys:
        metal = True
    else:
        metal = False

    if arguments.no_restraints:
        apply_restraints = False
    else:
        apply_restraints = True

    cold_meze = coldMeze(is_metal=metal,
                         restraints=apply_restraints,
                         group_name=protocol["group name"],
                         ligand_name=arguments.ligand_name,
                         afe_input_directory=protocol["afe input directory"],
                         outputs=protocol["outputs"],
                         log_directory=protocol["log directory"],
                         equilibration_directory=protocol["equilibration directory"],
                         input_protein_file=protocol["prepared protein file"],
                         protein_directory=protocol["protein directory"],
                         ligand_directory=protocol["ligand directory"],
                         min_steps=protocol["minimisation steps"],
                         short_nvt=protocol["short nvt"],
                         nvt=protocol["nvt"],
                         npt=protocol["npt"],
                         min_dt=protocol["minimisation stepsize"],
                         min_tol=protocol["minimisation tolerance"],
                         temperature=protocol["temperature"],
                         pressure=protocol["pressure"])
    
    cold_meze.heat_bound()  
    cold_meze.heat_unbound()


if __name__ == "__main__":
    main()