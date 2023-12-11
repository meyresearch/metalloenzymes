"""
Minimise and equilibrate bound and unbound stages.
"""
import shutil
from definitions import FEMTOSECOND, PICOSECOND, KELVIN, ATM
import functions
import BioSimSpace as bss
import argparse
import os 
from Meze import Meze


class coldMeze(Meze):
    def __init__(self, group_name, ligand_name, equilibration_directory, input_protein_file, protein_directory, ligand_directory, min_steps, short_nvt, nvt, npt, min_dt, min_tol, temperature, pressure, short_timestep=0.5, is_metal=True):
        
        self.is_metal = is_metal
        if self.is_metal:
            super().__init__(protein_file=input_protein_file, prepared=True, group_name=group_name)
        self.ligand_name = ligand_name
        self.equilibration_directory = equilibration_directory
        self.ligand_path = functions.path_exists(ligand_directory)
        self.protein_path = functions.path_exists(protein_directory)
        self.short_nvt = functions.convert_to_units(short_nvt, PICOSECOND)
        self.nvt = functions.convert_to_units(nvt, PICOSECOND)
        self.npt = functions.convert_to_units(npt, PICOSECOND)
        self.short_timestep = functions.convert_to_units(short_timestep, FEMTOSECOND)
        self.min_steps = min_steps
        self.min_dt = min_dt
        self.min_tol = min_tol
        self.temperature = functions.convert_to_units(temperature, KELVIN)
        self.pressure = functions.convert_to_units(pressure, ATM)
        

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
            extra options to pass to the process
        checkpoint: str
            path to a checkpoint file from a previous run; corresponds to the -t flag for gmx grompp

        Return:
        -------
        system: bss.System
            equilibrated or minimised system
        """
        if self.is_metal and restraints_file: # restraints for bound
            restraints_file = shutil.copy(restraints_file, working_directory).split("/")[-1]
            namelist = ["&wt TYPE='DUMPFREQ', istep1=1 /"]
            amber_path = os.environ["AMBERHOME"] + "/bin/pmemd.cuda"
            process = bss.Process.Amber(system=system, protocol=protocol, name=name, work_dir=working_directory, extra_options=configuration, extra_lines=namelist, exe=amber_path)
            config = working_directory + "/*.cfg"
            config_file = functions.read_files(config)[0]
            with open(config_file, "a") as file:
                file.write("\n")
                file.write(f"DISANG={restraints_file}\n")
                file.write(f"DUMPAVE=distances.out\n")
        elif self.is_metal and not restraints_file: # unbound doesn't have restraints
            amber_path = os.environ["AMBERHOME"] + "/bin/pmemd.cuda"
            process = bss.Process.Amber(system=system, protocol=protocol, name=name, work_dir=working_directory, extra_options=configuration, exe=amber_path)
        else:
            process = bss.Process.Gromacs(system, protocol, name=name, work_dir=working_directory, extra_options=configuration, checkpoint_file=checkpoint)
            process.setArg("-ntmpi", 1)
        process.start()
        process.wait()
        if process.isError():
            print(process.stdout())
            print(process.stderr())
            raise bss._Exceptions.ThirdPartyError("The process exited with an error!")
        system = process.getSystem()
        return system


    def minimise(self, system, working_directory, configuration={}, restraints_file=None):
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
        protocol = bss.Protocol.Minimisation(steps=self.min_steps)
        minimised_system = self.run(system, protocol, "min", working_directory, configuration=configuration, restraints_file=restraints_file)
        return minimised_system   


    def heat(self, system, process_name, working_directory, time, timestep=2, start_t=300, end_t=300, temperature=None, pressure=None, configuration={}, restraints=None, checkpoint=None, restraints_file=None):
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
                                              restraint=restraints)
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
        files = functions.read_files(f"{self.ligand_path}/{self.ligand_name}_solvated.*")
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
                                   restraints="all")
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
                                   restraints="heavy",
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
        files = functions.read_files(f"{self.protein_path}/{filename}" + ".*")
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
        if self.is_metal:
            configuration = {"nmropt": 1}
            restraints_file = self.write_restraints_file_0()

        else:
            configuration = {"emstep": self.min_dt, "emtol": self.min_tol}
        minimised_system = self.minimise(system=solvated_system, working_directory=min_dir, configuration=configuration, restraints_file=restraints_file)

        restrained_nvt = self.heat(system=minimised_system,
                                   working_directory=r_nvt_dir,
                                   process_name="r_nvt",
                                   time=self.short_nvt,
                                   start_t=start_temp, end_t=self.temperature,
                                   restraints="all",
                                   timestep=self.short_timestep, 
                                   restraints_file=restraints_file,
                                   configuration=configuration) 
        backbone_restrained_nvt = self.heat(system=restrained_nvt,
                                            process_name="bb_r_nvt",
                                            working_directory=bb_r_nvt_dir,
                                            time=self.nvt,
                                            temperature=self.temperature,
                                            restraints="backbone",
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
                                   restraints="heavy",
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
    

def main():

    parser = argparse.ArgumentParser(description="minimisation and equilibration for meze workflow")

    parser.add_argument("ligand_name",
                        help="ligand name",
                        type=str)
    
    parser.add_argument("protocol_file",
                        help="protocol file containing equilibration options",
                        type=str,
                        default=os.getcwd() + "/afe/protocol.dat")
    
    arguments = parser.parse_args()
    protocol = functions.input_to_dict(file=arguments.protocol_file)

    keys = list(protocol.keys())
    if "metal" in keys:
        metal = True
    else:
        metal = False

    cold_meze = coldMeze(is_metal=metal,
                         group_name=protocol["group name"],
                         ligand_name=arguments.ligand_name,
                         equilibration_directory=protocol["equilibration directory"],
                         input_protein_file=protocol["protein input file"],
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