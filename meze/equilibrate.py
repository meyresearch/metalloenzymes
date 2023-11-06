"""
Minimise and equilibrate bound and unbound stages.
"""
import Network
from definitions import FEMTOSECOND, PICOSECOND, KELVIN, ATM
import functions
import BioSimSpace as bss
import argparse
import os 
import Meze


class coldMeze(object):
    def __init__(self, parameters, equilibration_directory, prepared_protein, protein_directory, ligand_directory, min_steps, short_nvt, nvt, npt, min_dt, min_tol, temperature, pressure, short_timestep=0.5, is_metal=True):
        self.parameters = parameters
        self.is_metal = is_metal
        self.equilibration_directory = equilibration_directory
        self.ligand_path = functions.path_exists(ligand_directory)
        self.protein_path = functions.path_exists(protein_directory)
        self.prepared_protein = functions.read_files(prepared_protein + ".*")
        self.short_nvt = functions.convert_to_units(short_nvt, PICOSECOND)
        self.nvt = functions.convert_to_units(nvt, PICOSECOND)
        self.npt = functions.convert_to_units(npt, PICOSECOND)
        self.short_timestep = functions.convert_to_units(short_timestep, FEMTOSECOND)
        self.min_steps = min_steps
        self.min_dt = min_dt
        self.min_tol = min_tol
        self.temperature = functions.convert_to_units(temperature, KELVIN)
        self.pressure = functions.convert_to_units(pressure, ATM)
        

    def run(self, protocol, name, working_directory, configuration=None, checkpoint=None):
        """
        Run a minimisation or equilibration process 
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
        configuration: dict 
            extra options to pass to the process
        checkpoint: str
            path to a checkpoint file from a previous run; corresponds to the -t flag for gmx grompp

        Return:
        -------
        system: bss.System
            equilibrated or minimised system
        """
        if self.is_metal:
            # process = bss.Process.Amber()
            pass
        else:
            process = bss.Process.Gromacs(self.parameters, protocol, name=name, work_dir=working_directory, extra_options=configuration, checkpoint_file=checkpoint)
            process.setArg("-ntmpi", 1)
        process.start()
        process.wait()
        if process.isError():
            print(process.stdout())
            print(process.stderr())
            raise bss._Exceptions.ThirdPartyError("The process exited with an error!")
        system = process.getSystem()
        return system


    def minimise(self, working_directory):
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
        if self.is_metal:
            # configuration = {"nmropt:" 1}
            pass
        else:
            configuration = {"emstep": self.min_dt, "emtol": self.min_tol}
        minimised_system = self.run(self.parameters, protocol, "min", working_directory, configuration=configuration)
        return minimised_system   


    def heat(self, process_name, working_directory, time, timestep=2, start_t=300, end_t=300, temperature=None, pressure=None, configuration=None, restraints=None, checkpoint=None):
        """
        Run NVT or NPT equilibration

        Parameters:
        -----------
        system: bss.System
            system to be equilibrated
        name: str
            name of equilibration process
        workdir: str
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
        # need to add restraint name list
        equilibrated_system = self.run(self.parameters, protocol, process_name, working_directory, configuration, checkpoint) 
        return equilibrated_system


    def heat_unbound(self, ligand_name):
        """
        Perform minimisation and NVT and NPT equilibrations on ligand

        Parameters:
        -----------
        ligand_name: str
            ligand name
        solvated_network: Network
            solvated network object
    
        Return:
        -------

        """

        directory = functions.mkdir(self.equilibration_directory + f"/unbound/{ligand_name}/")
        files = functions.read_files(f"{self.ligand_path}/{ligand_name}_solvated.*")
        solvated_ligand = bss.IO.readMolecules(files)
        print(f"Equilibrating unbound ligand {ligand_name}")
        directories = lambda step: functions.mkdir(directory+step)
        min_directory = directories("min")
        r_nvt_directory = directories("r_nvt")
        nvt_directory = directories("nvt")
        r_npt_directory = directories("r_npt")
        npt_directory = directories("npt")

        minimised_ligand = self.minimise(system=solvated_ligand, workdir=min_directory)
        start_temp = functions.convert_to_units(0, KELVIN)
        restrained_nvt = self.heat(system=minimised_ligand,
                                   name="r_nvt",
                                   workdir=r_nvt_directory,
                                   time=self.short_nvt,
                                   start_t=start_temp, end_t=self.temperature,
                                   timestep=self.short_timestep,
                                   restraints="all")
        nvt = self.heat(system=restrained_nvt,
                        name="nvt",
                        workdir=nvt_directory,
                        time=self.nvt,
                        temperature=self.temperature,
                        checkpoint=r_nvt_directory + "/r_nvt.cpt")
        restrained_npt = self.heat(system=nvt,
                                   name="r_npt",
                                   workdir=r_npt_directory,
                                   time=self.npt,
                                   pressure=self.pressure,
                                   temperature=self.temperature,
                                   restraints="heavy",
                                   checkpoint=nvt_directory + "/nvt.cpt")
        equilibrated_molecule = self.heat(system=restrained_npt,
                                            name="npt",
                                            workdir=npt_directory,
                                            time=self.npt,
                                            pressure=self.pressure,
                                            temperature=self.temperature,
                                            checkpoint=r_npt_directory + "/r_npt.cpt")
        unbound_savename = npt_directory + f"/{ligand_name}"
        bss.IO.saveMolecules(filebase=unbound_savename, system=equilibrated_molecule, fileformat=["PRM7", "RST7"])        

    
    def heat_bound(self, ligand_name):
        """
        Perform minimisation and NVT and NPT equilibrations on bound ligand 

        Parameters:
        -----------
        ligand_name: str
            ligand name
        solvated_network: Network
            solvated network object

        Return:
        -------
        
        """ 
            
        directory = functions.mkdir(self.equilibration_directory+f"/bound/{ligand_name}/")
        files = functions.read_files(f"{self.protein_path}/bound_{ligand_name}_solvated.*")
        solvated_system = bss.IO.readMolecules(files)
        directories = lambda step: functions.mkdir(directory+step)
        min_dir = directories("min")
        r_nvt_dir = directories("r_nvt")
        bb_r_nvt_dir = directories("bb_r_nvt")
        nvt_dir = directories("nvt")
        r_npt_dir = directories("r_npt")
        npt_dir = directories("npt")     
        start_temp = functions.convert_to_units(0, KELVIN)
        print(f"Equilibrating bound ligand {ligand_name}")
        minimised_system = self.minimise(system=solvated_system, workdir=min_dir)
        restrained_nvt = self.heat(system=minimised_system,
                                   workdir=r_nvt_dir,
                                   name="r_nvt",
                                   time=self.short_nvt,
                                   start_t=start_temp, end_t=self.temperature,
                                   restraints="all",
                                   timestep=self.short_timestep) # need to be able to change
        backbone_restrained_nvt = self.heat(system=restrained_nvt,
                                            name="bb_r_nvt",
                                            workdir=bb_r_nvt_dir,
                                            time=self.nvt,
                                            temperature=self.temperature,
                                            restraints="backbone",
                                            checkpoint=r_nvt_dir + "/r_nvt.cpt")
        nvt = self.heat(system=backbone_restrained_nvt,
                        name="nvt",
                        workdir=nvt_dir,
                        time=self.nvt,
                        temperature=self.temperature,
                        checkpoint=bb_r_nvt_dir + "/bb_r_nvt.cpt")
        restrained_npt = self.heat(system=nvt,
                                   name="r_npt",
                                   workdir=r_npt_dir,
                                   time=self.npt,
                                   pressure=self.pressure,
                                   temperature=self.temperature,
                                   restraints="heavy",
                                   checkpoint=nvt_dir + "/nvt.cpt")
        equilibrated_protein = self.heat(system=restrained_npt,
                                         name="npt",
                                         workdir=npt_dir,
                                         time=self.npt,
                                         pressure=self.pressure,
                                         temperature=self.temperature,
                                         checkpoint=r_npt_dir + "/r_npt.cpt")
        bound_savename = npt_dir + f"/bound_{ligand_name}"
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
    
    if metal:
        network = Meze.Meze(prepared=True,
                            metal=protocol["metal"],
                            cut_off=protocol["cutoff"],
                            force_constant_0=protocol["force constant"],
                            workdir=protocol["protein directory"],
                            equilibration_path=protocol["equilibration directory"],
                            outputs=protocol["outputs"],
                            protein_file=protocol["prepared protein file"],
                            protein_path=protocol["protein directory"],
                            ligand_path=protocol["ligand directory"],
                            group_name=protocol["group name"],
                            engine=protocol["engine"],
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
        network = coldMeze(
                                  workdir=protocol["protein directory"],

                                  protein_file=protocol["prepared protein file"],
                                  protein_path=protocol["protein directory"],
                                  ligand_path=protocol["ligand directory"],
                                  group_name=protocol["group name"],

                                  min_steps=protocol["minimisation steps"],
                                  short_nvt=protocol["short nvt"],
                                  nvt=protocol["nvt"],
                                  npt=protocol["npt"],
                                  min_dt=protocol["minimisation stepsize"],
                                  min_tol=protocol["minimisation tolerance"],
                                  repeats=protocol["repeats"],
                                  temperature=protocol["temperature"],
                                  pressure=protocol["pressure"])
        
        heat_unbound(ligand_name=arguments.ligand_name, solvated_network=network)
        heat_bound(ligand_name=arguments.ligand_name, solvated_network=network)
    

if __name__ == "__main__":
    main()