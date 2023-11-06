"""
Minimise and equilibrate bound and unbound stages.
"""
import Network
from definitions import PICOSECOND, KELVIN, ATM
import functions
import BioSimSpace as bss
import argparse
import os 
import Meze


def run_process(system, protocol, process, working_directory, configuration=None, checkpoint=None):
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
    configuration: dict 
        extra options to pass to the Gromacs process
    checkpoint: str
        path to a checkpoint file forom a previous run; corresponds to the -t flag for gmx grompp

    Return:
    -------
    system: bss.System
        equilibrated or minimised system
    """
    process = bss.Process.Gromacs(system, protocol, name=process, work_dir=working_directory, extra_options=configuration, checkpoint_file=checkpoint)
    # config = process.getConfig()
    # if configuration:
    #     for setting in configuration:
    #         key = setting.split()[0]
    #         try:
    #             index = [i for i, string in enumerate(config) if key in string][0]
    #             config[index] = setting
    #             process.setConfig(config)
    #         except IndexError:
    #             process.addToConfig(setting)
    #             config = process.getConfig()
    process.setArg("-ntmpi", 1)
    process.start()
    process.wait()
    if process.isError():
        print(process.stdout())
        print(process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")
    system = process.getSystem()
    return system


def minimise(system, workdir, min_steps, min_dt, min_tol):
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
    protocol = bss.Protocol.Minimisation(steps=min_steps)
    configuration = {"emstep": min_dt, "emtol": min_tol}
    minimised_system = run_process(system, protocol, "min", workdir, configuration=configuration)
    return minimised_system   


def equilibrate(system, name, workdir, time, start_t=300, end_t=300, temperature=None, pressure=None, configuration=None, restraints=None, checkpoint=None):
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
    protocol = bss.Protocol.Equilibration(runtime=time,
                                          temperature_start=start_t,
                                          temperature_end=end_t,
                                          temperature=temperature,
                                          pressure=pressure,
                                          restraint=restraints)
    equilibrated_system = run_process(system, protocol, name, workdir, configuration, checkpoint)
    return equilibrated_system


def heat_unbound(ligand_name, solvated_network):
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

    directory = functions.mkdir(solvated_network.equilibration_directory + f"/unbound/{ligand_name}/")
    files = functions.read_files(f"{solvated_network.ligand_path}/{ligand_name}_solvated.*")
    solvated_ligand = bss.IO.readMolecules(files)
    print(f"Equilibrating unbound ligand {ligand_name}")
    directories = lambda step: functions.mkdir(directory+step)
    min_directory = directories("min")
    r_nvt_directory = directories("r_nvt")
    nvt_directory = directories("nvt")
    r_npt_directory = directories("r_npt")
    npt_directory = directories("npt")

    minimised_ligand = minimise(system=solvated_ligand, workdir=min_directory, min_steps=solvated_network.min_steps, min_dt=solvated_network.min_dt, min_tol=solvated_network.min_tol)
    start_temp = functions.convert_to_units(0, KELVIN)
    restrained_nvt = equilibrate(system=minimised_ligand,
                                 name="r_nvt",
                                 workdir=r_nvt_directory,
                                 time=solvated_network.short_nvt,
                                 start_t=start_temp, end_t=solvated_network.temperature,
                                 configuration={"dt": 0.0005}, # need to be able to change
                                 restraints="all")
    nvt = equilibrate(system=restrained_nvt,
                      name="nvt",
                      workdir=nvt_directory,
                      time=solvated_network.nvt,
                      temperature=solvated_network.temperature,
                      configuration={"gen-vel": "no"}, #https://github.com/OpenBioSim/biosimspace/issues/168
                      checkpoint=r_nvt_directory + "/r_nvt.cpt")
    restrained_npt = equilibrate(system=nvt,
                                 name="r_npt",
                                 workdir=r_npt_directory,
                                 time=solvated_network.npt,
                                 pressure=solvated_network.pressure,
                                 temperature=solvated_network.temperature,
                                 restraints="heavy",
                                 configuration={"gen-vel": "no"}, #https://github.com/OpenBioSim/biosimspace/issues/168
                                 checkpoint=nvt_directory + "/nvt.cpt")
    equilibrated_molecule = equilibrate(system=restrained_npt,
                                        name="npt",
                                        workdir=npt_directory,
                                        time=solvated_network.npt,
                                        pressure=solvated_network.pressure,
                                        temperature=solvated_network.temperature,
                                        configuration={"gen-vel": "no"}, #https://github.com/OpenBioSim/biosimspace/issues/168
                                        checkpoint=r_npt_directory + "/r_npt.cpt")
    unbound_savename = npt_directory + f"/{ligand_name}"
    bss.IO.saveMolecules(filebase=unbound_savename, system=equilibrated_molecule, fileformat=["PRM7", "RST7"])        

    
def heat_bound(ligand_name, solvated_network):
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
        
    directory = functions.mkdir(solvated_network.equilibration_directory+f"/bound/{ligand_name}/")
    files = functions.read_files(f"{solvated_network.protein_path}/bound_{ligand_name}_solvated.*")
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
    minimised_system = minimise(system=solvated_system, workdir=min_dir, min_steps=solvated_network.min_steps, min_dt=solvated_network.min_dt, min_tol=solvated_network.min_tol)
    restrained_nvt = equilibrate(system=minimised_system,
                                 workdir=r_nvt_dir,
                                 name="r_nvt",
                                 time=solvated_network.short_nvt,
                                 start_t=start_temp, end_t=solvated_network.temperature,
                                 restraints="all",
                                 configuration={"dt": 0.0005}) # need to be able to change
    backbone_restrained_nvt = equilibrate(system=restrained_nvt,
                                          name="bb_r_nvt",
                                          workdir=bb_r_nvt_dir,
                                          time=solvated_network.nvt,
                                          temperature=solvated_network.temperature,
                                          restraints="backbone",
                                          configuration={"gen-vel": "no"}, #https://github.com/OpenBioSim/biosimspace/issues/168
                                          checkpoint=r_nvt_dir + "/r_nvt.cpt")
    nvt = equilibrate(system=backbone_restrained_nvt,
                      name="nvt",
                      workdir=nvt_dir,
                      time=solvated_network.nvt,
                      temperature=solvated_network.temperature,
                      configuration={"gen-vel": "no"}, #https://github.com/OpenBioSim/biosimspace/issues/168
                      checkpoint=bb_r_nvt_dir + "/bb_r_nvt.cpt")
    restrained_npt = equilibrate(system=nvt,
                                 name="r_npt",
                                 workdir=r_npt_dir,
                                 time=solvated_network.npt,
                                 pressure=solvated_network.pressure,
                                 temperature=solvated_network.temperature,
                                 restraints="heavy",
                                 configuration={"gen-vel": "no"}, #https://github.com/OpenBioSim/biosimspace/issues/168
                                 checkpoint=nvt_dir + "/nvt.cpt")
    equilibrated_protein = equilibrate(system=restrained_npt,
                                       name="npt",
                                       workdir=npt_dir,
                                       time=solvated_network.npt,
                                       pressure=solvated_network.pressure,
                                       temperature=solvated_network.temperature,
                                       configuration={"gen-vel": "no"}, #https://github.com/OpenBioSim/biosimspace/issues/168
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
                            )
    elif not metal:
        network = Network.Network(prepared=True,
                                  workdir=protocol["protein directory"],
                                  afe_input_path=protocol["afe input directory"],
                                  equilibration_path=protocol["equilibration directory"],
                                  protein_file=protocol["prepared protein file"],
                                  protein_path=protocol["protein directory"],
                                  ligand_path=protocol["ligand directory"],
                                  group_name=protocol["group name"],
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
        
    heat_unbound(ligand_name=arguments.ligand_name, solvated_network=network)
    heat_bound(ligand_name=arguments.ligand_name, solvated_network=network)
    

if __name__ == "__main__":
    main()