"""
Minimise and equilibrate bound and unbound stages.
"""
from definitions import PICOSECOND, KELVIN, ATM
import functions
import BioSimSpace as bss


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
    process = bss.Process.Gromacs(system, protocol, name=process, work_dir=working_directory)
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
    process.setArg("-ntmpi", 1)
    process.start()
    process.wait()
    if process.isError():
        print(process.stdout())
        print(process.stderr())
        raise bss._Exceptions.ThirdPartyError("The process exited with an error!")
    system = process.getSystem()
    return system


def minimise(system, workdir, AFE):
    """
    Minimise the system using Gromacs

    Parameters:
    -----------
    system: bss.System
        system to be minimised
    workdir: str
        current working dir
    AFE: AlchemicalFreeEnergy
        AFE class object

    Return:
    -------
    minimsed_system: bss.System
        minimised system
    """
    print("Minimisation")
    protocol = bss.Protocol.Minimisation(steps=AFE.min_steps)
    minimised_system = run_process(system, protocol, "min", workdir, configuration=[f"emstep = {AFE.emstep}", f"emtol = {AFE.emtol}"])
    return minimised_system


def equilibrate(system, workdir, name, runtime, config=None, start_t=300, end_t=300, temperature=None, pressure=None, restraints=None):
    """
    Run NVT or NPT equilibration

    Parameters:
    -----------
    system: bss.System
        system to be equilibrated
    workdir: str
        working directory for equilibration step
    name: str
        process name
    runtime: str/float
        runtime for equilibration

    Return:
    -------
    equilibrated_system: bss.System
        equilibrated system
    """
    print(f"{name.upper()}")
    start_t = functions.convert_to_units(start_t, KELVIN)
    end_t = functions.convert_to_units(end_t, KELVIN)
    if temperature:
        temperature = functions.convert_to_units(temperature, KELVIN)
    if pressure:
        pressure = functions.convert_to_units(pressure, ATM)

    protocol = bss.Protocol.Equilibration(runtime=runtime,
                                          temperature_start=start_t,
                                          temperature_end=end_t,
                                          temperature=temperature,
                                          pressure=pressure,
                                          restraint=restraints)
    equilibrated_system = run_process(system, protocol, name, workdir, config)
    return equilibrated_system

def heat_meze(idx, Protein, Network, AFE):
    """
    Run minimisation and equilibration in unbound and bound stages

    Parameters:
    -----------
    idx: int
        ligand index
    Protein: Protein
        Protein class object
    Network: Network
        Network class object
    AFE: AlchemicalFreeEnergy
        AFE Class object

    Return:
    -------
    """
    ligand_number = Network.names[idx].split("_")[-1]
    unbound_directory = AFE.create_directory(AFE.equilibration_dir + f"/unbound/ligand_{ligand_number}/")
    solvated_ligand = Network.ligands[idx]

    unbound_directories = lambda step: functions.mkdir(unbound_directory+step)
    unbound_min = unbound_directories("min")
    unbound_r_nvt = unbound_directories("r_nvt")
    unbound_nvt = unbound_directories("nvt")
    unbound_r_npt = unbound_directories("r_npt")
    unbound_npt = unbound_directories("npt")    

    minimised_ligand = minimise(system=solvated_ligand, 
                                workdir=unbound_min, 
                                AFE=AFE)
    unbound_restrained_nvt = equilibrate(system=minimised_ligand, 
                                         workdir=unbound_r_nvt, 
                                         name="r_nvt",
                                         runtime=AFE.short_nvt,
                                         config=["dt = 0.0005"],
                                         start_t=0, end_t=300, 
                                         restraints="all")
    unbound_nvt = equilibrate(system=unbound_restrained_nvt,
                              workdir=unbound_nvt,
                              name="nvt",
                              runtime=AFE.nvt,
                              temperature=300)
    unbound_restrained_npt = equilibrate(system=unbound_nvt,
                                         workdir=unbound_r_npt,
                                         name="r_npt",
                                         runtime=AFE.npt,
                                         pressure=1, temperature=300,
                                         restraints="heavy")
    equilibrated_ligand = equilibrate(system=unbound_restrained_npt,
                                      workdir=unbound_npt,
                                      name="npt",
                                      runtime=AFE.npt,
                                      pressure=1, temperature=300)
    unbound_savename = unbound_npt + f"/ligand_{ligand_number}"
    bss.IO.saveMolecules(filebase=unbound_savename, system=equilibrated_ligand, fileformat=["PRM7", "RST7"]) 
    Network.ligands[idx] = equilibrated_ligand

    bound_directory = AFE.create_directory(AFE.equilibration_dir + f"/bound/ligand_{ligand_number}/")
    solvated_system = Protein.complex
    bound_directories = lambda step: functions.mkdir(bound_directory+step)
    bound_min = bound_directories("min")
    bound_r_nvt = bound_directories("r_nvt")
    bound_bb_r_nvt = bound_directories("bb_r_nvt")
    bound_nvt = bound_directories("nvt")
    bound_r_npt = bound_directories("r_npt")
    bound_npt = bound_directories("npt") 

    minimised_system = minimise(system=solvated_system, workdir=bound_min, AFE=AFE)
    bound_restrained_nvt = equilibrate(system=minimised_system,
                                       workdir=bound_r_nvt,
                                       name="r_nvt",
                                       runtime=AFE.short_nvt,
                                       config=["dt = 0.0005"],
                                       start_t=0, end_t=300,
                                       restraints="all")
    bound_backbone_restrained_nvt = equilibrate(system=bound_restrained_nvt,
                                                workdir=bound_bb_r_nvt,
                                                name="bb_r_nvt",
                                                runtime=AFE.nvt,
                                                temperature=300,
                                                restraints="backbone")
    bound_nvt_equilibration = equilibrate(system=bound_backbone_restrained_nvt,
                                          workdir=bound_nvt,
                                          name="nvt",
                                          runtime=AFE.nvt,
                                          temperature=300)
    bound_restrained_npt = equilibrate(system=bound_nvt_equilibration,
                                       workdir=bound_r_npt,
                                       name="r_npt",
                                       runtime=AFE.npt,
                                       pressure=1, temperature=300,
                                       restraints="heavy")
    equilibrated_system = equilibrate(system=bound_restrained_npt,
                                      workdir=bound_npt,
                                      name="npt",
                                      runtime=AFE.npt,
                                      pressure=1, temperature=300)
    bound_savename = bound_npt + f"/system_{ligand_number}"
    bss.IO.saveMolecules(filebase=bound_savename, system=equilibrated_system, fileformat=["PRM7", "RST7"])
    Network.bound[idx] = equilibrated_system
    return Protein, Network, AFE

def main():
    pass


if __name__ == "__main__":
    main()