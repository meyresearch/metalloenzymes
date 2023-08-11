"""
Minimise and equilibrate bound and unbound stages.
"""
from definitions import PICOSECOND, KELVIN, ATM
import functions
import BioSimSpace as bss


def run_process(system, protocol, process, working_directory, config=None):
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
    if config:
        for setting in config:
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
    minimised_system = run_process(system, protocol, "min", workdir, config=[f"emstep = {AFE.emstep}", f"emtol = {AFE.emtol}"])
    return minimised_system

# def equilibrate(AFE)
def equilibrate(runtime, t_start=300, t_end=300, temperature=None, pressure=None, restraints=None):
    # this was dumb
    # i think i just rewrote their code
    # lol
    # make these attributes of AFE: runtime, t_start=300, t_end=300, temperature=None, pressure=None, restraints=None

    # DO THESE CONVERSIONS IN THE CLASS CONSTRUCTOR
    t_start = functions.convert_to_units(t_start, KELVIN)
    t_end = functions.convert_to_units(t_start, KELVIN)
    time = functions.convert_to_units(runtime, PICOSECOND)

    if temperature:
        temperature = functions.convert_to_units(temperature, KELVIN)
    if pressure:
        pressure = functions.convert_to_units(pressure, ATM)

    protocol = bss.Protocol.Equilibration(runtime=time,
                                          temperature_start=t_start,
                                          temperature_end=t_end,
                                          temperature=temperature,
                                          pressure=pressure,
                                          restraint=restraints)
    equilibrated_system = run_process()


def heat_meze(idx, Protein, Network, AFE):

    # min_steps = AFE.min_steps
    # short_nvt = AFE.short_nvt
    # nvt, npt = AFE.nvt, AFE.npt
    ligand_number = Network.names[idx].split("_")[-1]
    unbound_directory = AFE.create_directory(AFE.equilibration_dir + f"/unbound/ligand_{ligand_number}/")
    solvated_ligand = Network.ligands[idx]

    unbound_directories = lambda step: functions.mkdir(unbound_directory+step)
    unbound_min = unbound_directories("min")
    unbound_r_nvt = unbound_directories("r_nvt")
    unbound_nvt = unbound_directories("nvt")
    unbound_r_npt = unbound_directories("r_nvt")
    unbound_npt = unbound_directories("nvt")    

    minimised_ligand = minimise(solvated_ligand, unbound_min, AFE)




    bound_directory = AFE.create_directory(AFE.equilibration_dir + f"/bound/ligand_{ligand_number}/")
    solvated_system = Protein.complex
    bound_directories = lambda step: functions.mkdir(bound_directory+step)
    bound_min = bound_directories("min")
    bound_r_nvt = bound_directories("r_nvt")
    bound_bb_r_nvt = bound_directories("bb_r_nvt")
    bound_nvt = bound_directories("nvt")
    bound_r_npt = bound_directories("r_nvt")
    bound_npt = bound_directories("nvt") 


def main():
    pass


if __name__ == "__main__":
    main()