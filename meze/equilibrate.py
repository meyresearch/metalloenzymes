"""
Minimise and equilibrate bound and unbound stages.
"""
from definitions import PICOSECOND, KELVIN, ATM
import functions
import BioSimSpace as bss
import argparse
import os 


def write_slurm_script(path, log_dir, project_dir, equil_dir, min_steps, min_dt, min_tol, short_nvt, nvt, npt):
    """
    Write a slurm script for running equilibration for ligands

    Parameters:
    -----------
    path: str
        path where slurm script will be saved
    log_dir: 
        directory for outputting slurm logs
    project_dir: str
        project working directory
    equil_dir: 
        full path to /equilibration/
    min_steps: int
        number of minimisation steps
    min_dt: float 
        step size for minimisation
    min_tol: float
        force tolerance for minimisation
    short_nvt: float
        runtime in ps for short NVT equilibration
    nvt: float
        runtime in ps for NVT equilibration
    npt: float
        runtime in ps for NPT equilibration

    Return:
    -------
    file: str
        slurm script 
    """
    file = path + "slurm_heat_meze.sh"
    with open(file, "w") as f:
        f.write("#!/bin/python\n")
        f.write("\n")
        f.write(f"#SBATCH -o {log_dir}/heat_%a.slurm.out\n")
        f.write(f"#SBATCH -e {log_dir}/heat_%a.slurm.err\n")
        f.write("#SBATCH -n 1\n")
        f.write("#SBATCH --gres=gpu:1\n")
        f.write("#SBATCH --cpus-per-gpu=10\n")
        f.write("#SBATCH --mem 4096\n")
        f.write("#SBATCH --job-name=heat_meze\n")
        f.write(f"export \"MEZEHOME\"={os.path.realpath(__file__)}\n") # installation?
        f.write(f"min_steps={min_steps}\n")
        f.write(f"min_dt={min_dt}\n")
        f.write(f"min_tol={min_tol}\n")
        f.write(f"short_nvt={short_nvt}\n")
        f.write(f"nvt={nvt}\n")
        f.write(f"npt={npt}\n")
        f.write(f"project_dir={project_dir}\n")
        f.write(f"equilibration_dir={equil_dir}\n")
        f.write("LIG_NUMBER=$SLURM_ARRAY_TASK_ID\n")
        f.write(f"python $MEZEHOME/equilibrate.py $LIG_NUMBER \
                                                  $equilibration_dir\n \
                                                  $project_dir\n \
                                                  $min_steps\n \
                                                  $min_dt\n \
                                                  $min_tol\n \
                                                  $short_nvt\n \
                                                  $nvt\n \
                                                  $npt\n")
    return file


def slurm_heat(n_ligands, script):
    """
    Run equilibration with slurm

    Parameters:
    -----------
    n_ligands: int
        number of ligands to equilibrate
    script: str
        full path to equilibration slurm script 

    Return:
    -------
    """
    os.system(f"sbatch --array=0-{n_ligands} {script}")
    

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
    process = bss.Process.Gromacs(system, protocol, name=process, work_dir=working_directory, #exe="/usr/local/gromacs/bin/gmx_mpi"
                                  )
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
    configuration = [f"emstep = {min_dt}", f"emtol = {min_tol}"]
    minimised_system = run_process(system, protocol, "min", workdir, configuration=configuration)
    return minimised_system   


def equilibrate(system, name, workdir, time, start_t=300, end_t=300, temperature=None, pressure=None, configuration=None, restraints=None):
    """
    Run NVT or NPT equilibration

    Parameters:
    -----------
    system: bss.System
        system to be equilibrated

    Return:
    -------
    equilibrated_system: bss.System
        equilibrated system
    """
    protocol = bss.Protocol.Equilibration(runtime=time,
                                            temperature_start=start_t,
                                            temperature_end=end_t,
                                            temperature=temperature,
                                            pressure=pressure,
                                            restraint=restraints)
    equilibrated_system = run_process(system, protocol, name, workdir, configuration)
    return equilibrated_system


def heat_unbound(ligand_number, equilibration_dir, project_dir, short_nvt, nvt, npt, temperature=300., pressure=1.):
    """
    Perform minimisation and NVT and NPT equilibrations on ligand

    Parameters:
    -----------
    ligand_number: int
        ligand number used in ligand file name
    equilibration_dir: str
        full path to /equilibration/
    project_dir: str
        full path to working directory for project
    short_nvt: float
        runtime in ps for short NVT equilibration
    nvt: float
        runtime in ps for NVT equilibration
    npt: float 
        runtime in ps for NPT equilibration
    temperature: bss.Units.Temperature
        system temperature in Kelvin
    pressure: bss.Unites.atm
        system pressure in ATM

    Return:
    -------

    """
    temperature = functions.convert_to_units(temperature, KELVIN)
    pressure = functions.convert_to_units(pressure, ATM)
    directory = functions.mkdir(equilibration_dir + f"/unbound/ligand_{ligand_number}/")
    files = functions.read_files(f"{project_dir}/inputs/ligands/ligand_{ligand_number}_solvated.*")
    solvated_ligand = bss.IO.readMolecules(files)
    print(f"Equilibrating unbound ligand {ligand_number}")
    directories = lambda step: functions.mkdir(directory+step)
    min_directory = directories("min")
    r_nvt_directory = directories("r_nvt")
    nvt_directory = directories("nvt")
    r_npt_directory = directories("r_npt")
    npt_directory = directories("npt")

    minimised_ligand = minimise(system=solvated_ligand, workdir=min_directory)
    start_temp = functions.convert_to_units(0, KELVIN)
    restrained_nvt = equilibrate(system=minimised_ligand,
                                        name="r_nvt",
                                        workdir=r_nvt_directory,
                                        time=short_nvt,
                                        start_t=start_temp, end_t=temperature,
                                        configuration=["dt = 0.0005"], # need to be able to change
                                        restraints="all")
    nvt = equilibrate(system=restrained_nvt,
                             name="nvt",
                             workdir=nvt_directory,
                             time=nvt,
                             temperature=temperature)
    restrained_npt = equilibrate(system=nvt,
                                        name="r_npt",
                                        workdir=r_npt_directory,
                                        time=npt,
                                        pressure=pressure,
                                        temperature=temperature,
                                        restraints="heavy")
    equilibrated_molecule = equilibrate(system=restrained_npt,
                                               name="npt",
                                               workdir=npt_directory,
                                               time=npt,
                                               pressure=pressure,
                                               temperature=temperature)
    unbound_savename = npt_directory + f"/ligand_{ligand_number}"
    bss.IO.saveMolecules(filebase=unbound_savename, system=equilibrated_molecule, fileformat=["PRM7", "RST7"])        

    
def heat_bound(ligand_number, equilibration_dir, project_dir, short_nvt, nvt, npt, temperature=300., pressure=1.):
    """
    Perform minimisation and NVT and NPT equilibrations on bound ligand 

    Parameters:
    -----------
    ligand_number: int
        ligand number used in ligand file name
    equilibration_dir: str
        full path to /equilibration/
    project_dir: str
        full path to working directory for project
    short_nvt: float
        runtime in ps for short NVT equilibration
    nvt: float
        runtime in ps for NVT equilibration
    npt: float 
        runtime in ps for NPT equilibration
    temperature: bss.Units.Temperature
        system temperature in Kelvin
    pressure: bss.Unites.atm
        system pressure in ATM

    Return:
    -------
    
    """ 
    temperature = functions.convert_to_units(temperature, KELVIN)
    pressure = functions.convert_to_units(pressure, ATM)           
    directory = functions.mkdir(equilibration_dir+f"/bound/ligand_{ligand_number}/")
    files = functions.read_files(f"{project_dir}/inputs/ligands/system_{ligand_number}_solvated.*")
    solvated_system = bss.IO.readMolecules(files)
    directories = lambda step: functions.mkdir(directory+step)
    min_dir = directories("min")
    r_nvt_dir = directories("r_nvt")
    bb_r_nvt_dir = directories("bb_r_nvt")
    nvt_dir = directories("nvt")
    r_npt_dir = directories("r_npt")
    npt_dir = directories("npt")     
    start_temp = functions.convert_to_units(0, KELVIN)
    print(f"Equilibrating bound ligand {ligand_number}")
    minimised_system = minimise(system=solvated_system, workdir=min_dir)
    restrained_nvt = equilibrate(system=minimised_system,
                                        workdir=r_nvt_dir,
                                        name="r_nvt",
                                        time=short_nvt,
                                        start_t=start_temp, end_t=temperature,
                                        restraints="all",
                                        configuration=["dt = 0.0005"]) # need to be able to change
    backbone_restrained_nvt = equilibrate(system=restrained_nvt,
                                                 name="bb_r_nvt",
                                                 workdir=bb_r_nvt_dir,
                                                 time=nvt,
                                                 temperature=temperature,
                                                 restraints="backbone")
    nvt = equilibrate(system=backbone_restrained_nvt,
                             name="nvt",
                             workdir=nvt_dir,
                             time=nvt,
                             temperature=temperature)
    restrained_npt = equilibrate(system=nvt,
                                        name="r_npt",
                                        workdir=r_npt_dir,
                                        time=npt,
                                        pressure=pressure,
                                        temperature=temperature,
                                        restraints="heavy")
    equilibrated_protein = equilibrate(system=restrained_npt,
                                              name="npt",
                                              workdir=npt_dir,
                                              time=npt,
                                              pressure=pressure,
                                              temperature=temperature)
    bound_savename = npt_dir + f"/system_{ligand_number}"
    bss.IO.saveMolecules(filebase=bound_savename, system=equilibrated_protein, fileformat=["PRM7", "RST7"])     
   

def main():

    parser = argparse.ArgumentParser(description="minimisation and equilibration for meze workflow")

    parser.add_argument("ligand-number",
                        help="ligand number used in ligand file name",
                        type=str)
    
    parser.add_argument("equil-dir",
                        dest="equil_dir",
                        help="full path to /equilibration/",
                        default=functions.path_exists(os.getcwd() + "/equilibration/"))
    
    parser.add_argument("project-dir",
                        dest="project_dir",
                        help="full path to project working directory",
                        default=os.getcwd())   

    parser.add_argument("minimisation-steps",
                        dest="min_steps",
                        help="number of minimisation steps for equilibration stage",
                        type=int,
                        default=500)

    parser.add_argument("short-nvt-runtime",
                        dest="short_nvt",
                        help="runtime in ps for short NVT equilibration",
                        type=float,
                        default=5) 

    parser.add_argument("nvt-runtime",
                        dest="nvt",
                        help="runtime in ps for NVT equilibration",
                        type=float,
                        default=50)

    parser.add_argument("npt-runtime",
                        dest="npt",
                        help="runtime in ps for NPT equilibration",
                        type=float,
                        default=200)

    parser.add_argument("em-step",
                        dest="emstep",
                        help="Step size for energy minimisation",
                        type=float,
                        default=0.01)

    parser.add_argument("em-tolerance",
                        dest="emtol",
                        help="kJ mol-1 nm-1, Maximum force tolerance for energy minimisation",
                        type=float,
                        default=1000)
    arguments = parser.parse_args()

    heat_unbound(ligand_number=arguments.ligand_number,
                 equilibration_dir=arguments.equil_dir,
                 project_dir=arguments.project_dir,
                 short_nvt=functions.convert_to_units(arguments.short_nvt, PICOSECOND),
                 nvt=functions.convert_to_units(arguments.nvt, PICOSECOND),
                 npt=functions.convert_to_units(arguments.npt, PICOSECOND))
    heat_bound(ligand_number=arguments.ligand_number,
               equilibration_dir=arguments.equil_dir,
               project_dir=arguments.project_dir,
               short_nvt=functions.convert_to_units(arguments.short_nvt, PICOSECOND),
               nvt=functions.convert_to_units(arguments.nvt, PICOSECOND),
               npt=functions.convert_to_units(arguments.npt, PICOSECOND))

if __name__ == "__main__":
    main()