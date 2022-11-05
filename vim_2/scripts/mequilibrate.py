import warnings
warnings.filterwarnings("ignore")
import BioSimSpace as bss
import os
import nglview as nv
import subprocess as sp
import argparse
import sys


def change_barostat(system: bss._SireWrappers._system.System,
                    protocol: bss.Protocol._equilibration.Equilibration,
                    work_directory: str,
                    process_name: str) -> None:
    """
    Change barostat in .mdp file for NPT runs
    @param: system
    @param: protocol
    @param: work_directory
    @param: process name
    @return: None
    """
    try:
        process = bss.Process.Gromacs(system,
                                      protocol,
                                      name=process_name,
                                      work_dir=work_directory)
    except RuntimeError:

        with open(f"{work_directory}/{process_name}.mdp", "r") as mdp:
            lines = mdp.read()

        new_lines = lines.replace("pcoupl = berendsen", "pcoupl = C-rescale")
        with open(f"{work_directory}/{process_name}.mdp", "w") as mdp:
            mdp.write(new_lines)


def write_mdp(work_dir: str, process_name: str, custom_options=[], default=True):
    """
    Write an mdp file for minimisation
    """
    default_options = ["nstlog = 100",
                       "integrator = steep",
                       "nstcgsteep = 1000",
                       "nsteps = 2000",
                       "pbc = xyz",
                       "cutoff-scheme = Verlet",
                       "ns-type = grid",
                       "nstlist = 20",
                       "coulombtype = PME",
                       "rvdw = 1.0",
                       "rcoulomb = 1.0",
                       "DispCorr = EnerPres",
                       "vdwtype = Cut-off",
                       "refcoord-scaling = all"]
    if default: 
        # print("Writing default mdp file.")
        with open(f"{work_dir}/{process_name}.mdp", "w") as mdp:
            for opt in default_options:
                mdp.write(f"{opt}\n")
        if len(custom_options) > 0:
            # print("Appending custom options.")
            with open(f"{work_dir}/{process_name}.mdp", "a") as mdp:
                for opt in custom_options:
                    mdp.write(f"{opt}\n")
    elif not default and custom_options != []:
        # print("Writing custom mdp file.")
        with open(f"{work_dir}/{process_name}.mdp", "w") as mdp:
            for opt in custom_options:
                mdp.write(f"{opt}\n")
    else:
        print("Custom options not specified. Not writing mdp file.")


def write_script(work_dir: str, process_name: str, restrained=False, previous="") -> str:
    script_name = work_dir + "/" + process_name + ".sh"
    with open(script_name, "w") as script:
        script.write("#!/bin/bash \n")
        script.write(f"cd {work_dir} \n")
        if not restrained and previous == "":
            script.write(f"gmx grompp -f {process_name}.mdp -c {process_name}.gro -p {process_name}.top -o {process_name}.tpr \n")
        elif restrained and previous == "min": 
            script.write(f"gmx grompp -f {process_name}.mdp -c ../{previous}/{previous}.gro -r ../{previous}/{previous}.gro -p {process_name}.top -o {process_name}.tpr \n")
        elif not restrained and previous != "":
            script.write(f"gmx grompp -f {process_name}.mdp -c ../{previous}/{previous}.gro -p {process_name}.top -t ../{previous}/{previous}.cpt -o {process_name}.tpr \n")
        elif restrained and previous != "":
            script.write(f"gmx grompp -f {process_name}.mdp -c ../{previous}/{previous}.gro -r ../{previous}/{previous}.gro -p {process_name}.top -t ../{previous}/{previous}.cpt -o {process_name}.tpr \n") 
        script.write(f"nohup gmx mdrun -v -deffnm {process_name} -nt 1 -nb gpu > nohup_{process_name}.out &")
    os.system(f"chmod +x {script_name}")
    executable = script_name.split("/")[-1]
    return executable


def write_log_file(work_dir: str, process_name: str, subprocess: sp.CompletedProcess) -> None:
    """
    Write custom log file for gromacs processes.
    """
    try:
        if "error" in subprocess.stderr or "error" in subprocess.stdout:
            print(f"There was an error with {process_name}.sh. Please check the log file.")
        with open(f"{work_dir}/{process_name}.out", "w") as log_file:
            log_file.writelines(subprocess.stdout)
            log_file.writelines(subprocess.stderr)
    except FileNotFoundError:
        print(f"File {work_dir}/{process_name}.out does not exist.")
    except NotADirectoryError as e:
        print(e)


def create_dirs(directory: str) -> str:
    """
    Check if directory exists, and if not, create the directory.
    """ 
    path_exists = os.path.exists(directory)
    if not path_exists:
        os.makedirs(directory)
    return directory


parser = argparse.ArgumentParser(description="perform minimisation and equilibration for AFE ligands")
parser.add_argument("ligands",
                    type=str,
                    help="file containing ligand names")

arguments = parser.parse_args()
input_file = arguments.ligands


if not os.path.isfile(input_file):
    print(f"The input file {input_file} does not exist")
    sys.exit()

minimisation_steps = 2000
picosecond = bss.Units.Time.picosecond
kelvin = bss.Units.Temperature.kelvin
atm = bss.Units.Pressure.atm
short_nvt_runtime = 5 * picosecond
nvt_runtime = 50 * picosecond
npt_runtime = 200 * picosecond

ligand_datafile = open("ligands.dat", "r")
ligand_lines = ligand_datafile.readlines()
ligand_datafile.close()
n_ligands = len(ligand_lines)  

print("Processing unbound stage.")
for i in range(n_ligands):
    ligand_name = ligand_lines[i].rstrip()
    ligand_work_dir = create_dirs(f"../runs/equilibration/unbound/{ligand_name}")
    ligand_min_unbound_dir = create_dirs(f"{ligand_work_dir}/min" )
    ligand_r_nvt_unbound_dir = create_dirs(f"{ligand_work_dir}/r_nvt")
    ligand_nvt_unbound_dir = create_dirs(f"{ligand_work_dir}/nvt")  
    ligand_r_npt_unbound_dir = create_dirs(f"{ligand_work_dir}/r_npt")
    ligand_npt_unbound_dir = create_dirs(f"{ligand_work_dir}/npt")
    
    # os.system(f"mkdir {ligand_work_dir}/min/")    
    # os.system(f"mkdir {ligand_work_dir}/r_nvt/")
    # os.system(f"mkdir {ligand_work_dir}/nvt/")  
    # os.system(f"mkdir {ligand_work_dir}/r_npt/")      
    # os.system(f"mkdir {ligand_work_dir}/npt/")      

    ligand_min_script = write_script(ligand_min_unbound_dir, "min")
    ligand_r_nvt_script = write_script(ligand_r_nvt_unbound_dir, "r_nvt", restrained=True, previous="min")
    ligand_nvt_script = write_script(ligand_nvt_unbound_dir, "nvt", previous="r_nvt")
    ligand_r_npt_script = write_script(ligand_r_npt_unbound_dir, "r_npt", restrained=True, previous="nvt")
    ligand_npt_script = write_script(ligand_npt_unbound_dir, "npt", previous="r_npt")

    solvated_ligand = bss.IO.readMolecules([f"../inputs/ligands/{ligand_name}_solv.prm7",
                                            f"../inputs/ligands/{ligand_name}_solv.rst7"])
    ligand_minimisation_protocol = bss.Protocol.Minimisation(steps=minimisation_steps)
    ligand_minimisation_process = bss.Process.Gromacs(solvated_ligand,
                                                      ligand_minimisation_protocol,
                                                      name="min",
                                                      work_dir=ligand_min_unbound_dir)
    write_mdp(ligand_min_unbound_dir, "min")
    min_sp = sp.run(["sh", f"{ligand_min_unbound_dir}/{ligand_min_script}"], capture_output=True, text=True)
    write_log_file(ligand_min_unbound_dir, "min", min_sp)

    minimised_ligand = bss.IO.readMolecules([f"{ligand_min_unbound_dir}/min.gro",
                                             f"{ligand_min_unbound_dir}/min.top"])
    ligand_r_nvt_protocol = bss.Protocol.Equilibration(runtime=short_nvt_runtime,
                                                       temperature_start=0*kelvin,
                                                       temperature_end=300*kelvin,
                                                       restraint="all")
    ligand_r_nvt_process = bss.Process.Gromacs(minimised_ligand,
                                               ligand_r_nvt_protocol,
                                               name="r_nvt",
                                               work_dir=ligand_r_nvt_unbound_dir)
    r_nvt_sp = sp.run(["sh", f"{ligand_r_nvt_unbound_dir}/{ligand_r_nvt_script}"], capture_output=True, text=True)
    write_log_file(ligand_r_nvt_unbound_dir, "r_nvt", r_nvt_sp)  

    restrained_nvt_ligand = bss.IO.readMolecules([f"{ligand_r_nvt_unbound_dir}/r_nvt.gro",
                                                  f"{ligand_r_nvt_unbound_dir}/r_nvt.top"])
    ligand_nvt_protocol = bss.Protocol.Equilibration(runtime=nvt_runtime,
                                                     temperature=300*kelvin)
    ligand_nvt_process = bss.Process.Gromacs(restrained_nvt_ligand,
                                             ligand_nvt_protocol,
                                             name="nvt",
                                             work_dir=ligand_nvt_unbound_dir) 
    nvt_sp = sp.run(["sh", f"{ligand_nvt_unbound_dir}/{ligand_nvt_script}"], capture_output=True, text=True)
    write_log_file(ligand_nvt_unbound_dir, "nvt", nvt_sp)        

    ligand_nvt = bss.IO.readMolecules([f"{ligand_nvt_unbound_dir}/nvt.gro",
                                       f"{ligand_nvt_unbound_dir}/nvt.top"])
    ligand_r_npt_protocol = bss.Protocol.Equilibration(runtime=npt_runtime,
                                                       pressure=1*atm,
                                                       temperature=300*kelvin,
                                                       restraint="heavy")
    change_barostat(ligand_nvt, ligand_r_npt_protocol, ligand_r_npt_unbound_dir, "r_npt")
    r_npt_sp = sp.run(["sh", f"{ligand_r_npt_unbound_dir}/{ligand_r_npt_script}"], capture_output=True, text=True)
    write_log_file(ligand_r_npt_unbound_dir, "r_npt", r_npt_sp)     

    restrained_npt_ligand = bss.IO.readMolecules([f"{ligand_r_npt_unbound_dir}/r_npt.gro",
                                                  f"{ligand_r_npt_unbound_dir}/r_npt.top"])
    ligand_npt_protocol = bss.Protocol.Equilibration(runtime=npt_runtime,
                                                     pressure=1*atm,
                                                     temperature=300*kelvin)    
    change_barostat(restrained_npt_ligand, ligand_npt_protocol, ligand_npt_unbound_dir, "npt")
    npt_sp = sp.run(["sh", f"{ligand_npt_unbound_dir}/{ligand_npt_script}"], capture_output=True, text=True)
    write_log_file(ligand_npt_unbound_dir, "npt", npt_sp)     

print("Processing bound stage.")
for i in range(n_ligands):
    ligand_name = "system_" + ligand_lines[i].rstrip().split("_")[-1]
    ligand_work_dir = create_dirs(f"../runs/equilibration/bound/{ligand_name}")
    system_min_bound_dir = create_dirs(f"{ligand_work_dir}/min") 
    system_r_nvt_bound_dir = create_dirs(f"{ligand_work_dir}/r_nvt")
    system_bb_r_nvt_bound_dir = create_dirs(f"{ligand_work_dir}/bb_r_nvt")  
    system_r_npt_bound_dir = create_dirs(f"{ligand_work_dir}/r_npt")
    system_npt_bound_dir = create_dirs(f"{ligand_work_dir}/npt")
    system_nvt_bound_dir = create_dirs(f"{ligand_work_dir}/nvt")  


    system_min_script = write_script(system_min_bound_dir, "min")
    system_r_nvt_script = write_script(system_r_nvt_bound_dir, "r_nvt", restrained=True, previous="min")
    system_bb_r_nvt_script = write_script(system_bb_r_nvt_bound_dir, "bb_r_nvt", restrained=True, previous="min")
    system_nvt_script = write_script(system_nvt_bound_dir, "nvt", previous="bb_r_nvt")
    system_r_npt_script = write_script(system_r_npt_bound_dir, "r_npt", restrained=True, previous="nvt")
    system_npt_script = write_script(system_npt_bound_dir, "npt", previous="r_npt")

    solvated_system = bss.IO.readMolecules([f"../inputs/protein/{ligand_name}_solv.prm7",
                                            f"../inputs/protein/{ligand_name}_solv.rst7"])
    system_minimisation_protocol = bss.Protocol.Minimisation(steps=minimisation_steps)
    system_minimisation_process = bss.Process.Gromacs(solvated_system,
                                                      system_minimisation_protocol,
                                                      name="min",
                                                      work_dir=system_min_bound_dir)
    write_mdp(system_min_bound_dir, "min")
    min_sp = sp.run(["sh", f"{system_min_bound_dir}/{system_min_script}"], capture_output=True, text=True)
    write_log_file(system_min_bound_dir, "min", min_sp)

    # minimised_ligand = bss.IO.readMolecules([f"{ligand_min_unbound_dir}/min.gro",
    #                                          f"{ligand_min_unbound_dir}/min.top"])
    # ligand_r_nvt_protocol = bss.Protocol.Equilibration(runtime=short_nvt_runtime,
    #                                                    temperature_start=0*kelvin,
    #                                                    temperature_end=300*kelvin,
    #                                                    restraint="all")
    # ligand_r_nvt_process = bss.Process.Gromacs(minimised_ligand,
    #                                            ligand_r_nvt_protocol,
    #                                            name="r_nvt",
    #                                            work_dir=ligand_r_nvt_unbound_dir)
    # r_nvt_sp = sp.run(["sh", f"{ligand_r_nvt_unbound_dir}/{ligand_r_nvt_script}"], capture_output=True, text=True)
    # write_log_file(ligand_r_nvt_unbound_dir, "r_nvt", r_nvt_sp)  

    # restrained_nvt_ligand = bss.IO.readMolecules([f"{ligand_r_nvt_unbound_dir}/r_nvt.gro",
    #                                               f"{ligand_r_nvt_unbound_dir}/r_nvt.top"])
    # ligand_nvt_protocol = bss.Protocol.Equilibration(runtime=nvt_runtime,
    #                                                  temperature=300*kelvin)
    # ligand_nvt_process = bss.Process.Gromacs(restrained_nvt_ligand,
    #                                          ligand_nvt_protocol,
    #                                          name="nvt",
    #                                          work_dir=ligand_nvt_unbound_dir) 
    # nvt_sp = sp.run(["sh", f"{ligand_nvt_unbound_dir}/{ligand_nvt_script}"], capture_output=True, text=True)
    # write_log_file(ligand_nvt_unbound_dir, "nvt", nvt_sp)        

    # ligand_nvt = bss.IO.readMolecules([f"{ligand_nvt_unbound_dir}/nvt.gro",
    #                                    f"{ligand_nvt_unbound_dir}/nvt.top"])
    # ligand_r_npt_protocol = bss.Protocol.Equilibration(runtime=npt_runtime,
    #                                                    pressure=1*atm,
    #                                                    temperature=300*kelvin,
    #                                                    restraint="heavy")
    # change_barostat(ligand_nvt, ligand_r_npt_protocol, ligand_r_npt_unbound_dir, "r_npt")
    # r_npt_sp = sp.run(["sh", f"{ligand_r_npt_unbound_dir}/{ligand_r_npt_script}"], capture_output=True, text=True)
    # write_log_file(ligand_r_npt_unbound_dir, "r_npt", r_npt_sp)     

    # restrained_npt_ligand = bss.IO.readMolecules([f"{ligand_r_npt_unbound_dir}/r_npt.gro",
    #                                               f"{ligand_r_npt_unbound_dir}/r_npt.top"])
    # ligand_npt_protocol = bss.Protocol.Equilibration(runtime=npt_runtime,
    #                                                  pressure=1*atm,
    #                                                  temperature=300*kelvin)    
    # change_barostat(restrained_npt_ligand, ligand_npt_protocol, ligand_npt_unbound_dir, "npt")
    # npt_sp = sp.run(["sh", f"{ligand_npt_unbound_dir}/{ligand_npt_script}"], capture_output=True, text=True)
    # write_log_file(ligand_npt_unbound_dir, "npt", npt_sp)   