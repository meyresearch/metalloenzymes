import warnings
warnings.filterwarnings("ignore")
import argparse
import BioSimSpace as bss
import os
import subprocess as sp
import sys


def check_positive(input):
    try:
        nsteps = int(input)
        if nsteps <= 0:
            raise argparse.ArgumentTypeError(f"{nsteps} is an invalid integer value")
    except ValueError:
        print("Error: number of minimisation steps and/or runtime should be integer.")
    except argparse.ArgumentTypeError as message:
        print(message)
    return nsteps


def check_integer(input):
    try:
        charge_float = float(input)
        if charge_float < 0:
            charge = - int(abs(charge_float))
        else:
            charge = int(charge_float)
    except ValueError:
        print("Error: charge should be integer.")
    return charge

# DEPRECATED
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
        bss.Process.Gromacs(system,
                            protocol,
                            name=process_name,
                            work_dir=work_directory)
    except RuntimeError:

        with open(f"{work_directory}/{process_name}.mdp", "r") as mdp:
            lines = mdp.read()

        new_lines = lines.replace("pcoupl = berendsen", "pcoupl = C-rescale")
        with open(f"{work_directory}/{process_name}.mdp", "w") as mdp:
            mdp.write(new_lines)

# COMBINE THIS AND edit_mdp_options TO ONE FUNCTION IN FUTURE
def lambda_mdps(afe_mdp_file: str, save_directory: str, process_name: str, options: dict):
    with open(afe_mdp_file) as afe:
        afe_lines = afe.readlines()
    
    get_mdp_options = lambda index: [line.split("=")[index].lstrip().rstrip().strip("\n") for line in afe_lines]
    mdp_keys = get_mdp_options(0)
    mdp_values = get_mdp_options(-1)

    change_indices = [mdp_keys.index(key) for key in options.keys() if key in mdp_keys]

    for i in change_indices:
        for key, value in options.items():
            if key in mdp_keys:
                mdp_values[i] = value

    for key, value in options.items():
        if key not in mdp_keys:
            formatted_key = key
            formatted_value = str(value)
            mdp_keys.append(formatted_key)
            mdp_values.append(formatted_value) 

    with open(f"{save_directory}/{process_name}.mdp", "w") as mdp:
        for i in range(len(mdp_keys)):
            mdp.write(f"{mdp_keys[i]} = {mdp_values[i]}\n")


def edit_mdp_options(work_dir: str, process_name: str, options: dict):
    """
    Write an mdp files
    """
    with open(f"{work_dir}/{process_name}.mdp", "r") as mdp:
        mdp_lines = mdp.readlines()

    get_mdp_options = lambda index: [line.split("=")[index].lstrip().rstrip().strip("\n") for line in mdp_lines]
    mdp_keys = get_mdp_options(0)
    mdp_values = get_mdp_options(-1)

    change_indices = [mdp_keys.index(key) for key in options.keys() if key in mdp_keys]

    for i in change_indices:
        for key, value in options.items():
            if key in mdp_keys:
                mdp_values[i] = value
    # print(mdp_values)
    for key, value in options.items():
        if key not in mdp_keys:
            formatted_key = key
            formatted_value = str(value)
            mdp_keys.append(formatted_key)
            mdp_values.append(formatted_value) 
    
    with open(f"{work_dir}/{process_name}.mdp", "w") as mdp:
        for i in range(len(mdp_keys)):
            mdp.write(f"{mdp_keys[i]} = {mdp_values[i]}\n")


def write_script(work_dir: str, process_name: str, gpu_id: str, restrained=False, previous="") -> str:
    script_name = work_dir + "/" + process_name + ".sh"
    with open(script_name, "w") as script:
        script.write("#!/bin/bash \n")
        script.write(f"cd {work_dir} \n")
        if not restrained and previous == "":
            script.write(f"gmx grompp -f {process_name}.mdp -c {process_name}.gro -p {process_name}.top -o {process_name}.tpr \n")
        elif restrained and previous == "min" or previous == "min_2": 
            script.write(f"gmx grompp -f {process_name}.mdp -c ../{previous}/{previous}.gro -r ../{previous}/{previous}.gro -p {process_name}.top -o {process_name}.tpr \n")
        elif not restrained and previous != "":
            script.write(f"gmx grompp -f {process_name}.mdp -c ../{previous}/{previous}.gro -p {process_name}.top -t ../{previous}/{previous}.cpt -o {process_name}.tpr \n")
        elif restrained and previous != "":
            script.write(f"gmx grompp -f {process_name}.mdp -c ../{previous}/{previous}.gro -r ../{previous}/{previous}.gro -p {process_name}.top -t ../{previous}/{previous}.cpt -o {process_name}.tpr \n") 
        script.write(f"gmx mdrun -v -deffnm {process_name} -nt 1 -nb gpu -gpu_id {gpu_id}")
    os.system(f"chmod +x {script_name}")
    executable = script_name.split("/")[-1]
    return executable


def write_log_file(work_dir: str, ligand_number: str, process_name: str, subprocess: sp.CompletedProcess) -> None:
    """
    Write custom log file for gromacs processes.
    """
    try:
        if "error" in subprocess.stderr or "error" in subprocess.stdout: # not if on the same line as GROMACS Reminds you!
            print(f"There was an error with {process_name}.sh. Please check the log file for ligand {ligand_number}.")
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


def create_complex(protein_path: str, protein_file: str, water_file: str, output=None) -> str:
    protein = bss.IO.readMolecules(protein_path+protein_file)
    if water_file is not None and output == "":
        output = "protein_complex"
        xtal_water = bss.IO.readMolecules(water_file)
        protein_complex = protein + xtal_water
        bss.IO.saveMolecules(protein_path+output, protein_complex, fileformat="pdb")
    elif water_file is not None and output != "":
        xtal_water = bss.IO.readMolecules(water_file)
        protein_complex = protein + xtal_water
        bss.IO.saveMolecules(protein_path+output, protein_complex, fileformat="pdb")
    elif water_file is None and output == "":
        output = "protein_complex"
        protein_complex = protein
        bss.IO.saveMolecules(protein_path+output, protein_complex, fileformat="pdb")
    elif water_file is None and output != "":
        protein_complex = protein
        bss.IO.saveMolecules(protein_path+output, protein_complex, fileformat="pdb")
    return output


def is_file(file: str) -> None:
    if not os.path.isfile(file):
        print(f"The file {file} does not exist")
        sys.exit()