import warnings
warnings.filterwarnings("ignore")
import BioSimSpace as bss
import os
import nglview as nv
import subprocess as sp
import argparse
import glob
import math


def check_positive(input):
    try:
        nsteps = int(input)
        if nsteps <= 0:
            raise argparse.ArgumentTypeError(f"{nsteps} is an invalid integer net charge")
    except ValueError:
        print("Error: number of minimisation steps should be integer.")
    except argparse.ArgumentTypeError as message:
        print(message)
    return nsteps


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
        elif restrained and previous == "min": 
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
    return directory, path_exists


parser = argparse.ArgumentParser(description="perform minimisation and equilibration for AFE ligands")

parser.add_argument("system",
                     type=str,
                     help="system name; this is used to find the folder containing input files")

parser.add_argument("-g",
                    "--gpu-id",
                    type=str,
                    help="GPU device ID availble for mdrun",
                    choices=["0", "1", "2", "3"],
                    default="0")

parser.add_argument("-m",
                    "--minimisation-steps",
                    type=check_positive,
                    default=250)

parser.add_argument("-s",
                    "--short-nvt-runtime",
                    help="runtime (in ps) for short NVT equilibration",
                    type=check_positive,
                    default=5)

parser.add_argument("-nvt",
                    "--nvt-runtime",
                    type=check_positive,
                    help="runtime (in ps) for NVT equilibration",
                    default=50)

parser.add_argument("-npt",
                    "--npt-runtime",
                    type=check_positive,
                    help="runtime (in ps) for NPT equilibration",
                    default=200)


picosecond = bss.Units.Time.picosecond
kelvin = bss.Units.Temperature.kelvin
atm = bss.Units.Pressure.atm

arguments = parser.parse_args()
system_name = arguments.system
minimisation_steps = arguments.minimisation_steps
runtime_short_nvt = arguments.short_nvt_runtime * picosecond
runtime_nvt = arguments.nvt_runtime * picosecond
runtime_npt = arguments.npt_runtime * picosecond
gpu_id = arguments.gpu_id

full_path = os.getcwd() + "/"
if "scripts" in full_path:
    full_path = full_path.replace("/scripts/", "/")

afe_folder_path = full_path + system_name + "/afe/"
ligands_file = afe_folder_path + "ligands.dat"
protocol_file = afe_folder_path + "protocol.dat"
network_file = afe_folder_path + "network.dat"
ligand_path = full_path + system_name + "/inputs/ligands/"
protein_path = full_path + system_name + "/inputs/protein/"
# DO THIS LATER
# if not os.path.isfile(input_file):
#     print(f"The input file {input_file} does not exist")
#     sys.exit()
with open(ligands_file) as file:
    ligand_lines = file.readlines()
ligand_names = [line.rstrip() for line in ligand_lines]
n_ligands = len(ligand_names)

print("Processing unbound stage.")
for i in range(n_ligands):
    ligand_number = ligand_names[i].split("_")[1]
    ligand_work_dir, _ = create_dirs(f"{full_path}/{system_name}/equilibration/unbound/ligand_{ligand_number}")
    ligand_min_unbound_dir, do_unbound = create_dirs(f"{ligand_work_dir}/min" )

    # check if min.log already exists # MAYBE DO THIS
        # if yes, skip ligand
        # if no, do equilibration

    ligand_r_nvt_unbound_dir, _ = create_dirs(f"{ligand_work_dir}/r_nvt")
    ligand_nvt_unbound_dir, _ = create_dirs(f"{ligand_work_dir}/nvt")  
    ligand_r_npt_unbound_dir, _ = create_dirs(f"{ligand_work_dir}/r_npt")
    ligand_npt_unbound_dir, _ = create_dirs(f"{ligand_work_dir}/npt")   

    if do_unbound:
        ligand_min_script = write_script(ligand_min_unbound_dir, "min", gpu_id)
        ligand_r_nvt_script = write_script(ligand_r_nvt_unbound_dir, "r_nvt", gpu_id, restrained=True, previous="min")
        ligand_nvt_script = write_script(ligand_nvt_unbound_dir, "nvt", gpu_id, previous="r_nvt", )
        ligand_r_npt_script = write_script(ligand_r_npt_unbound_dir, "r_npt", gpu_id, restrained=True, previous="nvt")
        ligand_npt_script = write_script(ligand_npt_unbound_dir, "npt", gpu_id, previous="r_npt")

        solvated_ligand_files = glob.glob(ligand_path + f"ligand_{ligand_number}_solvated.*")

        solvated_ligand = bss.IO.readMolecules(solvated_ligand_files)

        ligand_minimisation_protocol = bss.Protocol.Minimisation(steps=minimisation_steps)
        ligand_minimisation_process = bss.Process.Gromacs(solvated_ligand,
                                                          ligand_minimisation_protocol,
                                                          name="min",
                                                          work_dir=ligand_min_unbound_dir)
        edit_mdp_options(ligand_min_unbound_dir, "min", {"emstep": 0.001})
        min_sp = sp.run(["sh", f"{ligand_min_unbound_dir}/{ligand_min_script}"], capture_output=True, text=True)
        write_log_file(ligand_min_unbound_dir, ligand_number, "min", min_sp)

        minimised_ligand = bss.IO.readMolecules([f"{ligand_min_unbound_dir}/min.gro",
                                                f"{ligand_min_unbound_dir}/min.top"])

        with open(f"{ligand_min_unbound_dir}/min.log", "r") as file:
            min_log_lines = file.readlines()
        
        max_force_line = [line for line in min_log_lines if "Maximum force" in line][0]
        max_force = float(max_force_line.split()[3])
        max_force_magnitude = math.floor(math.log10(max_force))
        if max_force_magnitude >= 5:
            print(f"Running another minimisation for ligand {ligand_number}")
            new_protocol = bss.Protocol.Minimisation(steps=minimisation_steps*10)
            new_process = bss.Process.Gromacs(minimised_ligand, 
                                              new_protocol,
                                              name="min_2",
                                              work_dir=ligand_min_unbound_dir)

            edit_mdp_options(ligand_min_unbound_dir, "min_2", {"nsteps": 20000})

            new_script = write_script(ligand_min_unbound_dir, "min_2", gpu_id="2")
            new_min_sp = sp.run(["sh", f"{ligand_min_unbound_dir}/{new_script}"], capture_output=True, text=True)
            write_log_file(ligand_min_unbound_dir, ligand_number, "min_2", new_min_sp)

            minimised_ligand = bss.IO.readMolecules([f"{ligand_min_unbound_dir}/min_2.gro",
                                                    f"{ligand_min_unbound_dir}/min.top"])
            ligand_r_nvt_script = write_script(ligand_r_nvt_unbound_dir, "r_nvt", gpu_id, restrained=True, previous="min_2")
                                                    
        ligand_r_nvt_protocol = bss.Protocol.Equilibration(runtime=runtime_short_nvt,
                                                        temperature_start=0*kelvin,
                                                        temperature_end=300*kelvin,
                                                        restraint="all")
        ligand_r_nvt_process = bss.Process.Gromacs(minimised_ligand,
                                                ligand_r_nvt_protocol,
                                                name="r_nvt",
                                                work_dir=ligand_r_nvt_unbound_dir)

        edit_mdp_options(ligand_r_nvt_unbound_dir, "r_nvt", {"dt": 0.0005})
        r_nvt_sp = sp.run(["sh", f"{ligand_r_nvt_unbound_dir}/{ligand_r_nvt_script}"], capture_output=True, text=True)
        write_log_file(ligand_r_nvt_unbound_dir, ligand_number, "r_nvt", r_nvt_sp)  

        restrained_nvt_ligand = bss.IO.readMolecules([f"{ligand_r_nvt_unbound_dir}/r_nvt.gro",
                                                    f"{ligand_r_nvt_unbound_dir}/r_nvt.top"])
        ligand_nvt_protocol = bss.Protocol.Equilibration(runtime=runtime_nvt,
                                                        temperature=300*kelvin)
        ligand_nvt_process = bss.Process.Gromacs(restrained_nvt_ligand,
                                                ligand_nvt_protocol,
                                                name="nvt",
                                                work_dir=ligand_nvt_unbound_dir) 
        nvt_sp = sp.run(["sh", f"{ligand_nvt_unbound_dir}/{ligand_nvt_script}"], capture_output=True, text=True)
        write_log_file(ligand_nvt_unbound_dir, ligand_number, "nvt", nvt_sp)        

        ligand_nvt = bss.IO.readMolecules([f"{ligand_nvt_unbound_dir}/nvt.gro",
                                        f"{ligand_nvt_unbound_dir}/nvt.top"])
        ligand_r_npt_protocol = bss.Protocol.Equilibration(runtime=runtime_npt,
                                                        pressure=1*atm,
                                                        temperature=300*kelvin,
                                                        restraint="heavy")
        ligand_r_npt_process =  bss.Process.Gromacs(ligand_nvt,
                                                    ligand_r_npt_protocol,
                                                    name="r_npt",
                                                    work_dir=ligand_r_npt_unbound_dir)
        edit_mdp_options(ligand_r_npt_unbound_dir, "r_npt", {"pcoupl": "C-rescale"})
        # change_barostat(ligand_nvt, ligand_r_npt_protocol, ligand_r_npt_unbound_dir, "r_npt")
        r_npt_sp = sp.run(["sh", f"{ligand_r_npt_unbound_dir}/{ligand_r_npt_script}"], capture_output=True, text=True)
        write_log_file(ligand_r_npt_unbound_dir, ligand_number, "r_npt", r_npt_sp)     

        restrained_npt_ligand = bss.IO.readMolecules([f"{ligand_r_npt_unbound_dir}/r_npt.gro",
                                                    f"{ligand_r_npt_unbound_dir}/r_npt.top"])
        ligand_npt_protocol = bss.Protocol.Equilibration(runtime=runtime_npt,
                                                        pressure=1*atm,
                                                        temperature=300*kelvin)    
        ligand_r_npt_process =  bss.Process.Gromacs(restrained_npt_ligand,
                                                    ligand_npt_protocol,
                                                    name="npt",
                                                    work_dir=ligand_npt_unbound_dir)
        # change_barostat(restrained_npt_ligand, ligand_npt_protocol, ligand_npt_unbound_dir, "npt")
        edit_mdp_options(ligand_npt_unbound_dir, "npt", {"pcoupl": "C-rescale"})
        npt_sp = sp.run(["sh", f"{ligand_npt_unbound_dir}/{ligand_npt_script}"], capture_output=True, text=True)
        write_log_file(ligand_npt_unbound_dir, ligand_number, "npt", npt_sp)     
    else:
        print(f"Unbound stage already computed for ligand {ligand_number}")

print("Processing bound stage.")
for i in range(n_ligands):
    ligand_number = ligand_names[i].split("_")[1]
    ligand_work_dir, _ = create_dirs(f"{full_path}/{system_name}/equilibration/bound/system_{ligand_number}/")
    system_min_bound_dir, _ = create_dirs(f"{ligand_work_dir}/min") 
    system_r_nvt_bound_dir, _ = create_dirs(f"{ligand_work_dir}/r_nvt")
    system_bb_r_nvt_bound_dir, _ = create_dirs(f"{ligand_work_dir}/bb_r_nvt")  
    system_r_npt_bound_dir, _ = create_dirs(f"{ligand_work_dir}/r_npt")
    system_npt_bound_dir, _ = create_dirs(f"{ligand_work_dir}/npt")
    system_nvt_bound_dir, _ = create_dirs(f"{ligand_work_dir}/nvt")  

    system_min_script = write_script(system_min_bound_dir, "min", gpu_id)
    system_r_nvt_script = write_script(system_r_nvt_bound_dir, "r_nvt", gpu_id, restrained=True, previous="min")
    system_bb_r_nvt_script = write_script(system_bb_r_nvt_bound_dir, "bb_r_nvt", gpu_id, restrained=True, previous="r_nvt")
    system_nvt_script = write_script(system_nvt_bound_dir, "nvt", gpu_id, previous="bb_r_nvt")
    system_r_npt_script = write_script(system_r_npt_bound_dir, "r_npt", gpu_id, restrained=True, previous="nvt")
    system_npt_script = write_script(system_npt_bound_dir, "npt", gpu_id, previous="r_npt")

    solvated_system_files = glob.glob(protein_path + f"system_{ligand_number}_solvated.*")

    solvated_system = bss.IO.readMolecules(solvated_system_files)

    system_minimisation_protocol = bss.Protocol.Minimisation(steps=minimisation_steps*10)
    system_minimisation_process = bss.Process.Gromacs(solvated_system,
                                                      system_minimisation_protocol,
                                                      name="min",
                                                      work_dir=system_min_bound_dir)
    edit_mdp_options(system_min_bound_dir, "min", {"emstep": 0.001})
    min_sp = sp.run(["sh", f"{system_min_bound_dir}/{system_min_script}"], capture_output=True, text=True)
    write_log_file(system_min_bound_dir, ligand_number, "min", min_sp)

    minimised_system = bss.IO.readMolecules([f"{system_min_bound_dir}/min.gro",
                                             f"{system_min_bound_dir}/min.top"])

    with open(f"{system_min_bound_dir}/min.log", "r") as file:
        min_log_lines = file.readlines()
    
    max_force_line = [line for line in min_log_lines if "Maximum force" in line][0]
    max_force = float(max_force_line.split()[3])
    max_force_magnitude = math.floor(math.log10(max_force))
    if max_force_magnitude >= 5:
        print(f"Running another minimisation for ligand {ligand_number}")
        new_protocol = bss.Protocol.Minimisation(steps=minimisation_steps*10)
        new_process = bss.Process.Gromacs(minimised_system,
                                          new_protocol, 
                                          name="min_2",
                                          work_dir=system_min_bound_dir)
        
        edit_mdp_options(system_min_bound_dir, "min_2", {"nsteps": 20000})

        new_script = write_script(system_min_bound_dir, "min_2", gpu_id="2")
        new_min_sp = sp.run(["sh", f"{system_min_bound_dir}/{new_script}"], capture_output=True, text=True)
        write_log_file(system_min_bound_dir, ligand_number, "min_2", new_min_sp)

        minimised_system = bss.IO.readMolecules([f"{system_min_bound_dir}/min_2.gro",
                                                 f"{system_min_bound_dir}/min.top"])
        system_r_nvt_script = write_script(system_min_bound_dir, "r_nvt", gpu_id, restrained=True, previous="min_2")
                                                
    system_r_nvt_protocol = bss.Protocol.Equilibration(runtime=runtime_short_nvt,
                                                       temperature_start=0*kelvin,
                                                       temperature_end=300*kelvin,
                                                       restraint="all")
    system_r_nvt_process = bss.Process.Gromacs(minimised_system,
                                               system_r_nvt_protocol,
                                               name="r_nvt",
                                               work_dir=system_r_nvt_bound_dir)
    edit_mdp_options(system_r_nvt_bound_dir, "r_nvt", {"dt": 0.0005})
    r_nvt_sp = sp.run(["sh", f"{system_r_nvt_bound_dir}/{system_r_nvt_script}"], capture_output=True, text=True)
    write_log_file(system_r_nvt_bound_dir, ligand_number, "r_nvt", r_nvt_sp)  

    restrained_nvt_system = bss.IO.readMolecules([f"{system_r_nvt_bound_dir}/r_nvt.gro",
                                                  f"{system_r_nvt_bound_dir}/r_nvt.top"])

    system_bb_restrained_protocol = bss.Protocol.Equilibration(runtime=runtime_nvt,
                                                               temperature=300*kelvin,
                                                               restraint="backbone") 
    system_bb_restrained_process = bss.Process.Gromacs(restrained_nvt_system,
                                                       system_bb_restrained_protocol,
                                                       name="bb_r_nvt",
                                                       work_dir=system_bb_r_nvt_bound_dir)
                                                       
    bb_r_nvt_sp = sp.run(["sh", f"{system_bb_r_nvt_bound_dir}/{system_bb_r_nvt_script}"], capture_output=True, text=True)
    write_log_file(system_bb_r_nvt_bound_dir, ligand_number, "bb_r_nvt", bb_r_nvt_sp)

    bb_restrained_system = bss.IO.readMolecules([f"{system_bb_r_nvt_bound_dir}/bb_r_nvt.gro",
                                                 f"{system_bb_r_nvt_bound_dir}/bb_r_nvt.top"])
    
    system_nvt_protocol = bss.Protocol.Equilibration(runtime=runtime_nvt,
                                                     temperature=300*kelvin)
    system_nvt_process = bss.Process.Gromacs(bb_restrained_system,
                                             system_nvt_protocol,
                                             name="nvt",
                                             work_dir=system_nvt_bound_dir) 

    nvt_sp = sp.run(["sh", f"{system_nvt_bound_dir}/{system_nvt_script}"], capture_output=True, text=True)
    write_log_file(system_nvt_bound_dir, ligand_number, "nvt", nvt_sp)        

    nvt_system = bss.IO.readMolecules([f"{system_nvt_bound_dir}/nvt.gro",
                                       f"{system_nvt_bound_dir}/nvt.top"])
    
    system_r_npt_protocol = bss.Protocol.Equilibration(runtime=runtime_npt,
                                                       pressure=1*atm,
                                                       temperature=300*kelvin,
                                                       restraint="heavy")
    system_r_npt_process =  bss.Process.Gromacs(nvt_system,
                                                system_r_npt_protocol,
                                                name="r_npt",
                                                work_dir=system_r_npt_bound_dir)
    edit_mdp_options(system_r_npt_bound_dir, "r_npt", {"pcoupl": "C-rescale"})
#     change_barostat(nvt_system, system_r_npt_protocol, system_r_npt_bound_dir, "r_npt")
    r_npt_sp = sp.run(["sh", f"{system_r_npt_bound_dir}/{system_r_npt_script}"], capture_output=True, text=True)
    write_log_file(system_r_npt_bound_dir, ligand_number, "r_npt", r_npt_sp)     

    restrained_npt_system = bss.IO.readMolecules([f"{system_r_npt_bound_dir}/r_npt.gro",
                                                  f"{system_r_npt_bound_dir}/r_npt.top"])
    system_npt_protocol = bss.Protocol.Equilibration(runtime=runtime_npt,
                                                     pressure=1*atm,
                                                     temperature=300*kelvin)    
    system_r_npt_process =  bss.Process.Gromacs(restrained_npt_system,
                                                system_npt_protocol,
                                                name="npt",
                                                work_dir=system_npt_bound_dir)

    edit_mdp_options(system_npt_bound_dir, "npt", {"pcoupl": "C-rescale"})
#     change_barostat(restrained_npt_system, system_npt_protocol, system_npt_bound_dir, "npt")
    npt_sp = sp.run(["sh", f"{system_npt_bound_dir}/{system_npt_script}"], capture_output=True, text=True)
    write_log_file(system_npt_bound_dir, ligand_number, "npt", npt_sp)   