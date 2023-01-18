import warnings
warnings.filterwarnings("ignore")
import BioSimSpace as bss
import os
import nglview as nv
import subprocess as sp
import argparse
import glob
import math
import functions as fn


parser = argparse.ArgumentParser(description="perform minimisation and equilibration for AFE ligands")

parser.add_argument("system",
                     type=str,
                     help="system name; this is used to find the folder containing input files")

# parser.add_argument("ligand_number",
#                     type=str,
#                     help="number of ligand to minimise&equilibrate for AFE")

parser.add_argument("-g",
                    "--gpu-id",
                    type=str,
                    help="GPU device ID availble for mdrun",
                    choices=["0", "1", "2", "3"],
                    default="0")

parser.add_argument("-m",
                    "--minimisation-steps",
                    type=fn.check_positive_integer,
                    default=250)

parser.add_argument("-s",
                    "--short-nvt-runtime",
                    help="runtime (in ps) for short NVT equilibration",
                    type=fn.check_positive_float,
                    default=5)

parser.add_argument("-nvt",
                    "--nvt-runtime",
                    type=fn.check_positive_float,
                    help="runtime (in ps) for NVT equilibration",
                    default=50)

parser.add_argument("-npt",
                    "--npt-runtime",
                    type=fn.check_positive_float,
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
    ligand_work_dir = fn.create_dirs(f"{full_path}/{system_name}/equilibration/unbound/ligand_{ligand_number}")
    ligand_min_unbound_dir = fn.create_dirs(f"{ligand_work_dir}/min" )

    ligand_r_nvt_unbound_dir = fn.create_dirs(f"{ligand_work_dir}/r_nvt")
    ligand_nvt_unbound_dir = fn.create_dirs(f"{ligand_work_dir}/nvt")  
    ligand_r_npt_unbound_dir = fn.create_dirs(f"{ligand_work_dir}/r_npt")
    ligand_npt_unbound_dir = fn.create_dirs(f"{ligand_work_dir}/npt")   


    ligand_min_script = fn.write_script(ligand_min_unbound_dir, "min", gpu_id)
    ligand_r_nvt_script = fn.write_script(ligand_r_nvt_unbound_dir, "r_nvt", gpu_id, restrained=True, previous="min")
    ligand_nvt_script = fn.write_script(ligand_nvt_unbound_dir, "nvt", gpu_id, previous="r_nvt", )
    ligand_r_npt_script = fn.write_script(ligand_r_npt_unbound_dir, "r_npt", gpu_id, restrained=True, previous="nvt")
    ligand_npt_script = fn.write_script(ligand_npt_unbound_dir, "npt", gpu_id, previous="r_npt")

    solvated_ligand_files = glob.glob(ligand_path + f"ligand_{ligand_number}_solvated.*")

    solvated_ligand = bss.IO.readMolecules(solvated_ligand_files)

    ligand_minimisation_protocol = bss.Protocol.Minimisation(steps=minimisation_steps)
    ligand_minimisation_process = bss.Process.Gromacs(solvated_ligand,
                                                        ligand_minimisation_protocol,
                                                        name="min",
                                                        work_dir=ligand_min_unbound_dir)
    fn.edit_mdp_options(ligand_min_unbound_dir, "min", {"emstep": 0.001, "emtol": 1000})
    min_sp = sp.run(["sh", f"{ligand_min_unbound_dir}/{ligand_min_script}"], capture_output=True, text=True)
    fn.write_log_file(ligand_min_unbound_dir, ligand_number, "min", min_sp)

    minimised_ligand = bss.IO.readMolecules([f"{ligand_min_unbound_dir}/min.gro",
                                            f"{ligand_min_unbound_dir}/min.top"])

    with open(f"{ligand_min_unbound_dir}/min.log", "r") as file:
        min_log_lines = file.readlines()

    # max_force_line = [line for line in min_log_lines if "Maximum force" in line][0]
    # max_force = float(max_force_line.split()[3])
    # max_force_magnitude = math.floor(math.log10(max_force))
    # if max_force_magnitude >= 5:
    #     print(f"Running another minimisation for ligand {ligand_number}")
    #     new_directory = fn.create_dirs(f"{ligand_work_dir}/min_2")
    #     new_protocol = bss.Protocol.Minimisation(steps=minimisation_steps*10)
    #     new_process = bss.Process.Gromacs(minimised_ligand, 
    #                                     new_protocol,
    #                                     name="min_2",
    #                                     work_dir=new_directory)

    #     fn.edit_mdp_options(new_directory, "min_2", {"nsteps": 20000})

    #     new_script = fn.write_script(new_directory, "min_2", gpu_id="2")
    #     new_min_sp = sp.run(["sh", f"{new_directory}/{new_script}"], capture_output=True, text=True)
    #     fn.write_log_file(new_directory, ligand_number, "min_2", new_min_sp)

    #     minimised_ligand = bss.IO.readMolecules([f"{new_directory}/min_2.gro",
    #                                             f"{ligand_min_unbound_dir}/min.top"])          
    #     ligand_r_nvt_script = fn.write_script(ligand_r_nvt_unbound_dir, "r_nvt", gpu_id, restrained=True, previous="min_2")
                                                
    ligand_r_nvt_protocol = bss.Protocol.Equilibration(runtime=runtime_short_nvt,
                                                    temperature_start=0*kelvin,
                                                    temperature_end=300*kelvin,
                                                    restraint="all")
    ligand_r_nvt_process = bss.Process.Gromacs(minimised_ligand,
                                            ligand_r_nvt_protocol,
                                            name="r_nvt",
                                            work_dir=ligand_r_nvt_unbound_dir)

    fn.edit_mdp_options(ligand_r_nvt_unbound_dir, "r_nvt", {"dt": 0.0005})
    r_nvt_sp = sp.run(["sh", f"{ligand_r_nvt_unbound_dir}/{ligand_r_nvt_script}"], capture_output=True, text=True)
    fn.write_log_file(ligand_r_nvt_unbound_dir, ligand_number, "r_nvt", r_nvt_sp)  

    restrained_nvt_ligand = bss.IO.readMolecules([f"{ligand_r_nvt_unbound_dir}/r_nvt.gro",
                                                f"{ligand_r_nvt_unbound_dir}/r_nvt.top"])
    ligand_nvt_protocol = bss.Protocol.Equilibration(runtime=runtime_nvt,
                                                    temperature=300*kelvin)
    ligand_nvt_process = bss.Process.Gromacs(restrained_nvt_ligand,
                                            ligand_nvt_protocol,
                                            name="nvt",
                                            work_dir=ligand_nvt_unbound_dir) 
    nvt_sp = sp.run(["sh", f"{ligand_nvt_unbound_dir}/{ligand_nvt_script}"], capture_output=True, text=True)
    fn.write_log_file(ligand_nvt_unbound_dir, ligand_number, "nvt", nvt_sp)        

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
    fn.edit_mdp_options(ligand_r_npt_unbound_dir, "r_npt", {"pcoupl": "C-rescale"})
    r_npt_sp = sp.run(["sh", f"{ligand_r_npt_unbound_dir}/{ligand_r_npt_script}"], capture_output=True, text=True)
    fn.write_log_file(ligand_r_npt_unbound_dir, ligand_number, "r_npt", r_npt_sp)     

    restrained_npt_ligand = bss.IO.readMolecules([f"{ligand_r_npt_unbound_dir}/r_npt.gro",
                                                f"{ligand_r_npt_unbound_dir}/r_npt.top"])
    ligand_npt_protocol = bss.Protocol.Equilibration(runtime=runtime_npt,
                                                    pressure=1*atm,
                                                    temperature=300*kelvin)    
    ligand_r_npt_process =  bss.Process.Gromacs(restrained_npt_ligand,
                                                ligand_npt_protocol,
                                                name="npt",
                                                work_dir=ligand_npt_unbound_dir)
    fn.edit_mdp_options(ligand_npt_unbound_dir, "npt", {"pcoupl": "C-rescale"})
    npt_sp = sp.run(["sh", f"{ligand_npt_unbound_dir}/{ligand_npt_script}"], capture_output=True, text=True)
    fn.write_log_file(ligand_npt_unbound_dir, ligand_number, "npt", npt_sp)     


print("Processing bound stage.")
for i in range(n_ligands):
    ligand_number = ligand_names[i].split("_")[1]
    ligand_work_dir = fn.create_dirs(f"{full_path}/{system_name}/equilibration/bound/system_{ligand_number}/")
    system_min_bound_dir = fn.create_dirs(f"{ligand_work_dir}/min") 
    system_r_nvt_bound_dir = fn.create_dirs(f"{ligand_work_dir}/r_nvt")
    system_bb_r_nvt_bound_dir = fn.create_dirs(f"{ligand_work_dir}/bb_r_nvt")  
    system_r_npt_bound_dir = fn.create_dirs(f"{ligand_work_dir}/r_npt")
    system_npt_bound_dir = fn.create_dirs(f"{ligand_work_dir}/npt")
    system_nvt_bound_dir = fn.create_dirs(f"{ligand_work_dir}/nvt")  

    system_min_script = fn.write_script(system_min_bound_dir, "min", gpu_id)
    system_r_nvt_script = fn.write_script(system_r_nvt_bound_dir, "r_nvt", gpu_id, restrained=True, previous="min")
    system_bb_r_nvt_script = fn.write_script(system_bb_r_nvt_bound_dir, "bb_r_nvt", gpu_id, restrained=True, previous="r_nvt")
    system_nvt_script = fn.write_script(system_nvt_bound_dir, "nvt", gpu_id, previous="bb_r_nvt")
    system_r_npt_script = fn.write_script(system_r_npt_bound_dir, "r_npt", gpu_id, restrained=True, previous="nvt")
    system_npt_script = fn.write_script(system_npt_bound_dir, "npt", gpu_id, previous="r_npt")

    solvated_system_files = glob.glob(protein_path + f"system_{ligand_number}_solvated.*")

    solvated_system = bss.IO.readMolecules(solvated_system_files)

    system_minimisation_protocol = bss.Protocol.Minimisation(steps=minimisation_steps*10)
    system_minimisation_process = bss.Process.Gromacs(solvated_system,
                                                        system_minimisation_protocol,
                                                        name="min",
                                                        work_dir=system_min_bound_dir)
    fn.edit_mdp_options(system_min_bound_dir, "min", {"emstep": 0.001, "emtol": 1000})
    min_sp = sp.run(["sh", f"{system_min_bound_dir}/{system_min_script}"], capture_output=True, text=True)
    fn.write_log_file(system_min_bound_dir, ligand_number, "min", min_sp)

    minimised_system = bss.IO.readMolecules([f"{system_min_bound_dir}/min.gro",
                                                f"{system_min_bound_dir}/min.top"])

    with open(f"{system_min_bound_dir}/min.log", "r") as file:
        min_log_lines = file.readlines()

    max_force_line = [line for line in min_log_lines if "Maximum force" in line][0]
    max_force = float(max_force_line.split()[3])
    max_force_magnitude = math.floor(math.log10(max_force))
    if max_force_magnitude >= 5:
        print(f"Running another minimisation for ligand {ligand_number}")
        new_directory = fn.create_dirs(f"{ligand_work_dir}/min_2")
        new_protocol = bss.Protocol.Minimisation(steps=minimisation_steps*10)
        new_process = bss.Process.Gromacs(minimised_system,
                                            new_protocol, 
                                            name="min_2",
                                            work_dir=new_directory)
        
        fn.edit_mdp_options(new_directory, "min_2", {"nsteps": 20000})

        new_script = fn.write_script(new_directory, "min_2", gpu_id="2")
        new_min_sp = sp.run(["sh", f"{new_directory}/{new_script}"], capture_output=True, text=True)
        fn.write_log_file(new_directory, ligand_number, "min_2", new_min_sp)

        minimised_system = bss.IO.readMolecules([f"{new_directory}/min_2.gro",
                                                f"{system_min_bound_dir}/min.top"])
        system_r_nvt_script = fn.write_script(system_r_nvt_bound_dir, "r_nvt", gpu_id, restrained=True, previous="min_2")
                                                
    system_r_nvt_protocol = bss.Protocol.Equilibration(runtime=runtime_short_nvt,
                                                        temperature_start=0*kelvin,
                                                        temperature_end=300*kelvin,
                                                        restraint="all")
    system_r_nvt_process = bss.Process.Gromacs(minimised_system,
                                                system_r_nvt_protocol,
                                                name="r_nvt",
                                                work_dir=system_r_nvt_bound_dir)
    fn.edit_mdp_options(system_r_nvt_bound_dir, "r_nvt", {"dt": 0.0005})
    r_nvt_sp = sp.run(["sh", f"{system_r_nvt_bound_dir}/{system_r_nvt_script}"], capture_output=True, text=True)
    fn.write_log_file(system_r_nvt_bound_dir, ligand_number, "r_nvt", r_nvt_sp)  

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
    fn.write_log_file(system_bb_r_nvt_bound_dir, ligand_number, "bb_r_nvt", bb_r_nvt_sp)

    bb_restrained_system = bss.IO.readMolecules([f"{system_bb_r_nvt_bound_dir}/bb_r_nvt.gro",
                                                    f"{system_bb_r_nvt_bound_dir}/bb_r_nvt.top"])

    system_nvt_protocol = bss.Protocol.Equilibration(runtime=runtime_nvt,
                                                        temperature=300*kelvin)
    system_nvt_process = bss.Process.Gromacs(bb_restrained_system,
                                                system_nvt_protocol,
                                                name="nvt",
                                                work_dir=system_nvt_bound_dir) 

    nvt_sp = sp.run(["sh", f"{system_nvt_bound_dir}/{system_nvt_script}"], capture_output=True, text=True)
    fn.write_log_file(system_nvt_bound_dir, ligand_number, "nvt", nvt_sp)        

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
    fn.edit_mdp_options(system_r_npt_bound_dir, "r_npt", {"pcoupl": "C-rescale"})
    r_npt_sp = sp.run(["sh", f"{system_r_npt_bound_dir}/{system_r_npt_script}"], capture_output=True, text=True)
    fn.write_log_file(system_r_npt_bound_dir, ligand_number, "r_npt", r_npt_sp)     

    restrained_npt_system = bss.IO.readMolecules([f"{system_r_npt_bound_dir}/r_npt.gro",
                                                    f"{system_r_npt_bound_dir}/r_npt.top"])
    system_npt_protocol = bss.Protocol.Equilibration(runtime=runtime_npt,
                                                        pressure=1*atm,
                                                        temperature=300*kelvin)    
    system_r_npt_process =  bss.Process.Gromacs(restrained_npt_system,
                                                system_npt_protocol,
                                                name="npt",
                                                work_dir=system_npt_bound_dir)

    fn.edit_mdp_options(system_npt_bound_dir, "npt", {"pcoupl": "C-rescale"})
    npt_sp = sp.run(["sh", f"{system_npt_bound_dir}/{system_npt_script}"], capture_output=True, text=True)
    fn.write_log_file(system_npt_bound_dir, ligand_number, "npt", npt_sp)   