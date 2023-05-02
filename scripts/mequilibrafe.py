import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
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

parser.add_argument("path",
                     type=str,
                     help="path; this is used to find the folder containing input files")

parser.add_argument("ligand_number",
                    type=str,
                    help="number of ligand to minimise&equilibrate for AFE")

parser.add_argument("-m",
                    "--minimisation-steps",
                    type=fn.check_positive_integer,
                    default=500)

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
path = arguments.path
ligand_number = arguments.ligand_number

minimisation_steps = arguments.minimisation_steps
runtime_short_nvt = arguments.short_nvt_runtime * picosecond
runtime_nvt = arguments.nvt_runtime * picosecond
runtime_npt = arguments.npt_runtime * picosecond

full_path = os.getenv("HOME") + "/projects/metalloenzymes/"

afe_folder_path = full_path + system_name + "/afe/"
afe_folder_path = path + "/afe/"
protocol_file = afe_folder_path + "protocol.dat"
network_file = afe_folder_path + "network.dat"
ligand_path = full_path + system_name + "/inputs/ligands/"
ligand_path = path + "/inputs/ligands/"
print(ligand_path)
protein_path = full_path + system_name + "/inputs/protein/"
protein_path = path + "/inputs/protein/"
print(protein_path)

#TODO  
# if not os.path.isfile(input_file):
#     print(f"The input file {input_file} does not exist")
#     sys.exit()

print("Processing unbound stage.")

ligand_work_dir = fn.create_dirs(f"{path}/equilibration/unbound/ligand_{ligand_number}")
ligand_min_unbound_dir = fn.create_dirs(f"{ligand_work_dir}/min" )
ligand_r_nvt_unbound_dir = fn.create_dirs(f"{ligand_work_dir}/r_nvt")
ligand_nvt_unbound_dir = fn.create_dirs(f"{ligand_work_dir}/nvt")  
ligand_r_npt_unbound_dir = fn.create_dirs(f"{ligand_work_dir}/r_npt")
ligand_npt_unbound_dir = fn.create_dirs(f"{ligand_work_dir}/npt")   

solvated_ligand_files = glob.glob(ligand_path + f"ligand_{ligand_number}_solvated.*")
solvated_ligand = bss.IO.readMolecules(solvated_ligand_files)
print("Minimisation")
ligand_minimisation_protocol = bss.Protocol.Minimisation(steps=minimisation_steps)
minimised_ligand = fn.run_process(system=solvated_ligand, 
                                  protocol=ligand_minimisation_protocol, 
                                  process_name="min", 
                                  working_directory=ligand_min_unbound_dir,
                                  configuration=["emstep = 0.001", "emtol = 1000"]
                                  )
print("R-NVT")
ligand_r_nvt_protocol = bss.Protocol.Equilibration(runtime=runtime_short_nvt,
                                                   temperature_start=0*kelvin,
                                                   temperature_end=300*kelvin,
                                                   restraint="all")

restrained_nvt_ligand = fn.run_process(system=minimised_ligand,
                                       protocol=ligand_r_nvt_protocol,
                                       process_name="r_nvt",
                                       working_directory=ligand_r_nvt_unbound_dir,
                                       configuration=["dt = 0.0005"]
                                       )
print("NVT")
ligand_nvt_protocol = bss.Protocol.Equilibration(runtime=runtime_nvt, temperature=300*kelvin)
ligand_nvt = fn.run_process(system=restrained_nvt_ligand,
                            protocol=ligand_nvt_protocol,
                            process_name="nvt",
                            working_directory=ligand_nvt_unbound_dir)
print("R-NPT")
ligand_r_npt_protocol = bss.Protocol.Equilibration(runtime=runtime_npt,
                                                   pressure=1*atm,
                                                   temperature=300*kelvin,
                                                   restraint="heavy")
restrained_npt_ligand = fn.run_process(system=ligand_nvt,
                                       protocol=ligand_r_npt_protocol,
                                       process_name="r_npt",
                                       working_directory=ligand_r_npt_unbound_dir)

print("NPT")
ligand_npt_protocol = bss.Protocol.Equilibration(runtime=runtime_npt,
                                                 pressure=1*atm,
                                                 temperature=300*kelvin)   
equilibrated_ligand = fn.run_process(system=restrained_npt_ligand,
                                     protocol=ligand_npt_protocol,
                                     process_name="npt",
                                     working_directory=ligand_npt_unbound_dir)

bss.IO.saveMolecules(filebase=ligand_npt_unbound_dir + f"/ligand_{ligand_number}",
                     system=equilibrated_ligand,
                     fileformat=["PRM7", "RST7", "GroTop", "Gro87"])

print("Processing bound stage.")

ligand_work_dir = fn.create_dirs(f"{path}/equilibration/bound/system_{ligand_number}/")
system_min_bound_dir = fn.create_dirs(f"{ligand_work_dir}/min") 
system_r_nvt_bound_dir = fn.create_dirs(f"{ligand_work_dir}/r_nvt")
system_bb_r_nvt_bound_dir = fn.create_dirs(f"{ligand_work_dir}/bb_r_nvt")  
system_r_npt_bound_dir = fn.create_dirs(f"{ligand_work_dir}/r_npt")
system_npt_bound_dir = fn.create_dirs(f"{ligand_work_dir}/npt")
system_nvt_bound_dir = fn.create_dirs(f"{ligand_work_dir}/nvt")  

solvated_system_files = glob.glob(protein_path + f"system_{ligand_number}_solvated.*")
solvated_system = bss.IO.readMolecules(solvated_system_files)

system_minimisation_protocol = bss.Protocol.Minimisation(steps=minimisation_steps*10)
print("Minimisation")
minimised_system = fn.run_process(system=solvated_system,
                                  protocol=system_minimisation_protocol,
                                  process_name="min",
                                  working_directory=system_min_bound_dir,
                                  configuration=["emstep = 0.001", "emtol = 1000"]
                                  )

system_r_nvt_protocol = bss.Protocol.Equilibration(runtime=runtime_short_nvt,
                                                   temperature_start=0*kelvin,
                                                   temperature_end=300*kelvin,
                                                   restraint="all")
print("R-NVT")
restrained_nvt_system = fn.run_process(system=minimised_system,
                                       protocol=system_r_nvt_protocol,
                                       process_name="r_nvt",
                                       working_directory=system_r_nvt_bound_dir,
                                       configuration=["dt = 0.0005"]
                                       )

system_bb_restrained_protocol = bss.Protocol.Equilibration(runtime=runtime_nvt,
                                                           temperature=300*kelvin,
                                                           restraint="backbone") 
print("BB-R-NVT")
bb_restrained_system = fn.run_process(system=restrained_nvt_system,
                                      protocol=system_bb_restrained_protocol,
                                      process_name="bb_r_nvt",
                                      working_directory=system_bb_r_nvt_bound_dir)

system_nvt_protocol = bss.Protocol.Equilibration(runtime=runtime_nvt, temperature=300*kelvin)
print("NVT")
nvt_system = fn.run_process(system=bb_restrained_system,
                            protocol=system_nvt_protocol,
                            process_name="nvt",
                            working_directory=system_nvt_bound_dir)

system_r_npt_protocol = bss.Protocol.Equilibration(runtime=runtime_npt,
                                                   pressure=1*atm,
                                                   temperature=300*kelvin,
                                                   restraint="heavy")
print("R-NPT")
restrained_npt_system = fn.run_process(system=nvt_system,
                                       protocol=system_r_npt_protocol,
                                       process_name="r_npt",
                                       working_directory=system_r_npt_bound_dir)                                                   

system_npt_protocol = bss.Protocol.Equilibration(runtime=runtime_npt,
                                                 pressure=1*atm,
                                                 temperature=300*kelvin)    
print("NPT")                                  
equilibrated_system = fn.run_process(system=restrained_npt_system,
                                     protocol=system_npt_protocol,
                                     process_name="npt",
                                     working_directory=system_npt_bound_dir)

bss.IO.saveMolecules(filebase=system_npt_bound_dir + f"/system_{ligand_number}",
                     system=equilibrated_system,
                     fileformat=["PRM7", "RST7", "GroTop", "Gro87"])
