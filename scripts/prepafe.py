import argparse
import os
import sys
import pandas as pd
import glob
import BioSimSpace as bss
from BioSimSpace import _Exceptions
import functions as fn


parser = argparse.ArgumentParser(description="prepare AFE calculations")
parser.add_argument("system",
                     type=str,
                     help="system name; this is used to find the folder containing input files")
parser.add_argument("network",
                    type=str,
                    help="filename including full path to network file")

arguments = parser.parse_args()

system_name = arguments.system
network_file = arguments.network

full_path = os.getcwd() + "/"
if "scripts" in full_path:
    full_path = full_path.replace("/scripts/", "/")
elif system_name in full_path:
    full_path = full_path.replace(f"/{system_name}/", "/")
equilibration_path = full_path + system_name + "/equilibration/"

if not os.path.isfile(network_file):
    print(f"The input file {network_file} does not exist")
    sys.exit()

network = pd.read_csv(network_file, header=None, comment="#", names=["ligand_1", 
                                                                     "ligand_2", 
                                                                     "n_windows", 
                                                                     "windows",
                                                                     "engine"])

columns_to_list = lambda column: network[column].tolist()
first_ligands = columns_to_list("ligand_1")
second_ligands = columns_to_list("ligand_2")
n_windows = columns_to_list("n_windows")
windows = columns_to_list("windows")
engines = columns_to_list("engine")
n_transformations = len(first_ligands)

lambda_values_string = [lambdas.split() for lambdas in windows]
lambda_values = [[float(value) for value in lambda_list] for lambda_list in lambda_values_string]

afe_folder_path = full_path + system_name + "/afe/"
protocol_file = afe_folder_path + "protocol.dat"
with open(protocol_file, "r") as file:
    protocol = file.readlines()
runtime = protocol[6].rstrip().replace(" ","").split("=")[-1].split("*")[0]

input_runtime_unit = protocol[6].rstrip().replace(" ","").split("=")[-1].split("*")[1]
if input_runtime_unit not in ["ns", "ps"]:
    raise NameError("Input runtime unit not supported. Please use 'ns' or 'ps'" \
                   +" on the seventh line of protocol.dat in the shape of (e.g.):\nsampling = 2*ns")

if input_runtime_unit == "ns":
    runtime_unit = bss.Units.Time.nanosecond
elif input_runtime_unit == "ps":
    runtime_unit = bss.Units.Time.picosecond

try:
    runtime = int(runtime)
except ValueError:
    raise NameError("Input runtime value not supported. Please use an integer" \
                   +" on the seventh line of protocol.dat in the shape of (e.g.):\nsampling = 2*ns")

counter = 1
for i in range(n_transformations):

    get_ligand_number = lambda ligand_name: ligand_name.split("_")[-1]
    ligand_1_number = get_ligand_number(first_ligands[i])
    ligand_2_number = get_ligand_number(second_ligands[i])
    unbound_path = equilibration_path + "unbound/"
    bound_path = equilibration_path + "bound/"
    print(f"progress: {counter}/{n_transformations}")
    print(f"ligand 1: {ligand_1_number}, ligand 2: {ligand_2_number}")
    engine = engines[i].rstrip()

    ligand_1_system = bss.IO.readMolecules([f"{unbound_path}/ligand_{ligand_1_number}/npt/ligand_{ligand_1_number}.prm7",
                                            f"{unbound_path}/ligand_{ligand_1_number}/npt/ligand_{ligand_1_number}.rst7"])
    ligand_2_system = bss.IO.readMolecules([f"{unbound_path}/ligand_{ligand_2_number}/npt/ligand_{ligand_2_number}.prm7",
                                            f"{unbound_path}/ligand_{ligand_2_number}/npt/ligand_{ligand_2_number}.rst7"])

    system_1 = bss.IO.readMolecules([f"{bound_path}/system_{ligand_1_number}/npt/system_{ligand_1_number}.prm7",
                                        f"{bound_path}/system_{ligand_1_number}/npt/system_{ligand_1_number}.rst7"])
    system_2 = bss.IO.readMolecules([f"{bound_path}/system_{ligand_2_number}/npt/system_{ligand_2_number}.prm7",
                                        f"{bound_path}/system_{ligand_2_number}/npt/system_{ligand_2_number}.rst7"])

    ligand_1 = ligand_1_system.getMolecule(0)
    ligand_2 = ligand_2_system.getMolecule(0)
    mapping = bss.Align.matchAtoms(ligand_1, ligand_2, complete_rings_only=True) 
    inverse_mapping = {value:key for key,value in mapping.items()}
    ligand_2_aligned = bss.Align.rmsdAlign(ligand_2, ligand_1, inverse_mapping)

    merged_ligands = bss.Align.merge(ligand_1, ligand_2_aligned, mapping, allow_ring_breaking=True, allow_ring_size_change=True) 

    ligand_1_system.removeMolecules(ligand_1)
    ligand_1_system.addMolecules(merged_ligands)
    unbound_system = ligand_1_system

    system_1_ligand = None
    protein = None
    n_residues = [molecule.nResidues() for molecule in system_1]
    n_atoms = [molecule.nAtoms() for molecule in system_1]
    for j, (n_residues, n_atoms) in enumerate(zip(n_residues[:20], n_atoms[:20])):
        if n_residues == 1 and n_atoms > 5:
            system_1_ligand = system_1.getMolecule(j)
        elif n_residues > 1:
            protein = system_1.getMolecule(j)
        else:
            pass
    
    system_2_ligand = None
    n_residues = [molecule.nResidues() for molecule in system_2]
    n_atoms = [molecule.nAtoms() for molecule in system_2]
    for j, (n_residues, n_atoms) in enumerate(zip(n_residues, n_atoms)):
        if n_residues == 1 and n_atoms > 5:
            system_2_ligand = system_2.getMolecule(j)
        else:
            pass    

    if system_1_ligand and system_2_ligand and protein:
        print(f"got molecules")
    else:
        raise _Exceptions.AlignmentError("Could not extract ligands or protein from input systems.")
    
    mapping = bss.Align.matchAtoms(system_1_ligand, system_2_ligand, complete_rings_only=True)
    inverse_mapping = {value:key for key,value in mapping.items()}

    system_2_ligand_aligned = bss.Align.rmsdAlign(system_2_ligand, system_1_ligand, inverse_mapping)

    bound_merged_ligands = bss.Align.merge(system_1_ligand, system_2_ligand_aligned, mapping, allow_ring_breaking=True, allow_ring_size_change=True)

    system_1.removeMolecules(system_1_ligand)
    system_1.addMolecules(bound_merged_ligands)
    bound_system = system_1

    free_energy_protocol = bss.Protocol.FreeEnergy(lam_vals=lambda_values[i], runtime=runtime*runtime_unit)

    working_directory = f"{full_path}/{system_name}/outputs/{engines[i].strip()}/lig_{ligand_1_number}~lig_{ligand_2_number}"
    bound_directory = working_directory + "/bound/"
    unbound_directory = working_directory + "/unbound/"

    bss.FreeEnergy.Relative(bound_system, free_energy_protocol, engine=engines[i].strip(), work_dir=bound_directory, setup_only=True)
    bss.FreeEnergy.Relative(unbound_system, free_energy_protocol, engine=engines[i].strip(), work_dir=unbound_directory, setup_only=True)

    bound_minimisation_directory = fn.create_dirs(bound_directory + "minimisation/")
    unbound_minimisation_directory = fn.create_dirs(unbound_directory + "minimisation/")

    bss.FreeEnergy.Relative(bound_system, free_energy_protocol, engine=engines[i].strip(), work_dir=bound_minimisation_directory, setup_only=True)
    bss.FreeEnergy.Relative(unbound_system, free_energy_protocol, engine=engines[i].strip(), work_dir=unbound_minimisation_directory, setup_only=True)

    bound_configuration_files = glob.glob(bound_minimisation_directory + "/*/*.cfg")
    unbound_configuration_files = glob.glob(unbound_minimisation_directory + "/*/*.cfg")

    minimisation_config = ["minimise = True\n", "minimise maximum iterations = 10000\n"]

    for i in range(len(bound_configuration_files)):
        with open(bound_configuration_files[i], "r") as file:
            old_config = file.readlines()
        with open(bound_configuration_files[i], "w") as file:
            replaced_config = [setting.replace("ncycles = 5\n", "ncycles = 1\n").replace("nmoves = 200000", "nmoves = 50000") for setting in old_config]
            for min_setting in minimisation_config:
                replaced_config.append(min_setting)
            file.writelines(replaced_config)
        with open(unbound_configuration_files[i], "r") as file:
            old_config = file.readlines()
        with open(unbound_configuration_files[i], "w") as file:
            replaced_config = [setting.replace("ncycles = 5\n", "ncycles = 1\n").replace("nmoves = 200000", "nmoves = 50000") for setting in old_config]
            for min_setting in minimisation_config:
                replaced_config.append(min_setting)
            file.writelines(replaced_config)

    counter += 1

 
