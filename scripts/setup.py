import warnings
warnings.filterwarnings("ignore")
import BioSimSpace as bss
import argparse
import os
import sys
import glob
import pathlib
import csv
import numpy as np
import pandas as pd
import MDAnalysis as mda


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


parser = argparse.ArgumentParser(description="setup AFE calculations with tLEAP and LOMAP")

parser.add_argument("system",
                     type=str,
                     help="system name; this is used to find the folder containing input files")

parser.add_argument("protein",
                    type=str,
                    help="protein pdb file in ../system/inputs/protein/")

parser.add_argument("-a",
                    "--active-site",
                    type=str,
                    help="csv file containing active site residues for the zinc site, required for --forcefiled=zaff")

parser.add_argument("-w",
                    "--water",
                    type=str,
                    help="water pdb file in ../system/inputs/protein/")

parser.add_argument("-ff",
                    "--forcefield",
                    type=str,
                    help="forcefield for protein",
                    choices=["ff14SB", "zaff"],
                    default="ff14SB")

parser.add_argument("-s",
                    "--solvent",
                    type=str,
                    help="solvent",
                    default="tip3p")

parser.add_argument("-o",
                    "--output",
                    type=str,
                    help="output file",
                    default="")

parser.add_argument("-e",
                    "--engine",
                    type=str,
                    help="MD engine for AFE",
                    default="GROMACS")

parser.add_argument("-n", 
                    "--n-windows",
                    type=str,
                    help="number of lambda windows for regular transformations",
                    choices=["3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"],
                    default="11")

parser.add_argument("-d", 
                    "--difficult",
                    type=str,
                    help="number of lambda windows for difficult transformations",
                    choices=["4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"],
                    default="17")

parser.add_argument("-t",
                    "--threshold",
                    type=str,
                    help="lomap score threshold for defining difficult transformations", 
                    choices=["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"],
                    default="0.4")

arguments = parser.parse_args()
system = arguments.system

full_path = os.getcwd() + "/"
if "scripts" in full_path:
    full_path = full_path.replace("/scripts/", "/")

protein_path = full_path + system + "/inputs/protein/"

ligand_path = full_path + system + "/inputs/ligands/"
protein_file = arguments.protein
is_file(protein_path+protein_file)
water_file = arguments.water
output_file = arguments.output
forcefield = arguments.forcefield
solvent = arguments.solvent.lower()
tleap_file = protein_path + "tleap.in"
engine = arguments.engine.upper()
windows = arguments.n_windows
difficult_windows = arguments.difficult
threshold = arguments.threshold
active_site_file = arguments.active_site
# SYSTEM
if arguments.water is not None:
    water_file = protein_path + arguments.water
    is_file(water_file)
else: 
    water_file = None

complex_file = create_complex(protein_path, protein_file, water_file, output_file)
is_file(protein_path+complex_file+".pdb")


if forcefield.lower() == "zaff" and active_site_file is None:
    parser.error("--forcefield=zaff requires an active site file --active-site={filename}.csv")

elif forcefield.lower() == "zaff" and active_site_file:
    try:
        if not active_site_file.endswith(".csv"):
            raise IOError("Active site file should be a .csv file.")
        
        active_site = pd.read_csv(active_site_file, header=None, names=["resname", "resid"])
        active_site_residues = active_site["resname"].tolist()
        active_site_residue_ids = active_site["resid"].tolist()
        n_active_site_residues = len(active_site_residue_ids)
        universe = mda.Universe(protein_path+complex_file+".pdb")
        for i in range(n_active_site_residues):
            atom_group = universe.select_atoms(f"resid {active_site_residue_ids[i]}")
            atom_group.residues.resnames = active_site_residues[i]
        fixed_residues = protein_path+complex_file+"_fixed_residues.pdb"
        universe.atoms.write(fixed_residues)

    except FileNotFoundError as error:
        print(f"{error}: Running tleap with zaff requires an active site csv file")
    except IOError as error:
        print(error)
    with open(fixed_residues, "r") as pdb_input:
        lines = pdb_input.readlines()
    n_zn = 0
    for line in lines:
        if "ZN" in line:
            n_zn += 1
    protein_terminus = [i for i in range(len(lines)) if "ZN" in lines[i]][0] - 1
    zn_line = protein_terminus+1
    zn_terminus = protein_terminus + n_zn + 1
    conect = [i for i in range(len(lines)) if "CONECT" in lines[i]][0] 
    fixed_file = protein_path + complex_file + "_fixed.pdb"
    with open(fixed_file, "w") as pdb_output:
        pdb_output.writelines(lines[:protein_terminus])
        pdb_output.write("TER\n")
        pdb_output.writelines(lines[zn_line:zn_terminus])
        pdb_output.write("TER\n")
        pdb_output.writelines(lines[zn_terminus:conect])
        pdb_output.write("END\n")
    
    try:
        tleap_command = f"tleap -s -f {tleap_file} > {protein_path}" + "tleap.out"
        os.system(tleap_command)
    except FileNotFoundError:
        print("tleap input file for zaff does not exist and is required with --forcefield=zaff")
        sys.exit()

elif forcefield.lower() == "ff14sb":
    forcefield = "ff14SB"
    with open(tleap_file, "w") as tleap_in:
        tleap_in.write(f"source leaprc.protein.{forcefield}\n")
        tleap_in.write(f"source leaprc.water.{solvent}\n")
        tleap_in.write(f"complex = loadpdb {protein_path}{complex_file}.pdb\n")
        tleap_in.write(f"saveamberparm complex {protein_path+system}_tleap.prm7 {protein_path+system}_tleap.rst7\n")
        tleap_in.write("quit")
    tleap_command = f"tleap -s -f {tleap_file} > {protein_path}" + "tleap.out"
    os.system(tleap_command)


# LIGANDS
ligand_files = sorted(glob.glob(f"{ligand_path}docked_*.sdf"))
ligands = [bss.IO.readMolecules(filepath)[0] for filepath in ligand_files]
ligand_names = [filepath.split("/")[-1].replace(".sdf","") for filepath in ligand_files]
transformations, lomap_scores = bss.Align.generateNetwork(ligands, plot_network=True, names=ligand_names, work_dir=ligand_path)

perturbation_network_dict = {}
transformations_named = [(ligand_names[transf[0]], ligand_names[transf[1]]) for transf in transformations]
for transformation, score in zip(transformations_named, lomap_scores):
    perturbation_network_dict[transformation] = score

with open(ligand_path+f"/lomap_{system}.csv", "w") as lomap_out:
    for key, value in perturbation_network_dict.items():
        lomap_out.write(f"{key}: {value}\n")

print("The LOMAP-generated perturbation network is:\n")
dataframe = pd.DataFrame.from_dict(perturbation_network_dict, "index")
perturbation_dataframe = dataframe.reset_index().rename(columns={"index":"perturbations", 0: "score"})
perturbation_dataframe.index.name = "index"
print(perturbation_dataframe)

adjust_network = False
print("\nDo you want to edit this network? ([y]es/[n]o/[q]uit)")

while True:
    edit = input("> ")
    if edit.lower().strip() == "yes" or edit.lower().strip() == "y":
        adjust_network = True
        break
    elif edit.lower().strip() == "no" or edit.lower().strip() == "n":
        adjust_network = False
        break
    elif edit.lower().strip() == "quit" or edit.lower().strip() == "q":
        print("Quitting.")
        sys.exit()
    elif edit == "":
        continue
    else: 
        print("Invalid option.")
        continue

first_index = perturbation_dataframe.first_valid_index()
last_index = perturbation_dataframe.last_valid_index()
options = ["\tdel: delete perturbations by index",
           "\tadd: add perturbations in format 'add ligand_1 ligand_2 score'",
           "\tedit: edit existing LOMAP scores in format 'edit index score'",
           "\ts: save and continue preparing",
           "\tq: quit"]
edited_dataframe = perturbation_dataframe.copy()
while adjust_network == True:
    print("Choose an option:")
    for option in options:
        print(option)
    option = input("> ").lower()
    check_option = option.replace(",", "").split()
    n_inputs = len(check_option)
    if 1 < n_inputs <= 2 and check_option[0] == "del":
        try:
            delete_index = int(check_option[1])
            try:
                edited_dataframe = edited_dataframe.drop([delete_index], axis=0, inplace=False)
                print("\nThe NEW network is:\n")
                print(edited_dataframe)
            except KeyError:
                print(f"Error: delete index should be between {first_index} and {last_index}")
            continue
        except ValueError:
            print("Delete index should be an integer")
            continue

    elif n_inputs > 2 and check_option[0] == "del":
        try:
            delete_indices = [int(check_option[i]) for i in range(1, n_inputs)]
            try:
                edited_dataframe = edited_dataframe.drop(delete_indices, axis=0, inplace=False)
                print("\nThe NEW network is:\n")
                print(edited_dataframe)
            except KeyError:
                print(f"Error: delete index should be between {first_index} and {last_index}")
            continue
        except ValueError:
            print("Error: Delete index should be an integer")
            continue

    elif check_option[0] == "edit":
        try: 
            if n_inputs != 3:
                raise ValueError
            index = int(check_option[1])
            score = np.float64(check_option[2])
            edited_dataframe.at[index, "score"] = score 
            print("\nThe NEW network is:\n")
            print(edited_dataframe)
        except ValueError:
            print("edit existing LOMAP scores in format 'edit index score'")
        continue
    elif check_option[0] == "add":            
        try:
            if n_inputs != 4:
                raise ValueError
            transformation = tuple([check_option[1], check_option[2]])
            score = np.float64(check_option[3])
            added_row = pd.DataFrame([[transformation, score]], columns= ["perturbations", "score"], index=[last_index+1])
            edited_dataframe = edited_dataframe.append(added_row, ignore_index=False)
            print("\nThe NEW network is:\n")
            print(edited_dataframe)
            last_index = edited_dataframe.last_valid_index()
        except ValueError:
            print("Error: add perturbations in format 'add ligand_1, ligand_2, score'")
        continue

    elif option.lower() == "s":
        edited_dataframe.to_csv(ligand_path+"adjusted_network.csv", sep=":", header=None, index=False)
        break
    elif n_inputs == 1 and option.lower() != "q":
        print("Error: expected index")
        continue
    elif option.lower() == "q":
        print("Quitting.")
        sys.exit()
    elif option == "":
        continue
    else:
        print("Invalid option.")
        continue


adjusted_dict = dict(zip(edited_dataframe["perturbations"], edited_dataframe["score"]))
production_directory = full_path + system + "/afe/"
pathlib.Path(production_directory).mkdir(parents=True, exist_ok=True)

with open(production_directory + "ligands.dat", "w") as ligands_file:
    writer = csv.writer(ligands_file)
    for ligand in ligand_names:
        writer.writerow([ligand])

with open(production_directory + "network_fwd.dat", "w") as network_file:
    for perturbation, lomap_score in adjusted_dict.items():
        if lomap_score == None or lomap_score < float(threshold):
            n_windows = difficult_windows
        else:
            n_windows = windows
        lambda_list_numpy = list(np.around(np.linspace(0, 1, int(n_windows)), decimals=5))
        lambda_list = [str(item) for item in lambda_list_numpy]

        lambda_array_bash = " ".join(lambda_list)
        network_file.write(f"{perturbation[0]}, {perturbation[1]}, {len(lambda_list_numpy)}, {lambda_array_bash}, {engine}\n")

with open(production_directory + "network_bwd.dat", "w") as network_file:
    for perturbation, lomap_score in adjusted_dict.items():
        if lomap_score == None or lomap_score < float(threshold):
            n_windows = difficult_windows
        else:
            n_windows = windows
        lambda_list_numpy = list(np.around(np.linspace(0, 1, int(n_windows)), decimals=5))
        lambda_list = [str(item) for item in lambda_list_numpy]

        lambda_array_bash = " ".join(lambda_list)
        network_file.write(f"{perturbation[1]}, {perturbation[0]}, {len(lambda_list_numpy)}, {lambda_array_bash}, {engine}\n")

protocol = [f"ligand forcefield = gaff2", # CHANGE SO THAT USER CAN CHANGE
            f"protein forcefield = {forcefield}", 
            f"solvent = {solvent}", 
            f"box edges = 20*angstrom", # CHANGE SO THAT USER CAN CHANGE
            f"box shape = orthorhombic", # CHANGE SO THAT USER CAN CHANGE
            f"protocol = default",
            f"sampling = 2*ns", # CHANGE SO THAT USER CAN CHANGE
            f"engine = {engine}"]

with open(production_directory + "protocol.dat", "w") as protocol_file:
    writer = csv.writer(protocol_file)
    for protocol_line in protocol:
        writer.writerow([protocol_line])
 
